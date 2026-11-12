from pathlib import Path

import pandas as pd

from qsar_ad.domain import (ad_inside_outside, attach_errors, compute_split_ad, coverage_curve,
                            error_by_similarity, error_similarity_correlation, fluorine_stats,
                            load_pfas, reference_similarity)
from qsar_ad.features import morgan_matrix
from qsar_ad.plotting import (plot_coverage_curve, plot_error_vs_similarity,
                              plot_similarity_distributions)

pd.set_option("display.width", 200, "display.max_columns", None)
out = Path("results/tables"); fig_dir = Path("results/figures")
fig_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_parquet("data/processed/lipo_curated.parquet").reset_index(drop=True)
X = morgan_matrix(df["smiles"])
preds = pd.read_parquet(out / "predictions.parquet")

# ---- Part A: does the AD predict error? ----
ad = compute_split_ad(df, X)                      # same defaults as run_benchmark
d = attach_errors(preds, ad)
ad.to_parquet(out / "ad_scores.parquet")

binned = error_by_similarity(d)
corr = error_similarity_correlation(d)
inout = ad_inside_outside(d, "in_ad")
inout_fixed = ad_inside_outside(d, "near_train")
curve = coverage_curve(d)
for name, t in [("error_by_similarity", binned), ("error_similarity_correlation", corr),
                ("ad_inside_outside", inout), ("ad_inside_outside_sim04", inout_fixed),
                ("coverage_curve", curve)]:
    t.to_csv(out / f"{name}.csv", index=False)

print(corr.round(3)); print(inout.round(3))
plot_error_vs_similarity(binned, path=fig_dir / "error_vs_similarity.png")
plot_coverage_curve(curve, path=fig_dir / "coverage_curve.png")

# ---- Fluorine content of the training set ----
fs = fluorine_stats(df["smiles"])
print("curated set: any F:", int((fs.n_F > 0).sum()),
      "| OECD-style CF2/CF3:", int(fs.pfas_oecd.sum()), "of", len(df))

# ---- Part B: PFAS coverage  ----
pfas_path = Path("data/raw/Chemical List PFASSTRUCT-2026-11-12.csv")
if pfas_path.exists():
    pfas, counts = load_pfas(pfas_path, smiles_col="SMILES")   # adjust column name
    print(counts)
    res = reference_similarity(pfas["smiles"], morgan_matrix(pfas["smiles"]), X)
    res.insert(0, "inchikey", pfas["inchikey"])
    res.to_csv(out / "pfas_coverage.csv", index=False)

    print(res[["max_sim", "knn_sim"]].describe().round(3))
    print(res.groupby("cf_group", observed=True).agg(
        n=("max_sim", "size"), median_max_sim=("max_sim", "median"),
        frac_in_ad=("in_ad", "mean"), frac_near=("near_train", "mean")).round(3))

    ref = {"test, random split": d[(d.method == "random") & (d.subset == "test")].drop_duplicates(["seed", "row"]).max_sim,
           "test, cluster split": d[(d.method == "cluster") & (d.subset == "test")].drop_duplicates(["seed", "row"]).max_sim,
           "PFAS list": res["max_sim"]}
    plot_similarity_distributions(ref, path=fig_dir / "pfas_similarity.png")
else:
    print("No PFAS list at", pfas_path, "- skipping part B")