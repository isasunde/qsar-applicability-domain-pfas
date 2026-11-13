from pathlib import Path

import pandas as pd

from qsar_ad import (load_pfas, morgan_matrix, rdkit_descriptors,
                    coverage_by_similarity, coverage_summary, make_intervals,
                    pfas_uncertainty, spread_vs_error, plot_coverage_by_similarity)

# Set pandas display options and define output directories
pd.set_option("display.width", 200, "display.max_columns", None)
out = Path("results/tables")

# Read predictions and AD scores from 02_train_benchmark & 03_domain_analysis.py results
preds = pd.read_parquet(out / "predictions.parquet")
ad = pd.read_parquet(out / "ad_scores.parquet")

# ---- Part A: Coverage on the benchmark splits ----
# Compute conformal intervals for each split and seed
ivs = [make_intervals(preds, alpha=a, kind=k)
       for a in (0.1, 0.2) for k in ("absolute", "normalized")]
iv = pd.concat(ivs, ignore_index=True)

# Compute summary statistics and coverage by similarity
summ = coverage_summary(iv)
cond = coverage_by_similarity(iv, ad)
corr, binned = spread_vs_error(preds)

# Save results to CSV files
for name, t in [("conformal_summary", summ), ("conformal_by_similarity", cond),
                ("spread_error_correlation", corr), ("spread_error_bins", binned)]:
    t.to_csv(out / f"{name}.csv", index=False)

# Print summary statistics for alpha=0.1 and correlation results
print(summ[summ["alpha"] == 0.1].round(3))
print(corr.round(3))

# Plot coverage by similarity for alpha=0.1 and save the figure
plot_coverage_by_similarity(cond, path=Path("results/figures/coverage_by_similarity.png"))

# ---- Part B: PFAS ----
# Load PFAS list
pfas_path = Path("data/raw/Chemical List PFASSTRUCT-2026-11-12.csv")

if pfas_path.exists():
    # Load PFAS & training data and compute descriptors
    df = pd.read_parquet("data/processed/lipo_curated.parquet").reset_index(drop=True)
    pfas, _ = load_pfas(pfas_path, smiles_col="SMILES")

    # Compute descriptors for PFAS or load to save time
    cache = Path("data/processed/pfas_descriptors.parquet")
    if cache.exists():
        pdesc = pd.read_parquet(cache)
    else:
        pdesc = rdkit_descriptors(pfas["smiles"]); pdesc.to_parquet(cache)
    assert len(pdesc) == len(pfas)

    tdesc = pd.read_parquet("data/processed/rdkit_descriptors.parquet")  # cached by 02, if saved

    #Compute morgan fingerprints for training and PFAS sets
    X_pool = {"morgan": morgan_matrix(df["smiles"]), "descriptors": tdesc.to_numpy()}
    X_pf = {"morgan": morgan_matrix(pfas["smiles"]), "descriptors": pdesc.to_numpy()}

    # Fit models on each split's training set, calibrate on its calibration set, and report prediction spread and interval widths for test molecules and for PFAS
    res = pfas_uncertainty(df, X_pool, X_pf, method="cluster", seeds=range(3))
    res.to_csv(out / "pfas_uncertainty.csv", index=False)
    print(res.groupby(["feature", "set"])[["std", "width_absolute", "width_normalized"]]
             .median().round(3))
    
else:
    print("No PFAS list - skipping part B")