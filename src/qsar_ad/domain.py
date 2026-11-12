"""Applicability-domain (AD) analysis: similarity to the training set, error versus
similarity, coverage-error curves, and PFAS coverage. Uses the predictions saved by
run_benchmark, so nothing is retrained."""
from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.stats import spearmanr

from .curation import standardize_pfas
from .splits import compute_groups, split_from_groups

KEYS = ["method", "feature", "model", "seed"]
SIM_BINS = [0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0001]


# ---------- similarity ----------

def tanimoto_matrix(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Tanimoto similarity between rows of two binary matrices."""
    inter = A @ B.T
    union = A.sum(1)[:, None] + B.sum(1)[None, :] - inter
    return inter / np.maximum(union, 1.0)


def ad_scores(Q: np.ndarray, R: np.ndarray, k: int = 5, exclude_self: bool = False,
              chunk: int = 1000) -> tuple[np.ndarray, np.ndarray]:
    """For each query row: max similarity to the reference set, and the mean similarity
    of its k nearest reference molecules. exclude_self is for Q is R (leave-one-out)."""
    max_sim = np.empty(len(Q))
    knn_sim = np.empty(len(Q))
    for start in range(0, len(Q), chunk):
        S = tanimoto_matrix(Q[start:start + chunk], R)
        if exclude_self:
            r = np.arange(S.shape[0])
            S[r, start + r] = -1.0
        top = -np.partition(-S, k - 1, axis=1)[:, :k]
        max_sim[start:start + len(S)] = top.max(axis=1)
        knn_sim[start:start + len(S)] = top.mean(axis=1)
    return max_sim, knn_sim


def ad_threshold(train_X: np.ndarray, k: int = 5, quantile: float = 0.05) -> float:
    """Data-derived AD threshold from the TRAINING set only: the 5th percentile of
    leave-one-out kNN similarity. A molecule with a lower kNN similarity than 95% of
    training molecules have to the rest of the training set is 'outside'."""
    _, knn = ad_scores(train_X, train_X, k, exclude_self=True)
    return float(np.quantile(knn, quantile))


# ---------- AD scores for every split ----------

def compute_split_ad(df: pd.DataFrame, X: np.ndarray, methods=("random", "scaffold", "cluster"),
                     seeds=range(5), test_size: float = 0.2, cal_size: float = 0.1,
                     k: int = 5, quantile: float = 0.05) -> pd.DataFrame:
    """Recreates each split (same groups, same seed, so identical to the benchmark) and
    scores calibration and test molecules against that split's training set.
    test_size and cal_size MUST match run_benchmark."""
    frames = []
    for method in methods:
        groups = compute_groups(df, method)
        for seed in seeds:
            split = split_from_groups(groups, test_size, cal_size, seed)
            Xtr = X[split["train"]]
            thr = ad_threshold(Xtr, k, quantile)
            for subset in ("cal", "test"):
                idx = split[subset]
                mx, kn = ad_scores(X[idx], Xtr, k)
                frames.append(pd.DataFrame({
                    "method": method, "seed": seed, "subset": subset, "row": idx,
                    "max_sim": mx, "knn_sim": kn, "ad_threshold": thr,
                    "in_ad": kn >= thr, "near_train": mx >= 0.4}))
    return pd.concat(frames, ignore_index=True)


def attach_errors(preds: pd.DataFrame, ad: pd.DataFrame) -> pd.DataFrame:
    """Join AD scores onto saved predictions. Fails loudly if the splits don't match."""
    keys = ["method", "seed", "subset", "row"]
    m = preds.merge(ad, on=keys, how="left", validate="many_to_one")
    if m["max_sim"].isna().any():
        raise ValueError("Some predictions have no AD scores: split settings differ "
                         "from the benchmark run.")
    m["abs_err"] = (m["y_true"] - m["y_pred"]).abs()
    m["sq_err"] = m["abs_err"] ** 2
    return m


# ---------- error versus domain ----------

def error_by_similarity(d: pd.DataFrame, col: str = "max_sim", bins=SIM_BINS,
                        subset: str = "test") -> pd.DataFrame:
    """RMSE and MAE in similarity bins, pooled over seeds (n counts molecule-seed pairs)."""
    d = d[d["subset"] == subset].copy()
    d["sim_bin"] = pd.cut(d[col], bins, right=False)
    g = d.groupby(["method", "feature", "model", "sim_bin"], observed=True)
    return g.agg(n=("abs_err", "size"), mae=("abs_err", "mean"),
                 rmse=("sq_err", lambda s: float(np.sqrt(s.mean())))).reset_index()


def error_similarity_correlation(d: pd.DataFrame, col: str = "max_sim",
                                 subset: str = "test") -> pd.DataFrame:
    """Spearman correlation between similarity and absolute error (negative = AD works),
    per seed, then mean and SD over seeds."""
    d = d[d["subset"] == subset]
    r = (d.groupby(KEYS).apply(lambda g: spearmanr(g[col], g["abs_err"]).statistic)
           .rename("spearman").reset_index())
    return (r.groupby(["method", "feature", "model"])["spearman"]
              .agg(["mean", "std"]).reset_index())


def ad_inside_outside(d: pd.DataFrame, flag: str = "in_ad", subset: str = "test") -> pd.DataFrame:
    """RMSE inside vs outside the AD and the fraction outside, averaged over seeds.
    flag is 'in_ad' (data-derived threshold) or 'near_train' (max similarity >= 0.4)."""
    d = d[d["subset"] == subset].assign(outside=lambda x: ~x[flag])
    rm = (d.groupby(KEYS + ["outside"])["sq_err"].mean().pow(0.5)
            .unstack("outside").reindex(columns=[False, True]))
    rm.columns = ["rmse_inside", "rmse_outside"]
    rm["frac_outside"] = d.groupby(KEYS)["outside"].mean()
    return (rm.reset_index().groupby(["method", "feature", "model"])
              [["frac_outside", "rmse_inside", "rmse_outside"]].mean().reset_index())


def coverage_curve(d: pd.DataFrame, col: str = "max_sim", subset: str = "test",
                   grid=np.linspace(0.1, 1.0, 10)) -> pd.DataFrame:
    """RMSE when only the most-similar fraction of test molecules is kept ('coverage').
    If the AD is informative, RMSE falls as coverage falls."""
    d = d[d["subset"] == subset]
    rows = []
    for keys, g in d.groupby(KEYS):
        sq = g["sq_err"].to_numpy()[np.argsort(-g[col].to_numpy())]
        cum = np.sqrt(np.cumsum(sq) / np.arange(1, len(sq) + 1))
        for f in grid:
            rows.append((*keys, float(f), cum[max(int(round(f * len(sq))) - 1, 0)]))
    out = pd.DataFrame(rows, columns=KEYS + ["coverage", "rmse"])
    return (out.groupby(["method", "feature", "model", "coverage"])["rmse"]
               .agg(["mean", "std"]).reset_index())


# ---------- fluorine / PFAS ----------

# OECD-style definition (approximate): at least one fully fluorinated CF3 or CF2 carbon
# with no H, Cl, Br or I attached.
_CF3 = Chem.MolFromSmarts("[CX4;H0](F)(F)F")
_CF2 = Chem.MolFromSmarts("[CX4;H0](F)(F)([!#9;!#17;!#35;!#53])[!#9;!#17;!#35;!#53]")


def fluorine_stats(smiles) -> pd.DataFrame:
    rows = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        carbons = ({x[0] for x in m.GetSubstructMatches(_CF3)}
                   | {x[0] for x in m.GetSubstructMatches(_CF2)})
        rows.append({"n_F": sum(a.GetAtomicNum() == 9 for a in m.GetAtoms()),
                     "n_cf_carbons": len(carbons)})
    out = pd.DataFrame(rows)
    out["pfas_oecd"] = out["n_cf_carbons"] > 0
    return out


def load_pfas(path: str, smiles_col: str = "SMILES") -> tuple[pd.DataFrame, dict]:
    """Read a PFAS structure list, standardize it exactly like the training data, and
    report what was lost. Returns (table, counts)."""
    raw = pd.read_csv(path, encoding="utf-8-sig")     # explicit comma, strips the BOM
    raw["dtxsid"] = raw["DTXSID"].str.split("/").str[-1]
    counts = {"rows_in": len(raw)}
    smi = raw[smiles_col].tolist()
    mols = [standardize_pfas(s) if isinstance(s, str) and s.strip() else None for s in smi]
    ok = [m is not None for m in mols]
    counts["no_or_failed_structure"] = len(raw) - sum(ok)
    out = raw[ok].copy().reset_index(drop=True)
    mols = [m for m in mols if m is not None]
    out["smiles"] = [Chem.MolToSmiles(m) for m in mols]
    out["inchikey"] = [Chem.MolToInchiKey(m) for m in mols]
    out = out.drop_duplicates("inchikey").reset_index(drop=True)
    counts["after_dedup"] = len(out)
    return out, counts


def reference_similarity(smiles, X_query: np.ndarray, X_ref: np.ndarray,
                         k: int = 5, quantile: float = 0.05) -> pd.DataFrame:
    """Similarity of external molecules (e.g. PFAS) to a reference set, with AD flags
    using the threshold derived from that reference set."""
    mx, kn = ad_scores(X_query, X_ref, k)
    thr = ad_threshold(X_ref, k, quantile)
    out = fluorine_stats(smiles)
    out["max_sim"], out["knn_sim"] = mx, kn
    out["in_ad"], out["near_train"] = kn >= thr, mx >= 0.4
    out["cf_group"] = pd.cut(out["n_cf_carbons"], [-1, 0, 1, 3, 10_000],
                             labels=["0", "1", "2-3", "4+"])
    return out