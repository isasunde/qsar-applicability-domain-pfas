"""Split conformal prediction intervals from saved calibration/test predictions,
coverage diagnostics, and an uncertainty-vs-error check."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from .models import make_model
from .evaluate import predict_with_std
from .splits import compute_groups, split_from_groups

KEYS = ["method", "feature", "model", "seed"]
SIM_BINS = [0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0001]


def conformal_quantile(scores: np.ndarray, alpha: float) -> float:
    """Finite-sample split-conformal quantile: the ceil((n+1)(1-alpha))/n empirical
    quantile of the calibration scores."""
    n = len(scores)
    level = min(np.ceil((n + 1) * (1 - alpha)) / n, 1.0)
    return float(np.quantile(scores, level, method="higher"))


def make_intervals(preds: pd.DataFrame, alpha: float = 0.1, kind: str = "absolute",
                   eps_frac: float = 0.1) -> pd.DataFrame:
    """Calibrate on each fit's calibration rows and build intervals for its test rows.

    kind='absolute':   pred +/- q, with q the quantile of |residual| (constant width).
    kind='normalized': pred +/- q*(std+eps), with q the quantile of |residual|/(std+eps).
                       Needs the tree spread, so it exists only for random forests.
    eps is a fraction of the median calibration spread, which stops tiny spreads
    from producing near-zero widths.
    """
    out = []
    for keys, g in preds.groupby(KEYS):
        cal, test = g[g["subset"] == "cal"], g[g["subset"] == "test"]
        resid = (cal["y_true"] - cal["y_pred"]).abs().to_numpy()
        if kind == "absolute":
            q = conformal_quantile(resid, alpha)
            half = np.full(len(test), q)
        else:
            if cal["y_std"].isna().any():
                continue
            eps = eps_frac * cal["y_std"].median()
            q = conformal_quantile(resid / (cal["y_std"].to_numpy() + eps), alpha)
            half = q * (test["y_std"].to_numpy() + eps)
        t = test[["row", "y_true", "y_pred", "y_std"]].copy()
        t[["method", "feature", "model", "seed"]] = keys
        t["alpha"], t["kind"], t["half_width"] = alpha, kind, half
        t["covered"] = (t["y_true"] - t["y_pred"]).abs() <= half
        t["abs_err"] = (t["y_true"] - t["y_pred"]).abs()
        out.append(t)
    return pd.concat(out, ignore_index=True)


def coverage_summary(iv: pd.DataFrame) -> pd.DataFrame:
    """Empirical coverage and mean interval width per fit, then mean and SD over seeds."""
    per_seed = (iv.groupby(["kind", "alpha"] + KEYS)
                  .agg(coverage=("covered", "mean"), width=("half_width", lambda s: 2 * s.mean()))
                  .reset_index())
    return (per_seed.groupby(["kind", "alpha", "method", "feature", "model"])
                    [["coverage", "width"]].agg(["mean", "std"])
                    .pipe(lambda t: t.set_axis([f"{a}_{b}" for a, b in t.columns], axis=1))
                    .reset_index())


def coverage_by_similarity(iv: pd.DataFrame, ad: pd.DataFrame, col: str = "max_sim",
                           bins=SIM_BINS) -> pd.DataFrame:
    """Empirical coverage in bins of similarity to the training set, pooled over seeds."""
    ad_t = ad[ad["subset"] == "test"][["method", "seed", "row", col]]
    m = iv.merge(ad_t, on=["method", "seed", "row"], how="left", validate="many_to_one")
    if m[col].isna().any():
        raise ValueError("Some intervals have no AD scores: split settings differ.")
    m["sim_bin"] = pd.cut(m[col], bins, right=False)
    return (m.groupby(["kind", "alpha", "method", "feature", "model", "sim_bin"], observed=True)
              .agg(n=("covered", "size"), coverage=("covered", "mean"),
                   mean_width=("half_width", lambda s: 2 * s.mean()))
              .reset_index())


def spread_vs_error(preds: pd.DataFrame, n_bins: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Random forests only. (1) Spearman correlation between tree spread and absolute
    error per fit; (2) mean absolute error by spread quantile within each fit."""
    t = preds[(preds["subset"] == "test") & preds["y_std"].notna()].copy()
    t["abs_err"] = (t["y_true"] - t["y_pred"]).abs()
    corr = (t.groupby(KEYS).apply(lambda g: spearmanr(g["y_std"], g["abs_err"]).statistic)
              .rename("spearman").reset_index()
              .groupby(["method", "feature", "model"])["spearman"].agg(["mean", "std"])
              .reset_index())
    t["std_bin"] = t.groupby(KEYS)["y_std"].transform(
        lambda s: pd.qcut(s, n_bins, labels=False, duplicates="drop"))
    binned = (t.groupby(["method", "feature", "model", "std_bin"])
                .agg(n=("abs_err", "size"), mean_std=("y_std", "mean"),
                     mean_abs_err=("abs_err", "mean")).reset_index())
    return corr, binned

def pfas_uncertainty(df: pd.DataFrame, X_train_pool: dict, X_pfas: dict, method: str = "cluster",
                     seeds=range(3), alpha: float = 0.1, eps_frac: float = 0.1,
                     features=("morgan", "descriptors")) -> pd.DataFrame:
    """Fit a random forest on each split's training set, calibrate on its calibration
    set, and report prediction spread and interval widths for test molecules and for
    PFAS. X_train_pool / X_pfas: {feature name: matrix}; PFAS rows keep their order."""
    y = df["y"].to_numpy()
    groups = compute_groups(df, method)
    rows = []
    for seed in seeds:
        split = split_from_groups(groups, seed=seed)
        for feat in features:
            X, Xp = X_train_pool[feat], X_pfas[feat]
            model = make_model("rf", feat, seed).fit(X[split["train"]], y[split["train"]])
            pc, sc = predict_with_std(model, X[split["cal"]])
            pt, st = predict_with_std(model, X[split["test"]])
            pp, sp = predict_with_std(model, Xp)
            eps = eps_frac * np.median(sc)
            q_abs = conformal_quantile(np.abs(y[split["cal"]] - pc), alpha)
            q_norm = conformal_quantile(np.abs(y[split["cal"]] - pc) / (sc + eps), alpha)
            for name, std in (("test", st), ("pfas", sp)):
                rows.append(pd.DataFrame({
                    "set": name, "feature": feat, "seed": seed, "std": std,
                    "width_absolute": 2 * q_abs,
                    "width_normalized": 2 * q_norm * (std + eps)}))
    return pd.concat(rows, ignore_index=True)