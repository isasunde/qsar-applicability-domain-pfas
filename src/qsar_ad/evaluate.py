"""Metrics, benchmark loop and summaries. Saves predictions for the AD and uncertainty steps."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from tqdm import tqdm

from .models import combos, make_model
from .splits import compute_groups, split_from_groups


def regression_metrics(y_true, y_pred) -> dict:
    return {"rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
            "mae": float(mean_absolute_error(y_true, y_pred)),
            "r2": float(r2_score(y_true, y_pred))}


def predict_with_std(model, X) -> tuple[np.ndarray, np.ndarray]:
    """Returns prediction and, for random forests, the spread across trees (NaN otherwise)."""
    reg = model.steps[-1][1]
    if hasattr(reg, "estimators_") and hasattr(reg.estimators_[0], "tree_"):
        Xt = X
        for _, step in model.steps[:-1]:
            Xt = step.transform(Xt)
        per_tree = np.stack([t.predict(Xt) for t in reg.estimators_])
        return per_tree.mean(axis=0), per_tree.std(axis=0)
    pred = model.predict(X)
    return pred, np.full(len(pred), np.nan)

def run_benchmark(df: pd.DataFrame, features: dict[str, np.ndarray],
                  methods=("random", "scaffold", "cluster"), seeds=range(5),
                  test_size: float = 0.2, cal_size: float = 0.1):
    """Run the benchmark for all methods, seeds, models, and features. Returns metrics and predictions."""
    # Setup metrics, groups and target variable
    y = df["y"].to_numpy()
    groups_by_method = {m: compute_groups(df, m) for m in methods}  # once per method
    metric_rows, pred_frames = [], []

    # Run the benchmark for all methods, seeds, models, and features
    jobs = [(m, s, mod, feat) for m in methods for s in seeds for mod, feat in combos()]
    for method, seed, mod, feat in tqdm(jobs, desc="benchmark"):
        # Split the data
        split = split_from_groups(groups_by_method[method], test_size, cal_size, seed)
        X = features[feat]

        # Setup model pipeline and fit on training data
        model = make_model(mod, feat, seed)
        model.fit(X[split["train"]], y[split["train"]])

        # Evaluate on calibration and test sets, saving metrics and predictions
        for subset in ("cal", "test"):
            idx = split[subset]
            pred, std = predict_with_std(model, X[idx])
            metric_rows.append({"method": method, "seed": seed, "model": mod,
                                "feature": "none" if mod == "mean" else feat,
                                "subset": subset, **regression_metrics(y[idx], pred)})
            pred_frames.append(pd.DataFrame({
                "row": idx, "method": method, "seed": seed, "model": mod,
                "feature": "none" if mod == "mean" else feat, "subset": subset,
                "y_true": y[idx], "y_pred": pred, "y_std": std}))

    return pd.DataFrame(metric_rows), pd.concat(pred_frames, ignore_index=True)


def summarize(metrics: pd.DataFrame, subset: str = "test") -> pd.DataFrame:
    """Mean and SD over seeds, one row per method/model/feature."""
    m = metrics[metrics["subset"] == subset]
    out = (m.groupby(["method", "feature", "model"])[["rmse", "mae", "r2"]]
             .agg(["mean", "std"]))
    out.columns = [f"{a}_{b}" for a, b in out.columns]
    return out.reset_index()