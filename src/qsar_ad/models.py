"""Leakage-safe model pipelines. Hyperparameters are fixed and identical across splits."""
from __future__ import annotations

import numpy as np
from lightgbm import LGBMRegressor
from sklearn.pipeline import Pipeline
from sklearn.linear_model import RidgeCV
from sklearn.impute import SimpleImputer
from sklearn.dummy import DummyRegressor
from sklearn.preprocessing import RobustScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import FunctionTransformer
from sklearn.feature_selection import VarianceThreshold

MODELS = ("mean", "ridge", "rf", "lgbm")
FEATURES = ("morgan", "descriptors")


def make_model(name: str, feature: str, seed: int) -> Pipeline:
    """Build a scikit-learn pipeline for the given model (mean baseline, ridge, random forrest or LightGBM)
    and feature type (descriptor or morgan). 
    The seed is used for random forests and lightGBM."""
    if name == "mean":
        return Pipeline([("reg", DummyRegressor(strategy="mean"))])

    # Preprocessing: only descriptors need cleaning; only ridge needs scaling.
    steps = []
    if feature == "descriptors":
        steps.append(("inf", FunctionTransformer(lambda X: np.where(np.isfinite(X), X, np.nan))))
        steps.append(("impute", SimpleImputer(strategy="median")))
        steps.append(("var", VarianceThreshold(0.0)))
        if name == "ridge":
            steps.append(("scale", RobustScaler()))
            steps.append(("clip", FunctionTransformer(lambda X: np.clip(X, -5, 5))))

    if name == "ridge":
        # RidgeCV picks alpha by efficient leave-one-out CV inside the training set.
        reg = RidgeCV(alphas=np.logspace(-2, 4, 25))
    elif name == "rf":
        reg = RandomForestRegressor(
            n_estimators=500, max_features="sqrt" if feature == "morgan" else 0.33,
            min_samples_leaf=1, n_jobs=-1, random_state=seed)
    elif name == "lgbm":
        reg = LGBMRegressor(
            n_estimators=500, learning_rate=0.05, num_leaves=31, min_child_samples=20,
            subsample=0.8, subsample_freq=1, colsample_bytree=0.5,
            random_state=seed, n_jobs=-1, verbose=-1)
    else:
        raise ValueError(f"unknown model: {name}")
    steps.append(("reg", reg))
    return Pipeline(steps)


def combos() -> list[tuple[str, str]]:
    """List the (model, feature) pairs. The mean baseline ignores features, so it is run once."""
    out = [("mean", "morgan")]
    out += [(m, f) for f in FEATURES for m in ("ridge", "rf", "lgbm")]
    return out