import numpy as np
import matplotlib.pyplot as plt

def plot_split_performance(summary, metric="rmse", path=None):
    """Draw grouped bar chart of test performance for each split method, model, and feature."""
    # Create a label for each model/feature combination, preserving the order of first appearance
    s = summary.assign(label=lambda d: d["model"] + " / " + d["feature"])
    labels = list(dict.fromkeys(s["label"]))

    fig, ax = plt.subplots(figsize=(9, 4))
    w = 0.25

    # Draw a bar for each method, offsetting the x position by the method index
    for i, method in enumerate(["random", "scaffold", "cluster"]):
        d = s[s["method"] == method].set_index("label").loc[labels]
        ax.bar([j + (i - 1) * w for j in range(len(labels))], d[f"{metric}_mean"],
               w, yerr=d[f"{metric}_std"], capsize=2, label=method)
    
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel(f"test {metric.upper()} (mean ± SD, 5 seeds)")
    ax.legend(title="split")
    fig.tight_layout()

    # Option to save the figure
    if path:
        fig.savefig(path, dpi=200)
    return fig

def plot_error_vs_similarity(binned, model="lgbm", feature="descriptors", path=None):
    d = binned[(binned["model"] == model) & (binned["feature"] == feature)]
    fig, ax = plt.subplots(figsize=(6, 4))
    for method, g in d.groupby("method"):
        ax.plot([i.left for i in g["sim_bin"]], g["rmse"], marker="o", label=method)
    ax.set_xlabel("max Tanimoto similarity to training set (bin start)")
    ax.set_ylabel("test RMSE")
    ax.set_title(f"{model} / {feature}")
    ax.legend(title="split"); fig.tight_layout()
    if path: fig.savefig(path, dpi=200)
    return fig

def plot_coverage_curve(curve, model="lgbm", feature="descriptors", path=None):
    d = curve[(curve["model"] == model) & (curve["feature"] == feature)]
    fig, ax = plt.subplots(figsize=(6, 4))
    for method, g in d.groupby("method"):
        ax.errorbar(g["coverage"], g["mean"], yerr=g["std"], marker="o",
                    capsize=2, label=method)
    ax.set_xlabel("coverage (fraction of test molecules kept, most similar first)")
    ax.set_ylabel("RMSE"); ax.legend(title="split"); fig.tight_layout()
    if path: fig.savefig(path, dpi=200)
    return fig

def plot_similarity_distributions(series: dict, threshold=0.4, path=None):
    """series: {label: array of max similarities}."""
    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.linspace(0, 1, 41)
    for label, v in series.items():
        ax.hist(v, bins=bins, density=True, histtype="step", linewidth=1.8, label=label)
    ax.axvline(threshold, color="grey", linestyle="--")
    ax.set_xlabel("max Tanimoto similarity to training set")
    ax.set_ylabel("density"); ax.legend(); fig.tight_layout()
    if path: fig.savefig(path, dpi=200)
    return fig