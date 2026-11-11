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