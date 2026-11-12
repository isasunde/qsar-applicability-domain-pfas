
import numpy as np
import pandas as pd
from pathlib import Path
from qsar_ad import (compute_groups, split_from_groups, split_report, 
                     morgan_fps, morgan_matrix, rdkit_descriptors,
                     run_benchmark, summarize, plot_split_performance)

data = pd.read_parquet("data/processed/lipo_curated.parquet").reset_index(drop=True)
fps = morgan_fps(data["smiles"])

# Run splitting diagnostics for each method and seed, printing the first split report
for method in ["random", "scaffold", "cluster"]:
    groups = compute_groups(data, method)          # once per method
    sizes = np.bincount(groups)
    print(method, "groups:", len(sizes), "largest:", sizes.max(),
          "singletons:", (sizes == 1).sum())

    results = []
    for seed in range(5):
        split = split_from_groups(groups, seed=seed)
        results.append(split_report(fps, split))

    print(f"\nMean split report for {method} :")
    print(pd.DataFrame(results).mean().to_string())

# Run the benchmark for all methods, seeds, models, and features
desc = rdkit_descriptors(data["smiles"])
features = {"morgan": morgan_matrix(data["smiles"]),
            "descriptors": desc.to_numpy()}
metrics, preds = run_benchmark(data, features, seeds=range(5))

summary = summarize(metrics)
print(summary)

# Plot the summary of test set performance for each method/model/feature
plot_split_performance(summary, metric="rmse", path="results/figures/split_performance.png")

# Store metrics, predictions, and summary to disk
Path("results/tables").mkdir(parents=True, exist_ok=True)
metrics.to_csv("results/tables/metrics_by_seed.csv", index=False)
preds.to_parquet("results/tables/predictions.parquet", index=False)
summary.to_csv("results/tables/summary_test.csv", index=False)


                
                

        
