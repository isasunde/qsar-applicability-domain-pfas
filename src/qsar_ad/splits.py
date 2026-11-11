"""Random, scaffold and cluster splits, with leakage checks."""
from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina


# ---------- Fingerprints and groups ----------
def morgan_fps(smiles, radius: int = 2, n_bits: int = 2048) -> list:
    """Build RDKit Morgan fingerprints (ECFP) for a list of SMILES strings."""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    return [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in smiles]


def murcko_scaffold(smiles: str) -> str:
    """Bemis-Murcko scaffold SMILES based on ring systems and linkers of a molecule. Acyclic molecules get empty scaffold ('') and share one group."""
    scaf = MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(smiles))
    return Chem.MolToSmiles(scaf)

def butina_clusters(fps: list, distance_cutoff: float = 0.6) -> np.ndarray:
    """Cluster molecules with Butina alhorithm based on Tanimoto distance (1 - similarity). 
    I.e. similar molecules get the same cluster label & held-out proteins are far from the training data.
    Returns cluster labels pr moelcule."""
    n = len(fps)
    parts = []

    # Compute the lower triangle of the distance matrix (1 - Tanimoto similarity)
    for i in range(1, n): 
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        parts.append(1.0 - np.asarray(sims))
    dists = np.concatenate(parts)

    # Cluster with Butina algorithm
    clusters = Butina.ClusterData(dists, n, distance_cutoff, isDistData=True)
    
    # Assign cluster labels to each molecule
    labels = np.empty(n, dtype=int)
    for k, members in enumerate(clusters):
        labels[list(members)] = k
    return labels

def compute_groups(df: pd.DataFrame, method: str, distance_cutoff: float = 0.6) -> np.ndarray:
    """Returns one group label pr molecule for random, scaffold or cluster. 
    Expensive for 'cluster', so compute once and reuse."""
    if method == "random":
        return np.arange(len(df))
    if method == "scaffold":
        return pd.factorize(df["smiles"].map(murcko_scaffold))[0]
    if method == "cluster":
        return butina_clusters(morgan_fps(df["smiles"]), distance_cutoff)
    raise ValueError(f"unknown method: {method}")

# ---------- Grouped splitting ----------
def hold_out_groups(groups: np.ndarray, idx: np.ndarray, frac: float,
                     rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Move whole groups out of idx until held-out set reaches ~frac of full dataset. Returns (kept, held_out)."""
    # Get the group labels and their sizes for the molecules in idx
    g = groups[idx]
    sizes = pd.Series(g).value_counts()
    order = sizes.index.to_numpy().copy()
    rng.shuffle(order)

    # Determine which groups to hold out until we reach the target fraction
    target = frac * len(idx)
    held, count = [], 0
    for grp in order:
        if count >= target:
            break
        if count + sizes[grp] > target * 1.1:   # skip groups that overshoot a lot
            continue
        held.append(grp)
        count += sizes[grp]
    mask = np.isin(g, held)
    return idx[~mask], idx[mask]

def split_from_groups(groups: np.ndarray, test_size: float = 0.2, cal_size: float = 0.1,
                      seed: int = 42) -> dict[str, np.ndarray]:
    """Train/calibration/test index arrays. Groups never straddle two sets.
    cal_size and test_size are fractions of the full dataset."""
    # Use a random number generator for reproducibility
    rng = np.random.default_rng(seed)

    # Hold out test set first, then calibration set from the remaining training data
    all_idx = np.arange(len(groups))
    trainval, test = hold_out_groups(groups, all_idx, test_size, rng)
    if cal_size > 0:
        train, cal = hold_out_groups(groups, trainval, cal_size / (1 - test_size), rng)
    else:
        train, cal = trainval, np.array([], dtype=int)
    split = {"train": np.sort(train), "cal": np.sort(cal), "test": np.sort(test)}
    check_no_leakage(groups, split)
    return split


# ---------- Checks and diagnostics ----------
def check_no_leakage(groups: np.ndarray, split: dict[str, np.ndarray]) -> None:
    """Check that no group straddles two sets and that all indices are accounted for."""
    # Create sets of indices for each split
    sets = {k: set(v.tolist()) for k, v in split.items()}

    # Check that no index is missing or duplicated
    total = sum(len(s) for s in sets.values())
    assert total == len(groups), "indices missing or duplicated"

    # Check that no group straddles two sets
    names = list(sets)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            assert not sets[a] & sets[b], f"index overlap: {a}/{b}"
            ga, gb = set(groups[split[a]]), set(groups[split[b]])
            assert not ga & gb, f"group leakage between {a} and {b}"

def max_sim_to_train(fps: list, train_idx: np.ndarray, query_idx: np.ndarray) -> np.ndarray:
    """Highest Tanimoto similarity of each query molecule to any training molecule."""
    train_fps = [fps[i] for i in train_idx]
    return np.array([max(DataStructs.BulkTanimotoSimilarity(fps[i], train_fps))
                     for i in query_idx])

def split_report(fps: list, split: dict[str, np.ndarray]) -> dict:
    """Diagnostics for a single split. Returns number of molecules in each set and similarity of test to train."""
    sim = max_sim_to_train(fps, split["train"], split["test"])
    return {"n_train": len(split["train"]), "n_cal": len(split["cal"]),
            "n_test": len(split["test"]),
            "median_max_sim_test_to_train": float(np.median(sim)),
            "frac_test_below_0.4": float((sim < 0.4).mean())}