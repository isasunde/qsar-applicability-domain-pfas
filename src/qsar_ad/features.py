"""Featurization: Morgan fingerprints, RDKit descriptors, and a train-only descriptor cleaner."""
from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sklearn.base import BaseEstimator, TransformerMixin

RDLogger.DisableLog("rdApp.*")


def morgan_matrix(smiles, radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    """Binary Morgan fingerprints as a float32 matrix (same bits as the Tanimoto analysis)."""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    return np.vstack([gen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in smiles]
                     ).astype(np.float32)


def rdkit_descriptors(smiles) -> pd.DataFrame:
    """All RDKit 2D descriptors per molecule. Computing them is per-molecule, so doing
    it on the full dataset does not leak anything. Failures become NaN."""
    names = [n for n, _ in Descriptors._descList]
    rows = []
    for s in smiles:
        try:
            d = Descriptors.CalcMolDescriptors(Chem.MolFromSmiles(s), missingVal=np.nan)
            rows.append([d.get(n, np.nan) for n in names])
        except Exception:
            rows.append([np.nan] * len(names))
    return pd.DataFrame(rows, columns=names, dtype="float64")
