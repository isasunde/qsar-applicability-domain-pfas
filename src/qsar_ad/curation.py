"""Curation of property dataset: Standardize structures, Resolve duplicates, and Log every step."""
from __future__ import annotations

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

_chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
_uncharger = rdMolStandardize.Uncharger()


def standardize_mol(smiles: str):
    """Cleanup -> Largest fragment (strips salts) -> Neutralize. Returns Mol or None."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.Cleanup(mol)
        mol = _chooser.choose(mol)
        mol = _uncharger.uncharge(mol)
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None

class CurationLog:
    def __init__(self):
        self.rows = []

    def add(self, step: str, n_in: int, n_out: int, note: str = ""):
        self.rows.append(
            {"step": step, "records_in": n_in, "records_out": n_out,
             "removed": n_in - n_out, "note": note}
        )

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)

    def to_markdown(self) -> str:
        return self.to_frame().to_markdown(index=False)  # needs `tabulate`

def curate_chembl_logd(
    data: pd.DataFrame,
    standard_type: str = "LogD7.4",
    max_conflict_sd: float = 0.5,
) -> tuple[pd.DataFrame, CurationLog]:
    log = CurationLog()
    n = len(data)
    log.add("Raw export", n, n, "ChEMBL activity export")

    # 1. Small molecules only
    out = data
    if "Type" in out.columns:
        out = out[out["Type"] == "Small molecule"]
    log.add("Type == Small molecule", n, len(out)); n = len(out)

    # 2. Ensure only the desired experimental endpoint is kept
    out = out[out["Standard Type"] == standard_type].copy()
    log.add(f"Standard Type == {standard_type}", n, len(out)); n = len(out)

    # 3. Convert endpoint to numeric values
    out["y"] = pd.to_numeric(out["Standard Value"], errors="coerce")
    out = out[out["y"].notna()]
    log.add("Numeric Standard Value", n, len(out)); n = len(out)

    # 4. Keep identifiers and the target (deliberately drop calculated properties)
    keep = ["ChEMBL ID", "Smiles", "Assay ChEMBL ID", "Document ChEMBL ID"]
    if "Inchi Key" in out.columns:
        keep.append("Inchi Key")
    out = out[keep + ["y"]].rename(columns={
        "ChEMBL ID": "chembl_id", "Smiles": "smiles_raw", "Inchi Key": "inchikey_raw"})

    # 5. Standardize structures
    mols = out["smiles_raw"].map(standardize_mol)
    out = out[mols.notna()].copy()
    mols = mols[mols.notna()]
    log.add("RDKit parse + Standardize", n, len(out),
            "Unparsable or Failed standardization"); n = len(out)
    out["smiles"] = mols.map(Chem.MolToSmiles)
    out["inchikey"] = mols.map(Chem.MolToInchiKey)

    # 6. Resolve duplicates by standardized InChIKey
    g = out.groupby("inchikey")["y"]
    stats = pd.DataFrame({"y_median": g.median(), "y_sd": g.std(), "n_meas": g.size()})
    conflicting = stats[(stats["n_meas"] > 1) & (stats["y_sd"] > max_conflict_sd)].index
    first = out.drop_duplicates("inchikey").set_index("inchikey")
    first = first.drop(columns="y").join(stats)
    log.add("Merge duplicate structures", n, len(first),
            f"{(stats['n_meas'] > 1).sum()} structures had repeats")
    n = len(first)

    # 7. Drop structures whose repeats disagree
    first = first.drop(index=conflicting)
    log.add(f"Drop conflicting repeats (SD > {max_conflict_sd})", n, len(first)); n = len(first)

    curated = first.reset_index().rename(columns={"y_median": "y"})
    curated = curated[["chembl_id", "smiles", "inchikey", "y", "n_meas",
                       "Assay ChEMBL ID", "Document ChEMBL ID"]]
    log.add("Final curated set", n, len(curated))
    return curated, log