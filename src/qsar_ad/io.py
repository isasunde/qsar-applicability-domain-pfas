"""Load and merge compound-level and activity-level ChEMBL exports."""
from __future__ import annotations
import pandas as pd


import csv
import pandas as pd


def read_chembl_csv(path: str) -> pd.DataFrame:
    """Read a ChEMBL export where every field is wrapped as "a";"b";"c".

    Splitting on the exact sequence ";" (quote-semicolon-quote) 
    to avoid splitting on semicolons inside quoted fields.
    """
    df = pd.read_csv(path, sep='";"', engine="python", quoting=csv.QUOTE_NONE)

    # The first column header/values start with ", the last end with ".
    df.columns = df.columns.str.strip('"').str.strip()
    for col in (df.columns[0], df.columns[-1]):
        df[col] = df[col].astype("string").str.strip('"')
    return df


def build_raw_table(molecule_path: str, activity_path: str) -> pd.DataFrame:
    """Merge ChEMBL molecule and activity exports on ChEMBL ID.
    The molecule file has one row per compound, the activity file has one row per assay result"""
    mol = read_chembl_csv(molecule_path).rename(
        columns={"Compound ChEMBL ID": "ChEMBL ID"})
    act = read_chembl_csv(activity_path).rename(
        columns={"Molecule ChEMBL ID": "ChEMBL ID"})

    if not mol["ChEMBL ID"].is_unique:
        raise ValueError("Molecule file has repeated ChEMBL IDs")

    # From the molecule file take only columns the activity file lacks
    # to avoid duplicated columns.
    mol_cols = ["ChEMBL ID"] + [c for c in mol.columns
                                if c not in act.columns and c != "ChEMBL ID"]
    data = act.merge(mol[mol_cols], on="ChEMBL ID", how="left",
                     validate="many_to_one", indicator=True)

    n_unmatched = (data["_merge"] != "both").sum()
    print(f"activity rows: {len(act)}, molecule rows: {len(mol)}, "
          f"unmatched after merge: {n_unmatched}")
    data = data.drop(columns="_merge").dropna(axis=1, how="all")
    return data