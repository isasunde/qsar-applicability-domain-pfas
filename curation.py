from pathlib import Path
from qsar_ad import build_raw_table, curate_chembl_logd

raw = build_raw_table("data/raw/CHEMBL_30_10_2026.csv",
                      "data/raw/CHEMBL_01_11_2026.csv")
print(raw["Standard Type"].value_counts())
print(raw["Assay ChEMBL ID"].nunique(), raw["Document ChEMBL ID"].nunique())

curated, log = curate_chembl_logd(raw)
Path("data/processed").mkdir(parents=True, exist_ok=True)
Path("docs").mkdir(exist_ok=True)
curated.to_parquet("data/processed/lipo_curated.parquet", index=False)
Path("docs/curation_log.md").write_text(log.to_markdown())
print(log.to_frame())