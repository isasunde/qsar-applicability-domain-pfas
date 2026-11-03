# qsar-applicability-domain-pfas
Reproducible QSAR pipeline for lipophilicity prediction: RDKit curation, benchmarking under random/scaffold/cluster splits, applicability domain analysis, and conformal uncertainty, stress-tested on PFAS chemical space.

### Source
The dataset reproduces the AstraZeneca lipophilicity set used in MoleculeNet
(Lipophilicity), taken directly from ChEMBL rather than the MoleculeNet copy.

- **Document:** CHEMBL3301361, "Experimental in vitro DMPK and
  physicochemical data on a set of publicly disclosed compounds"
  (Wenlock and Tomkinson; source: AstraZeneca DMPK/physicochemical;
  DOI 10.6019/CHEMBL3301361). The document contains 15 assays (pKa,
  LogD7.4, solubility, plasma protein binding, clearance), of which only the
  LogD7.4 assay is used here.
- **Assay:** CHEMBL3301363, 4,200 measurements.
- **Endpoint:** experimental octanol/water distribution coefficient at
  pH 7.4 (logD7.4, unitless). This is logD, not logP, so ionizable compounds
  are measured at their pH 7.4 ionization state.

### Extraction
| File | Level | Query / filter | Extracted | Rows |
|---|---|---|---|---|
| `CHEMBL_30_10_2026.csv` | Compounds | `_metadata.related_documents.all_chembl_ids:CHEMBL3301361` | 2026-10-30 | 5,794 |
| `CHEMBL_01_11_2026.csv` | Activities | `document_chembl_id:CHEMBL3301361` + `standard_type:logD7.4` | 2026-11-01 | 4,200 |

The compound file covers every compound in the document (all 15 assays), so
it has more structures than there are LogD7.4 measurements. The activity file
is the base of the merge; all 4,200 rows matched a structure.

### What was used
Only the measured value (`Standard Value`, type `LogD7.4`) and the structure
(SMILES) were used. ChEMBL-calculated properties in the exports (AlogP,
polar surface area, QED and similar) were **not** used as features or
targets, to avoid circularity.

### Curation
Curation removed 10 of 4,200 records (0.24%), leaving **4,190** compounds.
The steps are documented in [`docs/curation_log.md`](docs/curation_log.md):

| Step | Records removed |
|---|---|
| Non-small-molecule entries | 5 |
| Non-numeric values / structure parsing failures | 0 |
| Duplicate structures merged (median value taken) | 5 |
| Conflicting repeats (SD > 0.5 log units) | 0 |

### Limitations
- All data come from one assay on drug-like compounds, so the chemical space
  is narrow. Few fluorinated compounds are present, which is the point of the
  PFAS applicability-domain analysis in this repo.
- Structures are standardized to neutral forms (largest fragment, uncharged)
  before featurization, while the measurement reflects the pH 7.4 state.
- The raw files can't be reproduced exactly from the queries alone, because
  the ChEMBL web interface changes between releases. Use the files in
  `data/raw/` for exact reproduction.