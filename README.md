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

  ## Data splits

Models are evaluated under three splitting strategies of increasing difficulty.
In each, the 4,190 curated compounds are divided into train (70%), calibration
(10%, used for conformal prediction) and test (20%) sets. Whole groups are
assigned to one set only, and the code asserts that no group is shared between
sets.

| Strategy | Grouping | Groups | Largest group | Singletons |
|---|---|---|---|---|
| Random | Each molecule is its own group | 4,190 | 1 | 4,190 |
| Scaffold | Bemis-Murcko scaffold | 2437 | 76 | 1,898 |
| Cluster | Butina clustering, Tanimoto distance cutoff 0.6 (Morgan radius 2, 2048 bits) | 1261 | 71 | 778 |

### How difficult is each split?

For each test molecule, I computed the highest Tanimoto similarity to any
training molecule (Morgan fingerprints, radius 2, 2048 bits) across 5 seeds.
The following values are the mean across all 5 seeds

| Split | Train / cal / test | Median max similarity (test to train) | Test molecules with max similarity < 0.4 |
|---|---|---|---|
| Random | 2,933 / 419 / 838 | 0.68 | 19% |
| Scaffold | 2,927 / 420 / 841 | 0.56 | 26% |
| Cluster | 2,928 / 421 / 840 | 0.37 | 65% |

Similarity to the training set falls steadily from random to scaffold to
cluster splits, so the three strategies probe progressively more novel
chemistry. Under the cluster split, two thirds of test molecules have no
training neighbour above 0.4 Tanimoto similarity, compared with about one fifth
under a random split.

### Notes
- Even the random split has 19% of test molecules below 0.4 similarity,
  because the dataset is structurally diverse (4,190 drug-like compounds,
  many without close analogues). Expect a smaller gap between random and
  scaffold performance than in more redundant QSAR datasets.
- Scaffold groups are uneven: acyclic molecules share one empty scaffold and
  are held out together, so the actual set sizes deviate slightly from the
  70/10/20 target.
- The cluster split is the closest proxy for PFAS, which are expected to sit
  far from drug-like training chemistry (this is tested in the applicability
  domain analysis).
- Splits are repeated over 5 random seeds (0-4) to give error bars.
- Calibration sets use the same grouping as the test set, so conformal
  intervals are calibrated under the same level of novelty as the test.
- Scalers and any tuning use the training set only.