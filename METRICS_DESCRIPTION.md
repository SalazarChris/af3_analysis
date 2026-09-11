# AF3 Confidence Metrics — Description for Analysis

This document describes every metric extracted from AlphaFold 3 (AF3)
output files by the analysis pipeline. It is written so that a
statistician or LLM can understand what each metric represents, what
values it can take, and how to interpret differences between experimental
conditions.

The document is experiment-agnostic. No biological assumptions are made
about what the chains represent, how many conditions exist, or what
combinations of experimental attributes are tested.

---

## 1. Experimental Design

**Prediction tool:** AlphaFold 3 (server or local)

**Replication:** Each experimental condition is predicted with multiple
independent seeds, each producing multiple stochastic samples. The
number of seeds and samples per seed is determined by the experimental
design (specified in `experiment_metadata.json`).

**Conditions:** The number of conditions and their defining attributes
vary by experiment. Conditions are defined in the experiment metadata
file, not inferred from file names.

**Unit of analysis:** One prediction (= one seed × one sample).

---

## 2. Raw AF3 Output Files

Each AF3 prediction produces two JSON files:

### 2a. `*_confidences.json` (per-atom / per-token confidence)

Contains the fine-grained confidence data extracted from AF3:

| Field | Type | Description |
|-------|------|-------------|
| `atom_plddts` | list[float] | Per-atom pLDDT score (one per atom in the structure) |
| `atom_chain_ids` | list[str] | Chain identifier for each atom (e.g. "A", "B", "C") |
| `token_chain_ids` | list[str] | Chain identifier for each token (= one per residue or nucleotide) |
| `token_res_ids` | list[int] | Residue/nucleotide index within each chain |
| `contact_probs` | list[float] | Predicted probability that two tokens are in contact |
| `pae` | list[float] | Predicted Aligned Error between pairs of tokens |

### 2b. `*_summary_confidences.json` (per-complex summary scores)

Contains scalar and per-chain summary metrics:

| Field | Type | Description |
|-------|------|-------------|
| `ptm` | float | predicted TM-score for the whole complex |
| `iptm` | float | interface predicted TM-score (measures inter-chain accuracy) |
| `ranking_score` | float | AF3's composite ranking score (higher = better prediction) |
| `fraction_disordered` | float | Fraction of residues predicted as disordered |
| `has_clash` | float | 1.0 if atomic clashes detected, 0.0 otherwise |
| `chain_ptm` | list[float] | Per-chain pTM scores |
| `chain_iptm` | list[float] | Per-chain ipTM scores |
| `chain_pair_iptm` | list[float] | Pairwise ipTM between all chain pairs (symmetric matrix) |
| `chain_pair_pae_min` | list[float] | Minimum PAE between all chain pairs (symmetric matrix) |

### 2c. `*_ranking_scores.csv` (tabular ranking scores)

A CSV with columns: `seed`, `sample`, `ranking_score`.

---

## 3. Extracted Metrics (Pipeline Output)

The extraction pipeline (`af3_condition_centric_extraction.py`) reads
the raw JSONs and computes summary statistics. These are the metrics
available in the analysis CSVs.

### 3a. Global Metrics (per-prediction scalars)

These are computed from the full structure regardless of chain:

| Metric Name | Computed From | Scale | Interpretation |
|------------|---------------|-------|----------------|
| **pLDDT_mean** | mean of `atom_plddts` | 0–100 | Average per-atom confidence in local structure. >90 = very high confidence; 70–90 = confident; 50–70 = low; <50 = very low / probably disordered |
| **pLDDT_max** | max of `atom_plddts` | 0–100 | Highest-confidence atom. Ceiling of structural confidence |
| **pLDDT_min** | min of `atom_plddts` | 0–100 | Lowest-confidence atom. Floor of structural confidence |
| **pLDDT_median** | median of `atom_plddts` | 0–100 | Robust center of pLDDT distribution (less sensitive to outlier atoms than mean) |
| **contact_prob_mean** | mean of `contact_probs` | 0–1 | Average predicted probability that token pairs are in contact. Low values = sparse contacts; high values = dense contacts |
| **contact_prob_max** | max of `contact_probs` | 0–1 | Highest pairwise contact probability |
| **contact_prob_min** | min of `contact_probs` | 0–1 | Lowest pairwise contact probability |
| **contact_prob_median** | median of `contact_probs` | 0–1 | Robust center of contact probability distribution |
| **pae_mean** | mean of `pae` | 0–~50 Å | Average predicted alignment error. Lower = higher confidence in relative positioning of token pairs |
| **pae_max** | max of `pae` | 0–~50 Å | Worst-case positional error across all token pairs |
| **pae_min** | min of `pae` | 0–~50 Å | Best-case positional error (often ≈ 0 for self-pairs) |
| **pae_median** | median of `pae` | 0–~50 Å | Robust center of PAE distribution |
| **ranking_score** | from `*_ranking_scores.csv` or summary JSON | 0–1 | AF3's composite score used to rank models. Higher = better prediction. Combines pTM, ipTM, and other factors |

**Important notes on global metrics:**
- `pLDDT_mean` is dominated by the largest chain (by atom count). For
  multi-chain complexes, a high `pLDDT_mean` may mask low confidence in
  a small chain.
- `contact_prob_mean` scales with the number of token pairs; it is
  lower for larger complexes simply because there are more non-contact
  pairs.
- `pae_mean` is sensitive to inter-chain positioning. Low mean PAE
  across the full matrix suggests confident overall fold and assembly.
- `ranking_score` is the metric AF3 uses internally to select the
  "best" model. It is the most holistic single-number summary.

### 3b. Chain-Level Metrics (per-chain scalars)

Computed separately for each chain (molecule) in the complex:

| Metric Name | Computed From | Scale | Interpretation |
|------------|---------------|-------|----------------|
| **chain_{ID}_plddt** | mean of `atom_plddts` for atoms belonging to chain | 0–100 | Per-chain average pLDDT. Identifies which chains are well-predicted vs disordered or uncertain |
| **chain_{ID}_residues** | count of tokens in chain | integer | Number of residues/nucleotides in the chain. Constant across seeds within a condition (same sequence) |

The chain identifiers (A, B, C, ...) are assigned by AF3 and correspond
to the individual molecules in the complex. Which chain letter
corresponds to which biological entity depends on the specific
experiment and is not inferred by the pipeline.

**Important:** Chain IDs are assigned by AF3 and may not be consistent
across conditions. The same biological chain may get different letter
IDs in different predictions. The pipeline uses the chain letter as a
label; it does not assume biological identity across conditions.

---

## 4. Derived Tables

### 4a. `metrics_replicates.csv` (one row per prediction)

**Unit of observation:** One prediction (seed × sample)

**Columns:** `condition_id`, `condition_name`, `replicate_id`, all 13 global metrics, plus `chain_{ID}_plddt` and `chain_{ID}_residues` for each chain present in that prediction.

This is the primary per-observation table. Every row represents one
independent AF3 prediction.

### 4b. `seed_aggregated.csv` (one row per condition × seed)

**Unit of observation:** One seed within one condition (samples averaged)

**Columns:** `condition_id`, `condition_name`, `seed`, plus the mean of
each numeric column from `metrics_replicates.csv` across the samples
within that seed.

This is the main table for seed-level analysis. Each value represents
the average of the samples from the same random seed.

**Note:** `chain_{ID}_residues` is constant across seeds (same
sequence), so its mean is identical to the raw value.

### 4c. `descriptive_stats.csv` (one row per condition)

**Unit of observation:** One experimental condition

**Columns:** `condition_id`, plus `{metric}_mean` and `{metric}_sd` for
each numeric metric in `seed_aggregated.csv`.

The `_mean` is the average across all seeds within that condition.
The `_sd` is the between-seed standard deviation (ddof=1).

---

## 5. Metric Properties Summary

| Metric | Scale | Type | Direction | Aggregation Level | Primary Use |
|--------|-------|------|-----------|-------------------|-------------|
| pLDDT_mean | 0–100 | continuous | higher = better | per-prediction | Overall structural confidence |
| pLDDT_max | 0–100 | continuous | higher = better | per-prediction | Ceiling of confidence |
| pLDDT_min | 0–100 | continuous | higher = better | per-prediction | Floor of confidence |
| pLDDT_median | 0–100 | continuous | higher = better | per-prediction | Robust central tendency |
| contact_prob_mean | 0–1 | continuous | context-dependent | per-prediction | Contact density |
| contact_prob_max | 0–1 | continuous | higher = stronger contact | per-prediction | Strongest predicted contact |
| contact_prob_min | 0–1 | continuous | lower = weaker contact | per-prediction | Weakest predicted contact |
| contact_prob_median | 0–1 | continuous | context-dependent | per-prediction | Robust contact density |
| pae_mean | 0–~50 Å | continuous | lower = better | per-prediction | Overall positional accuracy |
| pae_max | 0–~50 Å | continuous | lower = better | per-prediction | Worst positional error |
| pae_min | 0–~50 Å | continuous | lower = better | per-prediction | Best positional error |
| pae_median | 0–~50 Å | continuous | lower = better | per-prediction | Robust positional accuracy |
| ranking_score | 0–1 | continuous | higher = better | per-prediction | Composite model quality |
| chain_{ID}_plddt | 0–100 | continuous | higher = better | per-chain per-prediction | Per-molecule confidence |
| chain_{ID}_residues | integer | count | N/A | per-chain (constant) | Sequence length |

---

## 6. Key Statistical Considerations

### 6a. Dependence Structure

- **Within-seed:** Samples from the same seed share the same initial
  random state and are NOT independent. They represent stochastic
  variation around a common starting point.
- **Between-seed:** Seeds are independent draws. Seeds are the primary
  independent units for between-condition inference.
- **Between-condition:** Conditions are independent experimental
  setups. Comparisons across conditions are the primary scientific
  interest.

The exact number of seeds and samples per seed is recorded in the
experiment metadata and should be checked before analysis.

### 6b. Missingness

- Some conditions may have different chain structures (e.g., a complex
  with and without a ligand), producing different columns in
  `metrics_replicates.csv`. Chain columns that don't exist for a
  prediction are implicitly NaN.
- The `ranking_score` may be missing for aggregate files (no
  seed/sample pairing).

### 6c. Scale Heterogeneity

Metrics span very different scales:
- pLDDT: 0–100
- contact_prob: 0–1
- pae: 0–~50 Å
- ranking_score: 0–1 (but typically occupies a narrow range within any
  given dataset)

Any analysis comparing or combining metrics must account for scale
differences (standardization, separate analyses, etc.).

### 6d. Distribution Shape

- **pLDDT:** Right-skewed in well-predicted structures (most values
  high, tail toward low). Bimodal if the structure has both well-folded
  and disordered regions.
- **contact_prob:** Highly right-skewed (most pairs are non-contacts).
  The mean is heavily influenced by the sparse-contact majority.
- **pae:** Right-skewed. Most pairs have low PAE (well-positioned),
  with a tail of high-PAE pairs (poorly positioned).
- **ranking_score:** Distribution shape varies by dataset. Typically
  approximately normal within a condition, but range and spread depend
  on the prediction difficulty.

### 6e. Metric Correlations

Within a single prediction:
- pLDDT and PAE are negatively correlated (high pLDDT → low PAE)
- contact_prob and PAE are negatively correlated (high contact prob →
  low PAE between contacting residues)
- pLDDT and ranking_score are positively correlated
- chain_{ID}_plddt contributes to pLDDT_mean (weighted by chain size)

Across predictions within a condition:
- All metrics have between-seed variance. This variance is the primary
  measure of prediction stochasticity.
- Metrics may show correlated between-seed variation.

---

## 7. What These Metrics Are NOT

- **Not structural accuracy.** Confidence metrics measure AF3's
  *self-assessed* confidence, not ground-truth accuracy. High confidence
  does not guarantee correctness.
- **Not experimental data.** These are computational predictions, not
  measurements. Statistical tests on these metrics assess prediction
  *consistency*, not biological truth.
- **Not interchangeable.** pLDDT, PAE, and contact probability measure
  different aspects of prediction quality. They should not be treated as
  redundant measures of the same thing.
- **Not independent across chains.** In a complex, the confidence of one
  chain depends on the confidence of its interaction partners.

---

## 8. Recommended Analysis Directions

(Descriptive only — no methodology imposed.)

**Within-condition questions:**
- How variable are predictions across seeds? (between-seed SD)
- Which chains are consistently well-predicted vs variable?
- Are pLDDT and PAE concordant, or does one suggest higher confidence
  than the other?

**Between-condition questions:**
- How does each experimental attribute affect prediction confidence?
- Are there interaction effects between attributes?
- Are the ranking scores significantly different across conditions?

**Exploratory questions:**
- Which metrics are most correlated across seeds/conditions?
- Are there conditions where prediction quality is unexpectedly low?
- Do samples within a seed vary more for some conditions than others?

---

## 9. File Layout

```
{data_directory}/
├── {condition_name}/
│   ├── *_confidences.json           # Per-prediction confidence data
│   ├── *_summary_confidences.json   # Summary scores (not used by pipeline)
│   ├── *_ranking_scores.csv         # Tabular ranking scores
│   ├── *_data.json                  # AF3 input metadata
│   ├── *_model.cif                  # Predicted 3D structure
│   └── seed-{N}_sample-{M}/        # Per-prediction subdirectories
│       ├── *_confidences.json
│       └── *_summary_confidences.json
```

The extraction pipeline recursively finds all `*_confidences.json`
files, extracts metrics, and produces:

```
outputs/{condition}_condition_centric/
├── metrics_replicates.csv       # One row per prediction
├── metrics_conditions.csv       # One row per condition (aggregated)
├── condition_registry.csv       # Condition metadata
└── condition_manifest.csv       # Summary statistics per condition
```

After seed aggregation (`_stage_run_analysis`):

```
tables/
├── seed_aggregated.csv          # One row per condition × seed
├── descriptive_stats.csv        # One row per condition (mean + SD)
└── pairwise_comparisons.csv     # (currently empty — pending implementation)
```

---

*Document generated from the analysis pipeline and AF3 documentation.*
*See also: DATA_PIPELINE_AUDIT.md, VISUALIZATION_AUDIT.md, PIPELINE_ARCHITECTURE.md*
