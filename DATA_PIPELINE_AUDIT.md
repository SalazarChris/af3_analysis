# Data Pipeline Audit

**Date:** 2026-08-26
**Status:** Read-only audit — no files modified.

---

## 1. Data Flow Overview

```
Raw AF3 outputs (JSON/CSV/CIF per replicate)
    │
    ▼
Stage 1: Extraction (af3_condition_centric_extraction.py)
    │  Reads: *_confidences.json, *_ranking_scores.csv, *_data.json
    │  Writes: 4 CSVs to tables/
    │
    ▼
Stage 4: Analysis (pipeline.py _stage_run_analysis)
    │  Reads: metrics_replicates.csv
    │  Writes: seed_aggregated.csv, descriptive_stats.csv, pairwise_comparisons.csv
    │
    ▼
Stage 5: Visualization (orchestrator.py)
    │  Reads: seed_aggregated.csv, descriptive_stats.csv, pairwise_comparisons.csv
    │  Writes: 8 PNG figures
    │
    ▼
Output: 4 tables + 8 figures + manifest
```

---

## 2. Source Data (Raw AF3 Outputs)

**Location:** User-specified input directory (e.g., `testdata/pou2/`)

**Structure:**
```
<input_dir>/
├── <condition_name>/
│   ├── <condition>_seed-<N>_sample-<M>_confidences.json
│   ├── <condition>_seed-<N>_sample-<M>_data.json
│   ├── <condition>_seed-<N>_sample-<M>_summary_confidences.json
│   ├── <condition>_seed-<N>_sample-<M>_ranking_scores.csv
│   ├── <condition>_seed-<N>_sample-<M>_model.cif
│   ├── <condition>_confidences.json          (aggregate)
│   ├── <condition>_data.json                 (metadata)
│   ├── <condition>_summary_confidences.json  (aggregate)
│   └── <condition>_ranking_scores.csv        (all seeds/samples)
```

**Key fields in `*_confidences.json`:**
- `atom_plddts` — per-atom pLDDT scores
- `contact_probs` — contact probability matrix
- `pae` — predicted aligned error matrix
- `token_chain_ids` — chain assignments for tokens
- `atom_chain_ids` — chain assignments for atoms

**Key fields in `*_ranking_scores.csv`:**
- `seed`, `sample`, `ranking_score`

---

## 3. Stage 1: Extraction

**Entry function:** `extract_metrics()` in `af3inputbuilder/scripts/af3_condition_centric_extraction.py`

**Called by:** `pipeline.py` → `_stage_extract_metrics()` (lines 60-100), which loads the script via `importlib.util.spec_from_file_location` using a relative path from `pipeline.py` → `parents[3] / "af3inputbuilder" / "scripts" / ...`

### Input Files Read

| File pattern | How read |
|---|---|
| `*_confidences.json` | `input_dir.rglob("*_confidences.json")` — JSON parsed for `atom_plddts`, `contact_probs`, `pae`, `token_chain_ids`, `atom_chain_ids` |
| `*_ranking_scores.csv` | `input_dir.rglob("*_ranking_scores.csv")` — pandas `read_csv`, columns: `seed`, `sample`, `ranking_score` |
| `*_data.json` | `folder.glob("*_data.json")` — reads `name` field for condition display name |

### Output CSVs Written to `tables/`

#### `condition_registry.csv`
| Column | Type | Description |
|--------|------|-------------|
| `condition_id` | str | e.g. `cond_001` |
| `condition_name` | str | Human-readable name from `_data.json` or folder name |
| `n_replicates` | int | Number of replicates for this condition |
| `replicate_ids` | str (JSON list) | JSON-encoded list of replicate_id strings |

#### `metrics_replicates.csv` (PRIMARY TABLE)
| Column | Type | Description |
|--------|------|-------------|
| `condition_id` | str | e.g. `cond_001` |
| `condition_name` | str | Human-readable name |
| `replicate_id` | str | e.g. `pou_baseline_seed-1_sample-0` |
| `pLDDT_mean` | float | Mean of `atom_plddts` |
| `pLDDT_max` | float | Max of `atom_plddts` |
| `pLDDT_min` | float | Min of `atom_plddts` |
| `pLDDT_median` | float | Median of `atom_plddts` |
| `contact_prob_mean` | float | Mean of `contact_probs` |
| `contact_prob_max` | float | Max of `contact_probs` |
| `contact_prob_min` | float | Min of `contact_probs` |
| `contact_prob_median` | float | Median of `contact_probs` |
| `pae_mean` | float | Mean of `pae` |
| `pae_max` | float | Max of `pae` |
| `pae_min` | float | Min of `pae` |
| `pae_median` | float | Median of `pae` |
| `ranking_score` | float | From ranking_scores.csv (if available) |
| `chain_<X>_residues` | int | Token count for chain X (dynamic) |
| `chain_<X>_plddt` | float | Mean pLDDT for chain X (dynamic) |

**Note:** Chain columns are dynamic — their count and labels depend on the input data. Chain letters are derived from `sorted(token_counts.keys())`.

#### `metrics_conditions.csv`
| Column | Type | Description |
|--------|------|-------------|
| `condition_id` | str | |
| `condition_name` | str | |
| `n_replicates` | int | |
| `n_seeds` | int | |
| `<metric>_mean` | float | Mean across replicates |
| `<metric>_sd` | float | Std dev across replicates |
| `<metric>_median` | float | Median across replicates |
| `<metric>_min` | float | Min across replicates |
| `<metric>_max` | float | Max across replicates |
| `<metric>_cv` | float | Coefficient of variation |
| `chain_<X>_<metric>_mean` | float | Chain-level aggregated means |
| `chain_<X>_<metric>_sd` | float | Chain-level aggregated SDs |

#### `condition_manifest.csv`
| Column | Type | Description |
|--------|------|-------------|
| `condition_id` | str | |
| `condition_name` | str | |
| `replicates` | int | |
| `seeds` | int | |
| `status` | str | `Complete` or `Empty` |

---

## 4. Stage 4: Analysis (Seed Aggregation)

**Entry function:** `_stage_run_analysis()` in `pipeline.py` (lines 130-177)

### Input

Reads `tables/metrics_replicates.csv` via `pd.read_csv()`.

### Processing

1. Extracts `seed` column from `replicate_id` by parsing `seed-<N>` pattern
2. Groups by `(condition_id, condition_name, seed)`, taking `.mean()` of all numeric columns
3. Writes `seed_aggregated.csv`
4. Groups by `condition_id`, computing `mean` and `std` for condition-level summaries
5. Writes `descriptive_stats.csv`
6. Creates empty `pairwise_comparisons.csv` (placeholder)

### Output CSVs

#### `seed_aggregated.csv` (KEY TABLE FOR VISUALIZATION)
| Column | Type | Description |
|--------|------|-------------|
| `condition_id` | str | |
| `condition_name` | str | |
| `seed` | int | Extracted from replicate_id |
| `pLDDT_mean` | float | Seed-level mean of pLDDT_mean |
| `pLDDT_min` | float | Seed-level mean of pLDDT_min |
| `pae_mean` | float | Seed-level mean of PAE mean |
| `contact_prob_mean` | float | Seed-level mean of contact probability |
| `ranking_score` | float | Seed-level mean of ranking score |
| `chain_<X>_plddt` | float | Seed-level mean of chain X pLDDT |
| ... (all numeric columns from metrics_replicates) | | |

**Row count:** One row per (condition × seed). For 8 conditions × 10 seeds = 80 rows.

#### `descriptive_stats.csv`
| Column | Type | Description |
|--------|------|-------------|
| `condition_id` | str | |
| `<metric>_mean` | float | Mean across seeds |
| `<metric>_std` | float | SD across seeds |
| ... (all numeric metrics × {mean, std}) | | |

**Note:** `condition_name` is NOT included in this CSV. The orchestrator adds it back by mapping from `seed_aggregated.csv`.

#### `pairwise_comparisons.csv` (EMPTY)
| Column | Type | Description |
|--------|------|-------------|
| `reference` | str | Empty |
| `condition` | str | Empty |
| `metric` | str | Empty |
| `diff_mean` | float | Empty |

---

## 5. Stage 5: Visualization

**Entry function:** `generate_all_figures()` in `visualization/orchestrator.py`

### Input

Reads three CSVs from `tables/`:
```python
seed_aggregated_df = pd.read_csv(tables_dir / "seed_aggregated.csv")
descriptive_stats_df = pd.read_csv(tables_dir / "descriptive_stats.csv")
pairwise_comparisons_df = pd.read_csv(tables_dir / "pairwise_comparisons.csv")
```

### Schema Inference

The orchestrator builds a `schema` dict dynamically:

```python
# Global metrics — HARD-CODED LIST
global_metrics = ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']

# Chain metrics — DYNAMIC (from column names)
chain_metrics = [c for c in seed_agg.columns if c.startswith('chain_') and c.endswith('_plddt')]

valid_metrics = [m for m in global_metrics + chain_metrics if m in seed_agg.columns]

# Factors — DYNAMIC with hard-coded fallback
factors = [c for c in seed_agg.columns if c.startswith('factor_')]
if not factors:
    factors = ['factor_DNA', 'factor_PTM']  # fallback
```

### Figure Mapping

| Plot # | Function | Input DataFrame | Key columns used |
|--------|----------|-----------------|------------------|
| F1 | `plot_qc_completeness` | `seed_aggregated` | `condition_name`, all metrics (notna fraction) |
| F2 | `plot_seed_distributions` | `seed_aggregated` | `condition_name`, all metrics (box + strip) |
| F3 | `plot_factorial_interaction` | `descriptive_stats` | `condition_name`, `{metric}_mean`, `factor_PTM`, `factor_DNA` |
| F4 | `plot_effect_size_forest` | `pairwise_comparisons` | `metric`, `condition`, `reference`, `diff_mean` |
| F5 | `plot_variability` | `descriptive_stats` | `condition_name`, `{metric}_std` |
| F6 | `plot_seed_trajectories` | `seed_aggregated` | `condition_name`, `seed`, all metrics (line per seed) |
| F7 | `plot_ecdf_overlay` | `seed_aggregated` | `condition_name`, all metrics (ECDF per condition) |
| F8 | `plot_metric_relationships` | `seed_aggregated` | `condition_name`, subset of metrics (pairplot) |

### Display Labels

Hard-coded in `visualization/utils.py`:
- Condition labels: `pou_baseline` → "Baseline (POU)", etc. (8 conditions)
- Metric labels: `plddt_mean` → "Mean pLDDT", etc.
- Chain labels: `chain_A_plddt` → "Chain A pLDDT"

### Condition Ordering

Hard-coded semantic rank in `utils.py`:
```python
ptm_rank = {
    'pou_baseline': 0, 'pou_sep102': 1, 'pou_tpo101': 2,
    'pou_tpo101_sep102': 3, 'pou_dna': 4, 'pou_sep102_dna': 5,
    'pou_tpo101_dna': 6, 'pou_tpo101_sep102_dna': 7
}
```

---

## 6. Downstream Modules NOT Currently Wired

The following modules exist in `af3_analysis/` but are **not called** by the current pipeline:

### `statistics/` — Requires long-format analysis_seed DataFrame

| Module | Function | Expected input |
|--------|----------|----------------|
| `descriptive.py` | `compute_condition_summaries(analysis_seed_df)` | Columns: `condition_id, seed, metric_id, scope_type, scope_id, value, sample_sd` |
| `descriptive.py` | `compute_variability_summaries(analysis_seed_df)` | Same as above |
| `resampling.py` | `Resampler.bootstrap_mean(values)` | numpy array of seed-level values |
| `resampling.py` | `Resampler.permutation_test(values_a, values_b)` | Two numpy arrays |
| `comparisons.py` | `compare_two_conditions(cond_a, cond_b, values_a, values_b)` | Two numpy arrays of seed-level values |
| `variance.py` | `decompose_variance(analysis_seed_df, metric_id)` | Long-format with `condition_id, seed, metric_id, scope_type, scope_id, value, sample_sd` |
| `multiplicity.py` | `holm_correction(p_values, ...)` | Lists of p-values |
| `eligibility.py` | `check_stage_eligibility(stage, ...)` | conditions_df, replicates_df, measurements_df |
| `factorial_models.py` | `run_factorial_analysis(seed_data)` | DataFrame with `SeedID` column |

### `exploratory/` — Requires long-format measurements DataFrame

| Module | Function | Expected input |
|--------|----------|----------------|
| `inventory.py` | `DesignInventory` | Factor-level metadata |
| `variance.py` | `seed_means(long_df)` | Columns: `condition_id, seed, sample, metric_id, scope_type, scope_id, value` |
| `variance.py` | `variance_components(long_df)` | Same as above |
| `variance.py` | `degenerate_metrics(seed_means_df)` | Seed-level means |
| `distributions.py` | `summarise_distributions(long_df)` | Long-format measurements |
| `factors.py` | `factor_marginals(df, metric_col, factor_cols)` | Wide or long with factor columns |

### `preprocessing/` — Requires raw AF3 artifacts or legacy CSVs

| Module | Function | Expected input |
|--------|----------|----------------|
| `metadata.py` | `load_condition_metadata(legacy_csv_root)` | Directory with condition_manifest.csv |
| `canonicalize.py` | `Canonicalizer(artifacts)` | List of ArtifactRecord from raw reader |
| `quality_control.py` | `run_quality_control(replicates_df, measurements_df, conditions_df)` | Canonical tables |
| `aggregation.py` | `aggregate_to_seed_level(measurements_df, replicates_df)` | Long-format measurements + replicates |
| `duplicate_resolution.py` | `resolve_duplicates(records)` | List of record dicts |
| `mappings.py` | `MappingBuilder` | Chain/residue mapping data |
| `definedness.py` | `classify_definedness(measurements_df)` | Long-format measurements |

### `io/` — Requires raw AF3 files or legacy CSVs

| Module | Function | Expected input |
|--------|----------|----------------|
| `raw_af3_reader.py` | `RawAF3Reader(raw_af3_root)` | Root directory with AF3 outputs |
| `phase1_loader.py` | `Phase1Loader(legacy_csv_root)` | Directory with 4 legacy CSVs |
| `parquet_store.py` | `ParquetStore(store_dir)` | Directory for Parquet files |
| `artifact_inventory.py` | `ArtifactInventory(root, dataframes)` | Legacy CSVs + root path |

### `registry/` — Requires registry JSON or DataFrame

| Module | Function | Expected input |
|--------|----------|----------------|
| `metric_registry.py` | `MetricRegistry.load(source)` | JSON file with `{"metrics": [...]}` |
| `resolution.py` | `resolve_metrics(registry, columns)` | Registry + column names |

---

## 7. Data Compatibility Analysis

### Current pipeline produces (wide format):

```
seed_aggregated.csv:
  condition_id | condition_name | seed | pLDDT_mean | pae_mean | chain_A_plddt | ...
  cond_001     | pou_baseline   | 1    | 78.5       | 3.2      | 82.1         | ...
```

### Downstream modules expect (long format):

```
analysis_seed DataFrame:
  condition_id | seed | metric_id     | scope_type | scope_id | value | sample_sd
  cond_001     | 1    | pLDDT_mean    | global     |          | 78.5  | 2.1
  cond_001     | 1    | chain_A_plddt | chain      | A        | 82.1  | 1.5
```

### Compatibility verdict

| Interface | Compatible? | Notes |
|-----------|-------------|-------|
| Current extraction → pipeline analysis | ✅ | Both use `metrics_replicates.csv` with same columns |
| Current pipeline → visualization | ✅ | Orchestrator reads the 3 CSVs it expects |
| Current pipeline → `statistics/` | ❌ | Wide vs long format mismatch |
| Current pipeline → `exploratory/` | ❌ | Wide vs long format mismatch |
| Current pipeline → `preprocessing/` | ❌ | Different data structures entirely |
| Current pipeline → `registry/` | ❌ | No registry JSON exists; hardcoded in `schemas/tables.py` |

---

## 8. Problems Identified

### P1: Hard-coded metric list in orchestrator (MEDIUM)

**File:** `visualization/orchestrator.py`, line 46:
```python
global_metrics = ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']
```

This excludes `pLDDT_max`, `pLDDT_median`, `contact_prob_max/min/median`, `pae_max/min/median` from all plots. If new metrics are added, this list must be manually updated.

**Impact:** Some available metrics are silently excluded from visualization.

### P2: Hard-coded factor fallback (LOW)

**File:** `visualization/orchestrator.py`, lines 52-54:
```python
factors = [c for c in seed_agg.columns if c.startswith('factor_')]
if not factors:
    factors = ['factor_DNA', 'factor_PTM']
```

The fallback assumes the POU experimental design. If conditions don't match the `*_dna` / `*_tpo101` naming pattern, the factorial plot will fail or produce meaningless results.

**Impact:** Visualization is coupled to a specific experimental design.

### P3: Hard-coded condition ordering and labels (LOW)

**File:** `visualization/utils.py` — `ptm_rank` dict and `labels['condition']` dict.

The 8 POU conditions are hard-coded. Any new condition gets rank 99 (sorted last) and a generic fallback label.

**Impact:** Works for current data but not extensible.

### P4: pairwise_comparisons.csv is always empty (HIGH)

**File:** `pipeline.py`, line 177:
```python
pd.DataFrame({'reference': [], 'condition': [], 'metric': [], 'diff_mean': []}).to_csv(...)
```

This means F4 (Effect Size Forest) never produces a figure — `plot_effect_size_forest` returns immediately on `len(pairwise_comparisons_df) == 0`.

**Impact:** One of the 8 planned figures is permanently empty.

### P5: fragile path to extraction script (HIGH)

**File:** `pipeline.py`, lines 71-76:
```python
script_path = (
    Path(__file__).resolve().parents[3]
    / "af3inputbuilder"
    / "scripts"
    / "af3_condition_centric_extraction.py"
)
```

This navigates up 3 parent directories from `af3_analysis/pipeline.py`. If the directory structure changes (e.g., af3_analysis is moved, or the repo is restructured), this breaks silently.

**Impact:** Pipeline stage 1 fails if directory nesting changes.

### P6: descriptive_stats.csv loses condition_name (MEDIUM)

**File:** `pipeline.py`, lines 168-172. The `groupby('condition_id')` drops `condition_name`. The orchestrator patches this by re-joining from `seed_aggregated.csv` (line 38-39), but this is a workaround, not a fix.

**Impact:** Works but fragile — if `seed_aggregated.csv` changes, the join breaks.

### P7: seed extraction from replicate_id is fragile (MEDIUM)

**File:** `pipeline.py`, lines 147-152:
```python
def extract_seed(rep_id):
    for part in str(rep_id).split('_'):
        if part.startswith('seed-'):
            try: return int(part.split('-')[1])
            except: pass
    return 1  # Fallback
```

The fallback `return 1` silently assigns all non-matching replicate_ids to seed 1, which would corrupt seed-level aggregation.

**Impact:** Silent data corruption if replicate_id format varies.

### P8: chain column detection is implicit (LOW)

The extraction script creates `chain_<X>_residues` and `chain_<X>_plddt` columns dynamically. The orchestrator filters for `chain_*_plddt` but not `chain_*_residues`. This works but the contract is implicit.

**Impact:** Low risk — but not documented.

### P9: No `condition_name` in `metrics_replicates.csv` join to `seed_aggregated` (LOW)

The pipeline's `_stage_run_analysis` reads `metrics_replicates.csv` which includes `condition_name`. The `groupby` preserves it because `condition_name` is in `available_groupby`. This works correctly but is not explicitly verified.

**Impact:** None currently — but implicit coupling.

---

## 9. What CAN Be Reused by Future Modules

### Immediate reuse (no changes needed):

| Data | Format | Consumers |
|------|--------|-----------|
| `seed_aggregated.csv` | Wide, one row per condition × seed | All current plots (F1, F2, F6, F7, F8) |
| `descriptive_stats.csv` | Wide, one row per condition | Plots F3, F5 |
| `metrics_replicates.csv` | Wide, one row per replicate | Can be melted to long format for statistics modules |

### Requires transformation (melt/pivot):

| Source | Target | Transformation |
|--------|--------|----------------|
| `metrics_replicates.csv` (wide) | `analysis_seed_df` (long) | Melt on metric columns → `(condition_id, seed, metric_id, value)` |
| `seed_aggregated.csv` (wide) | `long_df` (long) | Melt on metric columns → `(condition_id, seed, metric_id, value)` |

The melt operation is straightforward:
```python
id_vars = ['condition_id', 'condition_name', 'seed']
metric_vars = [c for c in df.columns if c not in id_vars]
long_df = seed_aggregated.melt(id_vars=id_vars, value_vars=metric_vars,
                                var_name='metric_id', value_name='value')
```

### Can be added without touching existing pipeline:

- `statistics/descriptive.py` — compute condition summaries from seed_aggregated
- `statistics/resampling.py` — bootstrap CIs from seed-level values
- `statistics/comparisons.py` — pairwise comparisons (would populate pairwise_comparisons.csv)
- `statistics/variance.py` — variance decomposition
- `exploratory/variance.py` — seed means, variance components
- `exploratory/distributions.py` — distribution summaries
- `exploratory/factors.py` — factor marginals

---

## 10. Recommendation

### **B. REUSE WITH SMALL MODIFICATIONS**

**Rationale:**

1. **The core data flow works.** Raw AF3 → extraction CSVs → seed aggregation → visualization is functional and produces correct results for the current dataset.

2. **The visualization data is reusable as-is.** The three CSVs consumed by the orchestrator (`seed_aggregated.csv`, `descriptive_stats.csv`, `pairwise_comparisons.csv`) can be consumed directly by new visualization and statistics modules after a simple wide-to-long melt.

3. **The problems are fixable without refactoring:**
   - P1 (hard-coded metrics): Add dynamic metric discovery from columns
   - P4 (empty pairwise comparisons): Implement real comparisons in pipeline
   - P5 (fragile path): Use a config or relative import instead of `parents[3]`
   - P6 (missing condition_name): Include it in the groupby output

4. **The downstream modules are well-designed.** They expect a clean long-format interface that can be produced from the existing CSVs with a one-line melt operation.

**Recommended modifications (in priority order):**

1. Fix P5 (fragile extraction path) — critical for reliability
2. Fix P4 (implement pairwise comparisons) — enables F4 figure
3. Add a wide-to-long converter function (5 lines) for downstream modules
4. Make metric discovery dynamic in orchestrator (replace hard-coded list)
5. Make factor inference configurable (pass experimental design metadata)

**Do NOT:**
- Replace the extraction script
- Change the CSV output format
- Modify the orchestrator's load_data() interface
- Change the visualization plotting functions
