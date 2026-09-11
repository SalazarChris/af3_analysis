# Pipeline Architecture

**Date:** 2026-08-26

---

## Overview

```
Raw AF3 data (JSON/CSV/CIF per replicate)
          │
          ▼
    ┌─────────────────────────────┐
    │   EXTRACTION                │
    │   (af3_condition_centric_   │
    │    extraction.py)           │
    │                             │
    │   Reads: *_confidences.json │
    │          *_ranking_scores   │
    │          *_data.json        │
    │   Writes: 4 CSVs to tables/ │
    └─────────────┬───────────────┘
                  │
                  ▼
    ┌─────────────────────────────┐
    │   CANONICAL ANALYSIS TABLES │
    │   (pipeline.py)             │
    │                             │
    │   metrics_replicates.csv    │  ← raw per-replicate metrics
    │   seed_aggregated.csv       │  ← one row per condition × seed
    │   descriptive_stats.csv     │  ← condition-level mean/SD
    │   pairwise_comparisons.csv  │  ← placeholder (empty)
    └─────────────┬───────────────┘
                  │
          ┌───────┴───────┐
          │               │
          ▼               ▼
    ┌──────────┐   ┌──────────────────┐
    │ METADATA │   │ VISUALIZATION    │
    │ (JSON)   │   │ (orchestrator.py │
    │          │   │  + core_plots.py)│
    │          │   │                  │
    │ Loads:   │   │ Reads:           │
    │ 3 CSVs   │   │  3 CSVs          │
    │ +        │   │  + metadata      │
    │ metadata │   │                  │
    └──────────┘   │ Writes:          │
                   │  8 PNG figures   │
                   └──────────────────┘
                          │
                          ▼
                   ┌──────────────────┐
                   │ FUTURE:          │
                   │ Statistical      │
                   │ Analysis         │
                   │ (long format)    │
                   └──────────────────┘
```

---

## Component Categories

### 1. Data Processing

| Component | Location | Responsibility |
|-----------|----------|---------------|
| Extraction | `af3inputbuilder/scripts/af3_condition_centric_extraction.py` | Reads raw AF3 outputs, produces 4 CSV tables |
| Seed aggregation | `af3_analysis/pipeline.py` → `_stage_run_analysis()` | Aggregates replicates to seed level, computes descriptive stats |
| Wide-to-long adapter | `af3_analysis/io/wide_to_long.py` | Converts wide CSVs to long format for statistics modules |

**Key principle:** These components process DATA VALUES only. They do not know about experimental design.

### 2. Experimental Design Metadata

| Component | Location | Responsibility |
|-----------|----------|---------------|
| Metadata file | `testdata/pou2/experiment_metadata.json` (per experiment) | Defines conditions and their attributes |
| Metadata module | `af3_analysis/experiment_metadata.py` | Loads, validates, inspects experiment design |
| Design inspection | `inspect_design()` in `experiment_metadata.py` | Describes design structure (factorial, incomplete, etc.) |

**Key principle:** Metadata describes the EXPERIMENTAL DESIGN, not the data. It answers: "What conditions exist and what attributes define them?"

### 3. Visualization

| Component | Location | Responsibility |
|-----------|----------|---------------|
| Orchestrator | `af3_analysis/visualization/orchestrator.py` | Loads data + metadata, calls plotting functions |
| Core plots | `af3_analysis/visualization/core_plots.py` | 8 plotting functions (F1–F8) |
| Utilities | `af3_analysis/visualization/utils.py` | Display labels, condition ordering, color palettes |

**Key principle:** Visualization consumes BOTH data AND metadata. The metadata provides labels, ordering, and factor structure. The data provides numerical values.

### 4. Statistical Infrastructure (future)

| Component | Location | Responsibility |
|-----------|----------|---------------|
| Descriptive | `af3_analysis/statistics/descriptive.py` | Condition-level summaries |
| Resampling | `af3_analysis/statistics/resampling.py` | Bootstrap, permutation tests |
| Comparisons | `af3_analysis/statistics/comparisons.py` | Pairwise condition comparisons |
| Variance | `af3_analysis/statistics/variance.py` | Variance decomposition |
| Factorial | `af3_analysis/statistics/factorial_models.py` | Design matrix, model fitting |

**Key principle:** Statistical modules consume LONG-FORMAT DataFrames. They do not need metadata directly — they operate on `(condition_id, seed, metric_id, value)` tuples.

### 5. Statistical Methodology (future)

Not implemented. When implemented, this category will contain:
- Hypothesis test selection
- Effect size computation
- Confidence interval computation
- Multiple comparison correction

**Key principle:** Statistical methodology choices are made by the researcher, not inferred from data.

---

## Data Flow Details

### Stage 1: Extraction

**Input:** Raw AF3 output directory containing per-replicate JSON/CSV files.

**Output:** 4 CSVs in `tables/`:

| File | Row grain | Key columns |
|------|-----------|-------------|
| `condition_registry.csv` | 1 per condition | `condition_id`, `condition_name`, `n_replicates` |
| `metrics_replicates.csv` | 1 per replicate | `condition_id`, `condition_name`, `replicate_id`, all metrics |
| `metrics_conditions.csv` | 1 per condition | `condition_id`, `condition_name`, aggregated metrics |
| `condition_manifest.csv` | 1 per condition | `condition_id`, `condition_name`, `seeds`, `status` |

**Design-agnostic:** The extraction knows nothing about experimental design. It discovers conditions from folder structure and names from `_data.json`.

### Stage 2: Analysis (Seed Aggregation)

**Input:** `metrics_replicates.csv`

**Processing:**
1. Parse seed from `replicate_id` (raises on invalid patterns)
2. Group by `(condition_id, condition_name, seed)`, take mean of numeric columns
3. Group by `condition_id`, compute mean/SD of seed means

**Output:**

| File | Row grain | Key columns |
|------|-----------|-------------|
| `seed_aggregated.csv` | 1 per condition × seed | `condition_id`, `condition_name`, `seed`, all metrics |
| `descriptive_stats.csv` | 1 per condition | `condition_id`, `<metric>_mean`, `<metric>_std` |
| `pairwise_comparisons.csv` | (empty placeholder) | `reference`, `condition`, `metric`, `diff_mean` |

**Design-agnostic:** This stage computes summary statistics. It does not know about experimental design.

### Stage 3: Visualization

**Input:** 3 CSVs + optional experiment metadata JSON

**Processing:**
1. Load CSVs into DataFrames
2. Load metadata (if provided)
3. Build schema dict with metrics, condition order, and design metadata
4. Call each plotting function with data + schema

**Output:** 8 PNG figures

**Design-aware:** The visualization layer uses metadata for:
- Condition labels (from `metadata.conditions[name].label`)
- Condition ordering (from `metadata.condition_names`)
- Factor structure (from `metadata.conditions[name].attributes`)
- Design description (from `inspect_design()`)

---

## Data vs. Metadata Separation

### What is DATA

- Condition identifiers (`condition_id`, `condition_name`)
- Seed identifiers
- Metric values (`pLDDT_mean`, `pae_mean`, etc.)
- Chain-level metrics
- Ranking scores

### What is METADATA

- Condition labels ("Baseline (POU)", "POU + DNA")
- Attribute definitions (DNA: binary, pTPO101: binary)
- Attribute values per condition (pou_baseline: {DNA: false, ...})
- Design structure (complete factorial, incomplete, etc.)

### What is METHODOLOGY (future)

- Statistical test selection
- Effect size computation
- Confidence interval parameters
- Multiple comparison correction
- Significance thresholds

---

## Metadata Schema

The experiment metadata is a JSON file with this structure:

```json
{
  "experiment_id": "unique_id",
  "description": "Human-readable description",
  "attributes": {
    "attribute_name": {
      "type": "binary" | "categorical",
      "description": "Optional description"
    }
  },
  "conditions": {
    "condition_name": {
      "label": "Display label",
      "attributes": {
        "attribute_name": value
      }
    }
  }
}
```

### Rules

1. **All defined attributes must be present** in every condition
2. **No undefined attributes** may appear in conditions
3. **Labels must be unique** across conditions
4. **Binary attributes** must have `true`/`false` values
5. **Categorical attributes** must have string values
6. **No inference from names** — all structure is explicit

---

## Visualization Design Principles

### F3 (Factorial Interaction)

F3 adapts to the experimental design via metadata:

| Design | F3 behavior |
|--------|-------------|
| 2+ attributes | x-axis = first attribute, hue = second attribute |
| 1 attribute | x-axis = attribute, bars for each level |
| No metadata | Alphabetical condition display, no factor interpretation |

F3 does NOT assume:
- Binary attributes
- Two factors
- Complete factorial design
- Any specific biological entities

### All other figures (F1, F2, F5, F6, F7, F8)

These figures consume data only (CSVs) and use metadata for:
- Condition labels (display names)
- Condition ordering (consistent across figures)
- Color palette assignment

They do NOT use factor/attribute information.

---

## File Modification Summary

### Files created in this task

| File | Purpose |
|------|---------|
| `EXPERIMENT_DESIGN_AUDIT.md` | Phase 1 audit of hard-coded design dependencies |
| `testdata/pou2/experiment_metadata.json` | POU experiment metadata |
| `af3_analysis/experiment_metadata.py` | Metadata loading, validation, inspection |
| `af3_analysis/tests/test_experiment_metadata.py` | 40 tests for metadata system |

### Files modified in this task

| File | Change |
|------|--------|
| `af3_analysis/visualization/orchestrator.py` | Accept `metadata_path`/`experiment_design` parameter, pass metadata to schema |
| `af3_analysis/visualization/utils.py` | Remove POU-specific hard-coding; use metadata for labels/ordering; use colorblind-friendly palette |
| `af3_analysis/visualization/core_plots.py` | F3 rewritten to use metadata attributes; all functions pass `schema` to label utilities |

### Files NOT modified

| File | Reason |
|------|--------|
| `af3inputbuilder/scripts/af3_condition_centric_extraction.py` | Extraction is design-agnostic |
| `af3_analysis/pipeline.py` | Analysis tables are design-agnostic |
| `af3_analysis/io/wide_to_long.py` | Data adapter is design-agnostic |
| `af3_analysis/statistics/` | Statistical modules are design-agnostic |
| `af3_analysis/exploratory/` | Exploratory modules are design-agnostic |
| `af3_analysis/schemas/` | Schema system is independent |

---

## Known Issues

### 1. Statistics module naming conflict (pre-existing)

`af3_analysis/statistics/` shadows Python's stdlib `statistics` module, causing seaborn imports to fail. This is documented in `VISUALIZATION_AUDIT.md` and must be resolved separately.

### 2. Pairwise comparisons always empty (pre-existing)

`pairwise_comparisons.csv` is always empty because real comparisons are not implemented. F4 (forest plot) never renders.

### 3. Descriptive stats loses condition_name (pre-existing)

`descriptive_stats.csv` drops `condition_name` during groupby. The orchestrator patches this with a join.

### 4. Hard-coded metric list in orchestrator (pre-existing)

`global_metrics` in orchestrator.py excludes some available metrics from visualization.

---

## Future Directions

### Near-term

1. Resolve the `statistics` package naming conflict
2. Implement real pairwise comparisons (populate F4)
3. Add `condition_name` to `descriptive_stats.csv` output

### Medium-term

4. Wire the `statistics/` and `exploratory/` modules to the pipeline
5. Add design-aware statistical models (using metadata)
6. Update F6 (seed trajectories) to use a heatmap instead of line connections

### Long-term

7. Add a CLI flag to specify metadata path
8. Auto-discover metadata from data directory
9. Support multi-experiment analyses
