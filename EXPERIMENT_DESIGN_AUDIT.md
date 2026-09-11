# Experiment Design Audit

**Date:** 2026-08-26
**Status:** Read-only audit — no files modified.

---

## 1. Current Representation

### 1.1 Where condition information originates

Condition identity enters the system through **filesystem structure**:

```
testdata/pou2/
├── pou_baseline/
├── pou_dna/
├── pou_sep102/
├── pou_sep102_dna/
├── pou_tpo101/
├── pou_tpo101_dna/
├── pou_tpo101_sep102/
├── pou_tpo101_sep102_dna/
└── metrics_example/
```

Each top-level subdirectory is treated as one condition. The extraction script (`af3_condition_centric_extraction.py`) discovers conditions by iterating `input_dir.iterdir()` and matching subdirectories.

**Condition name resolution** (extraction, lines 54-68):
1. Looks for `*_data.json` inside the condition folder
2. Reads the `name` field from that JSON
3. Falls back to the folder name if no JSON or `name` field exists

**There is no metadata file describing what each condition represents.** The condition name is treated as a opaque string. The only structure extracted is:
- `condition_id`: auto-assigned as `cond_001`, `cond_002`, … (sorted by folder name)
- `condition_name`: the string from `_data.json` or folder name

### 1.2 How condition names propagate

```
filesystem folders
    ↓ extraction
condition_registry.csv     (condition_id, condition_name, n_replicates, replicate_ids)
metrics_replicates.csv     (condition_id, condition_name, replicate_id, metrics...)
    ↓ seed aggregation
seed_aggregated.csv        (condition_id, condition_name, seed, metrics...)
descriptive_stats.csv      (condition_id, <metric>_mean, <metric>_std, ...)  ← condition_name DROPPED
    ↓ visualization patching
descriptive_stats_df gets condition_name re-joined from seed_aggregated.csv
    ↓
plot functions receive condition_name
```

### 1.3 What is NOT tracked

- No attribute/modification/component metadata per condition
- No experimental design description
- No factor definitions
- No mapping from condition → components
- No indication of which combinations are observed vs. absent

---

## 2. Hard-Coded Design Dependencies

### 2.1 `visualization/utils.py` — Display labels (lines 20-32)

Hard-coded condition-name → display-name mapping:

```python
labels = {
    'condition': {
        'pou_baseline': 'Baseline (POU)',
        'pou_dna': 'POU + DNA',
        'pou_sep102': 'POU (pSEP102)',
        'pou_sep102_dna': 'POU (pSEP102) + DNA',
        'pou_tpo101': 'POU (pTPO101)',
        'pou_tpo101_dna': 'POU (pTPO101) + DNA',
        'pou_tpo101_sep102': 'POU (pTPO101 + pSEP102)',
        'pou_tpo101_sep102_dna': 'POU (pTPO101 + pSEP102) + DNA',
    }
}
```

**Impact:** Any condition not in this dict gets a generic fallback label (line 64: `identifier.replace('pou_', '').replace('_', ' ').title()`). This works but produces less meaningful labels.

**Files affected:** `core_plots.py` (all 8 plot functions call `get_display_label`)

### 2.2 `visualization/utils.py` — Condition ordering (lines 68-80)

Hard-coded semantic rank:

```python
ptm_rank = {
    'pou_baseline': 0,
    'pou_sep102': 1,
    'pou_tpo101': 2,
    'pou_tpo101_sep102': 3,
    'pou_dna': 4,
    'pou_sep102_dna': 5,
    'pou_tpo101_dna': 6,
    'pou_tpo101_sep102_dna': 7
}
```

**Impact:** Unknown conditions get rank 99 (sorted last). This is acceptable as a fallback but produces an undefined ordering for new experiments.

**Files affected:** All 8 plot functions (via `get_condition_order()`)

### 2.3 `visualization/orchestrator.py` — Factor inference fallback (lines 52-54)

```python
factors = [c for c in seed_agg.columns if c.startswith('factor_')]
if not factors:
    factors = ['factor_DNA', 'factor_PTM']  # ← POU-specific
```

**Impact:** If no `factor_` columns exist in the data (they never do in current output), the code assumes `factor_DNA` and `factor_PTM` exist. These are never actually added to the CSV, so the fallback is the only path taken.

**Files affected:** `orchestrator.py`, which passes `schema['factors']` to `plot_factorial_interaction`

### 2.4 `visualization/core_plots.py` — F3 factor inference from condition names (lines 110-130)

```python
if not factor_cols:
    for name in df['condition_name']:
        dna = 1 if name.endswith('_dna') else 0
        clean_name = name.replace('_dna', '')
        if clean_name == 'pou_baseline':
            ptm = 'baseline'
        elif clean_name == 'pou':
            ptm = 'pou'
        else:
            ptm = clean_name.replace('pou_', '')
        ptms.append(ptm)
        dnas.append(dna)
    df['factor_PTM'] = ptms
    df['factor_DNA'] = dnas
```

This derives two binary factors from the condition name:
- `factor_DNA`: whether the name ends with `_dna`
- `factor_PTM`: the string after stripping `pou_` prefix and `_dna` suffix

**Impact:** This is the most severe hard-coding. It assumes:
1. Conditions contain `pou_` prefix
2. DNA presence is indicated by `_dna` suffix
3. Everything else is a PTM modification
4. All factors are binary
5. There are exactly 2 factors

**Files affected:** `core_plots.py` (F3 only)

### 2.5 `visualization/core_plots.py` — F3 PTM pretty labels (lines 140-146)

```python
ptm_pretty = {
    'baseline': 'Baseline (POU)',
    'pou': 'POU Only',
    'sep102': 'pSEP102 modification',
    'tpo101': 'pTPO101 modification',
    'tpo101_sep102': 'pTPO101 + pSEP102 modification'
}
```

**Impact:** Hard-coded display labels for PTM states. These are POU-specific.

### 2.6 `visualization/core_plots.py` — F8 metric list (lines 182-186)

```python
f8_metrics = ['pLDDT_mean', 'pae_mean', 'contact_prob_mean', 'ranking_score']
```

**Impact:** Ignores metrics from the schema and uses a hard-coded subset. Not design-related but coupled to the current metric names.

### 2.7 `visualization/orchestrator.py` — Global metrics list (line 46)

```python
global_metrics = ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']
```

**Impact:** Excludes several available metrics (pLDDT_max, pLDDT_median, contact_prob_max/min/median, pae_max/min/median) from visualization. Not design-related but reduces information content.

---

## 3. Summary of Hard-Coding by Category

### A. Condition identity (severity: HIGH)

| Location | What is hard-coded | Impact |
|----------|-------------------|--------|
| `utils.py` labels dict | 8 POU condition display names | New conditions get generic labels |
| `utils.py` ptm_rank dict | Semantic ordering for 8 POU conditions | New conditions get rank 99 |
| `orchestrator.py` factors fallback | `['factor_DNA', 'factor_PTM']` | Passes nonexistent columns to schema |

### B. Factor inference (severity: HIGH)

| Location | What is hard-coded | Impact |
|----------|-------------------|--------|
| `core_plots.py` F3 factor inference | `_dna` suffix → DNA factor, `pou_` prefix → PTM factor | Only works for POU naming convention |
| `core_plots.py` F3 ptm_pretty | 5 specific PTM state labels | Breaks for any other PTM names |

### C. Metric selection (severity: MEDIUM)

| Location | What is hard-coded | Impact |
|----------|-------------------|--------|
| `orchestrator.py` global_metrics | 5 specific metric column names | Silently excludes 8 available metrics |
| `core_plots.py` F8 f8_metrics | 4 specific metric column names | Ignores schema-driven metric list |

### D. NOT hard-coded (GOOD)

| What | Where | Notes |
|------|-------|-------|
| Condition discovery | Extraction script | Dynamic from filesystem |
| Chain detection | Extraction script + orchestrator | Dynamic from column names |
| Condition count | All plots | Reads from DataFrame |
| Seed count | All plots | Reads from DataFrame |
| Metric count | All plots (except F8) | Reads from schema/columns |

---

## 4. Proposed Minimal Architecture

### 4.1 Core concept: Experiment Metadata File

A single JSON file per experiment that describes:

```json
{
  "experiment_id": "pou_2024",
  "description": "POU protein with PTM and DNA modifications",
  "conditions": {
    "pou_baseline": {
      "label": "Baseline (POU)",
      "attributes": {
        "DNA": false,
        "pTPO101": false,
        "pSEP102": false
      }
    },
    "pou_dna": {
      "label": "POU + DNA",
      "attributes": {
        "DNA": true,
        "pTPO101": false,
        "pSEP102": false
      }
    },
    "pou_tpo101_sep102_dna": {
      "label": "POU (pTPO101 + pSEP102) + DNA",
      "attributes": {
        "DNA": true,
        "pTPO101": true,
        "pSEP102": true
      }
    }
  }
}
```

Key properties:
- **Arbitrary attribute names** — no assumption about what attributes exist
- **Arbitrary condition names** — no assumption about naming convention
- **Arbitrary number of attributes** — not limited to 2 or 3
- **Binary or categorical** — `true`/`false` for binary, strings for categorical
- **No assumption of complete factorial** — missing combinations are simply absent
- **No inference from names** — all structure is explicit

### 4.2 What this replaces

| Current | New |
|---------|-----|
| `ptm_rank` dict in utils.py | Condition ordering from metadata (alphabetical or explicit order) |
| Display labels dict in utils.py | `label` field in metadata |
| Factor inference in F3 | `attributes` field in metadata |
| `factor_DNA`, `factor_PTM` column inference | Attribute columns derived from metadata |
| Hard-coded `ptm_pretty` dict | Attribute labels from metadata |

### 4.3 What this does NOT replace

| Current | Stays |
|---------|-------|
| `condition_id` assignment | Still auto-assigned by extraction |
| `condition_name` | Still comes from filesystem/JSON |
| Metric names | Still from extraction |
| Statistical methodology | Unchanged |
| Pipeline CSV format | Unchanged |

---

## 5. Files Requiring Modification

### Must modify

| File | Change | Reason |
|------|--------|--------|
| `visualization/utils.py` | Remove `ptm_rank` dict, remove POU labels, add metadata-driven ordering | Condition ordering and labels must come from metadata |
| `visualization/core_plots.py` (F3) | Remove POU factor inference, use metadata attributes | Factor inference is entirely POU-specific |
| `visualization/orchestrator.py` | Remove `factor_DNA`/`factor_PTM` fallback, pass metadata to plots | Factor columns never exist in CSVs |

### Should NOT modify

| File | Reason |
|------|--------|
| `af3inputbuilder/scripts/af3_condition_centric_extraction.py` | Extraction is working correctly; metadata is separate |
| `af3_analysis/pipeline.py` | Core data flow is correct; metadata is orthogonal |
| `visualization/core_plots.py` (F1, F2, F4, F5, F6, F7, F8) | These don't use experimental design directly |
| `io/wide_to_long.py` | Data adapter is design-agnostic |
| `schemas/` | Schema system is independent |
| `statistics/` | Statistical modules are design-agnostic |

### May need minor updates

| File | Change | Reason |
|------|--------|--------|
| `visualization/orchestrator.py` | Load metadata file and pass to plot functions | F3 needs metadata |
| `visualization/core_plots.py` (F2, F5, F7) | Use metadata for condition ordering | Currently uses hard-coded fallback |

---

## 6. Files That Should NOT Be Modified

- **Extraction pipeline** — `af3_condition_centric_extraction.py`
- **Core analysis** — `pipeline.py` (seed aggregation, descriptive stats)
- **Data adapters** — `io/wide_to_long.py`
- **Statistics modules** — `statistics/descriptive.py`, `statistics/resampling.py`, etc.
- **Exploratory modules** — `exploratory/variance.py`, `exploratory/distributions.py`, etc.
- **Schema definitions** — `schemas/`
- **Test infrastructure** — `tests/conftest.py`, existing test files

---

## 7. Metadata Format Decision

### Options considered

| Option | Pros | Cons |
|--------|------|------|
| **JSON file** | Simple, human-readable, no dependencies | Manual editing |
| YAML file | More readable for complex structures | Requires pyyaml dependency |
| CSV | Easy to inspect | Poor for nested attributes |
| Python dict in code | No file I/O | Not portable, couples metadata to code |
| Existing `_data.json` | Already exists per condition | Doesn't contain attribute info |

### Recommended: JSON file

Rationale:
1. **No new dependencies** — `json` is in the standard library
2. **Human-readable** — easy to verify and edit
3. **Portable** — works across environments
4. **Consistent** — the project already uses JSON extensively (`*_confidences.json`, `*_data.json`)
5. **Simple** — a single file per experiment, loaded once at pipeline start

### File location

Place the metadata file alongside the data directory:

```
testdata/pou2/
├── experiment_metadata.json    ← NEW
├── pou_baseline/
├── pou_dna/
└── ...
```

This keeps metadata co-located with the data it describes, and the path can be passed to the pipeline as a parameter or auto-discovered.

---

## 8. Design Inspection Component

The metadata system should include a small inspection component that can answer:

- How many conditions are defined?
- What attributes exist?
- What are the levels of each attribute?
- Which attribute combinations are observed?
- Which combinations are absent?
- Is this a complete factorial design?
- What is the attribute structure (binary, categorical, mixed)?

This is **descriptive**, not **inferential**. It describes the design without selecting statistical models.

### Integration point

The inspection component should be callable from:
1. The pipeline (logged to manifest)
2. The visualization layer (to adapt plots)
3. The user (via CLI or function call)

---

## 9. Migration Path

### Phase 1 (this task): Create metadata system + audit
- Create `experiment_metadata.json` for POU dataset
- Implement metadata loading/validation
- Write tests
- Document

### Phase 2 (future): Update visualization
- Replace hard-coded labels/ordering/factors with metadata
- Update F3 to use metadata attributes
- Verify all 8 figures still produce correct output

### Phase 3 (future): Pipeline integration
- Auto-discover metadata file from data directory
- Log metadata to manifest
- Pass to visualization stage

### NOT in scope
- Statistical methodology changes
- Extraction changes
- New visualization types
- Factorial model selection
