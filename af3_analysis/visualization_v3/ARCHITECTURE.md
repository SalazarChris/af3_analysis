# V3 Visualization Pipeline — Architecture Map

**Status:** additive extension. V3 does not modify, rewrite, or replace any existing
module. The exception clause (V3 cannot function without a change) was never triggered:
V3 consumes existing functionality exclusively through public imports.

---

## Data flow

```
EXISTING PIPELINE (unchanged)
  ├─ extraction (af3inputbuilder/scripts/af3_condition_centric_extraction.py)
  ├─ pipeline.py stages 1-6  →  run_YYYYMMDD_HHMMSS/{tables,figures,manifest,reports,logs}
  ├─ structural Stage 4b     →  tables/structural_*.csv  (scalar metrics only)
  └─ raw AF3 output root     →  condition/seed-*/..._model.cif  + experiment_metadata.json
                    │
                    ▼
            ┌───────────────┐
            │  V3 ADAPTER    │  adapter.py   (reads; never writes outside run_dir/v3/)
            └──────┬────────┘
                   ▼
        V3 NORMALIZED DATA MODEL  (model.py)
        Dataset → Condition → Seed → Prediction → StructureData → Residue/Atom classes
                   │
     ┌─────────────┼──────────────────┐
     ▼             ▼                  ▼
 CONFIDENCE    GEOMETRY           INTERFACES
 (tables CSVs  (structural/       (structural/
 + per-residue   rmsd, displacement,  interfaces.py)
 pLDDT)          contacts, regions,
                 pairwise, clustering)
     └─────────────┼──────────────────┘
                   ▼
        INTEGRATED V3 ANALYSIS  (analysis/)
        matched-seed aggregation, confidence x geometry,
        effect sizes (reuses statistical/comparisons.py),
        factorial contrasts (metadata-driven), seed reproducibility
                   │
                   ▼
        V3 FIGURE GENERATION  (figures/, V3-F01 .. V3-F20)
                   │
        ┌──────────┼──────────────┐
        ▼          ▼              ▼
     FIGURES     TABLES         REPORT
     (run_dir/v3/figures/)  (run_dir/v3/tables/)  (run_dir/v3/report/V3_VISUALIZATION_REPORT.md)
                             + figure_manifest.json (run_dir/v3/metadata/)
```

## Module map

| Module | Responsibility | Phases |
|---|---|---|
| `config.py` | V3Config dataclass, JSON load/validate, figure switches, thresholds | 5 |
| `model.py` | Normalized hierarchy + compact StructureData + QC records | 6 |
| `adapter.py` | Builds Dataset from run tables + raw root (reuses `discover_structures`, `parse_mmcif`, `load_experiment_design`) | 6 |
| `validation.py` | Condition/seed/prediction mapping validation, row classification (summary rows), structural QC, coverage | 7, 8 |
| `reference.py` | Explicit reference resolution from configuration only | 9 |
| `structural/comparability.py` | Common structural space (chains, residues, atoms), coverage, LOW_COVERAGE flagging | 10 |
| `structural/rmsd.py` | Global Cα RMSD (+ configurable atom), coverage; reuses `alignment._kabsch` | 11 |
| `analysis/matched_seed.py` | Matched-seed (and matched-sample) pairing, seed-level aggregation | 12, 28 |
| `structural/displacement.py` | Per-residue displacement vectors, magnitudes, seed direction consistency | 13, 14 |
| `structural/contacts.py` | Residue contact maps, Δcontact encoding (+1/0/-1), gain/loss summaries | 15, 16 |
| `structural/interfaces.py` | Entity-aware interface contacts (protein-DNA/RNA/ligand/ion/protein) | 17, 18 |
| `structural/regions.py` | Config-driven local sites/regions: local RMSD, local displacement, local confidence | 19, 20 |
| `structural/pairwise.py` | All-prediction pairwise RMSD matrix (cached) | 21 |
| `structural/clustering.py` | Hierarchical clustering + MDS embedding of predicted structural space | 22, 23 |
| `analysis/confidence.py` | Confidence x geometry integration (kept conceptually separate) | 24, 25 |
| `analysis/effects.py` | Structural effect sizes via existing `statistical.comparisons` framework | 26 |
| `analysis/factorial.py` | Metadata-driven one-attribute factorial contrasts | 27 |
| `figures/` | V3-F01..F20 generators, failure-safe runner, pathological-figure guards | 29, 30, 35, 36 |
| `reporting/manifest.py` | `figure_manifest.json` (SUCCESS/SKIPPED/FAILED, never hides) | 31 |
| `reporting/tables.py` | structural_summary.csv, per_residue_displacement.csv, contact_changes.csv, interface_contacts.csv, structural_clusters.csv, confidence_geometry.csv | 33 |
| `reporting/report.py` | `V3_VISUALIZATION_REPORT.md` | 43 |
| `cache.py` | Disk cache for pairwise matrices / contact maps, config-hash keys | 37 |
| `runner.py` | Orchestration with `[V3]` stage logging | 34 |

## Protected surface (must remain byte-identical)

Everything outside `af3_analysis/visualization_v3/`. In particular:

- `af3_analysis/pipeline.py`, `af3_analysis/config.py`, `af3_analysis/cli.py`
- `af3_analysis/structural/**` (imported, never modified)
- `af3_analysis/statistical/**` (effect sizes reused via `compare_two_conditions`)
- `af3_analysis/visualization/**` (V1 + V2 figures untouched)
- `af3_analysis/experiment_metadata.py`, `af3_analysis/io/**`, `af3_analysis/schemas/**`
- `af3inputbuilder/**`
- Existing run directories (V3 writes only under `<run_dir>/v3/`)

Verification: `git status --porcelain` (tracked tree must stay clean), md5 checksums of
`run_20260830_214934/tables/*.csv` recorded before V3 (in `/tmp/v3_baseline_checksums.txt`),
and re-running the pre-V3 existing test baseline (120 tests) after V3.

## Unit-of-analysis rules

| Level | Definition | Used by |
|---|---|---|
| prediction | one (condition, seed, sample) CIF | pairwise matrix, contact maps, displacement vectors |
| seed | aggregation over samples of one (condition, seed) | matched-seed comparisons, effect sizes |
| condition | aggregation over seeds | summaries, factorial contrasts |

Seeds are prediction-process robustness samples (not biological replicates, not a
physical ensemble). Sample indices are treated as exchangeable draws; the default
`matched_sample` pairing pairs identical (seed, sample) indices and is documented as
arbitrary-but-deterministic; `seed_pool` pairing compares all sample combinations.
