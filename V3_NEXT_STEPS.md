# V3 Pipeline — Status & Next Steps

Status: **2026-09-06** — V3 now runs end-to-end on real data
(17/20 figures pass, 6/6 tables written, 0 failed).

## What V3 is

`af3_analysis/visualization_v3/` is an additive visualization layer
consuming existing pipeline run outputs (`<run_dir>/tables/`) plus raw
AF3 outputs (CIF + confidences JSON). It writes **only** to
`<run_dir>/v3/` and never modifies existing pipeline outputs. It is
invoked with:

```bash
python -m af3_analysis.visualization_v3 <run_dir> \
    --raw-af3-root <raw_af3_root> \
    [--metadata experiment_metadata.json] [--figures F01,F02,...] [--force]
```

Design and module map: `af3_analysis/visualization_v3/ARCHITECTURE.md`.

## Audit findings (what was broken and is now fixed)

All fixes below were made 2026-08-30 … 2026-09-06. Each is covered by a
unit test (`af3_analysis/tests/test_v3_runner.py`,
`test_v3_adapter.py`) and validated end-to-end on
`af3_analysis_output_struct_fig_test/test_struct_fig` + `testdata/pou2`.

| # | Bug | Fix |
|---|-----|-----|
| 1 | `figures/__init__.py` never exported `f11_domain_motion`; `runner.py` imported it → `ImportError` on any import | export added |
| 2 | `validation.py` used `@dataclass` without importing it → `NameError` | import added (plus mutable-default fix in `ValidationResult`) |
| 3 | `structural/comparability.py` + `contacts.py` same missing `dataclass` import | imports added |
| 4 | `f12_structural_clustering.py` had a non-default parameter after a default one → `SyntaxError` | reordered |
| 5 | Runner wired data only for F01–F03; F04–F20 returned `not yet implemented` | all 20 figures wired with correct data shapes |
| 6 | `_generate_v3_tables` was a stub (`pass`) — no CSVs written | all 6 tables implemented (structural_summary, per_residue_displacement, contact_changes, interface_contacts, structural_clusters, confidence_geometry) |
| 7 | No CLI entry point | `__main__.py` added |
| 8 | No V3 tests | 16 V3 tests added |
| 9 | `adapter._convert_to_structure_data` used `structure.prediction_id`, which does not exist on `NormalisedStructure` → `AttributeError` on every real CIF | prediction_id built from condition/seed/sample (project convention); real CIF smoke test added |
| 10 | Adapter keyed predictions by CIF stem (`pou_baseline`) but conditions/seeds by run `condition_id` (`cond_001`) → matched-seed pairing silently empty | `_build_condition_stem_map` bridges via `condition_registry.csv` replicate_ids + `seed_aggregated.condition_name`; unmapped stems logged as warnings |
| 11 | Adapter left pLDDT/PAE/contact_prob `None` → F15/F20 dead | `_load_confidence_metrics` reads the sibling `_confidences.json` (per-prediction; missingness preserved, no imputation) |
| 12 | `hierarchical_clustering` condensed the distance matrix with `squareform`, but `AgglomerativeClustering(metric="precomputed")` needs the full (N,N) matrix → `ValueError: Expected 2D array` | pass the full matrix |
| 13 | F01's generator has no `reference_condition` parameter; runner passed it via `**common_params` → `TypeError` | F01 called with only its accepted params |
| 14 | `V3Config.reference` defaults to `{}` (not `None`), so `resolve_reference` never used the "first condition" fallback → reference resolution returned `None` | empty dict treated as no explicit reference |
| 15 | F15 used `struct.has_dna`, which `StructureData` lacked | `has_dna` property added (derived from polymer types) |

## Validation results (real data)

Run: `af3_analysis_output_struct_fig_test/test_struct_fig`
Raw AF3: `testdata/pou2` (8 conditions × 50 predictions = 400 CIFs)

```
Figures: 17 pass, 3 skipped (F10 sites, F11 regions, F19 metadata — all config-driven), 0 failed
Tables:  6/6 written
Pairwise: 400 structures, mean RMSD 9.67 Å, median 11.13 Å
```

Sample observations per figure: F01 n=400, F02 n=350, F03 n=70,
F04/F05 n=53550, F06 n=4M, F12–F14 n=400, F15/F20 n=350,
F16–F18 n=7 (see note below).

**Note on F16/F17/F18 n=7:** in this particular test run,
`seed_aggregated.csv` contains only **1 seed for cond_001** (10 for
others), while its CIFs carry 10 seeds. Seed-level confidence pairing
therefore has only 7 (condition, seed) rows with both deltas present.
This is faithful missingness handling, not a bug. Runs whose seed tables
are complete will pair all seeds.

## Known limitations / remaining work

1. **F10/F11 need configured sites/regions.** They are intentionally
   config-driven (`v3_config.sites`, `v3_config.regions`). For the POU
   thesis, the DNA-binding domain (e.g. POU-specific + POU-homeodomain)
   and any post-translational-modification sites must be defined in a V3
   config JSON and passed with `--config`. *Scientific decision required:
   which residues define the domains/sites.*

2. **F19 needs experiment metadata.** `--metadata
   testdata/pou2/experiment_metadata.json` enables factor contrasts
   (DNA × PTM). The metadata file exists — rerun with `--metadata`.

3. **Reference condition default is alphabetical first**
   (`cond_001` = pou_baseline). For the thesis, verify this matches the
   intended reference (baseline POU) or pass `--reference`.

4. **`contact_prob_mean` from confidences JSON** uses the off-diagonal
   upper triangle of the (tokens × tokens) matrix. Confirm this matches
   the V2 pipeline's definition (which derived it from
   `contact_probs` similarly) before using in the thesis.

5. **Structural figures are compute-heavy**: F06 (contact map
   difference) touches ~4M residue pairs for 400 predictions; the full
   run takes ~5 min (largely pairwise matrix + per-residue work).
   Pairwise matrix is cached in `<run_dir>/v3/cache/`; the per-prediction
   loops (F04/F05/F06/F08/F09/F15/F20 + tables) are not cached and
   recompute displacement/contacts on each run.

6. **No integration test with full seed tables.** The unit tests use
   synthetic structures; the end-to-end validation used a test run whose
   reference condition has a single table seed. A validation run on a
   production run dir (e.g. `run_20260830_214934`) with complete seed
   tables and raw AF3 root present is still pending.

7. **`_build_seed_summaries` collapses samples** (n_samples=1) — the V3
   seed summaries are built from `seed_aggregated.csv`, which is already
   sample-aggregated by the pipeline. Fine for seed-level analyses;
   prediction-level analyses read `dataset.predictions` directly.

## Suggested next steps (in order)

1. **Configure sites/regions + metadata, re-run on the real run:**
   create a `v3_config.json` with the POU DNA-binding domains
   (`regions`) and any PTM sites (`sites`), then run with
   `--metadata testdata/pou2/experiment_metadata.json` on a full run
   dir. (Needs the scientific definition of the domains.)
2. **Validate F16–F18 seed pairing on a run with complete seed
   tables** (`run_20260830_214934` when its raw AF3 root is available).
3. **Cache per-prediction structural results** (displacement, contacts,
   interfaces) keyed like the pairwise matrix to cut the ~5 min runtime.
4. **Cross-check `contact_prob_mean` and PAE definitions** against the
   V2 metric pipeline (`METRICS_DESCRIPTION.md`) to guarantee identical
   metric semantics before thesis figures are generated.
5. **Write a short scientific-methods section** describing V3 outputs
   (unit of analysis per figure, matched-sample vs matched-seed pairing)
   for the thesis appendix — most of the material is already in
   `ARCHITECTURE.md` and the runner docstring.