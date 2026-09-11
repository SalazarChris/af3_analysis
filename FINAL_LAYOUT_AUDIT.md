# Final Layout Audit — FINALVERSIONTHESIS

**Date:** 2026-09-11
**Scope:** verify the expected layout — `af3inputbuilder/` (main scripts, unchanged
except the `af3.py` main menu) plus `af3_analysis/` (3 versions of analysis).
**Method:** read-only inspection (git state, file inventory, imports, test run).
**Changes made by this audit:** none.

> **Update 2026-09-11 — recommendations 1 & 2 executed** (commit `03e1d06`,
> `.gitignore` follow-up `4782f8d`, local only, not pushed):
> - The root duplicate is archived at `.archive/root_package_duplicate/`
>   (kept, not deleted), so `af3_analysis/` is the single canonical package.
> - `af3inputbuilder/` is now tracked by the outer repo — its embedded `.git`
>   was demoted to `.git.nested-backup/` (ignored, reversible).
> - The core `af3_analysis/` package (92 files) is now tracked; 138 total.
> - One bare `from experiment_metadata import` in `core_plots.py` was repointed
>   to `af3_analysis.experiment_metadata`.
>
> Findings 3–10 below remain open.

---

## Verdict

The layout **matches the intent in substance**: `af3inputbuilder/` holds the
input-building toolkit with `af3.py` as the reworked main menu, and
`af3_analysis/` holds the analysis package with **three** analysis versions
(V1, V2, V3), all of which import cleanly and pass their tests.

Two structural deviations remain, both predating this audit:

- **A.** An obsolete **duplicate of the entire core package still sits at the
  repository root** (91 tracked files) and is what the outer git repo actually
  tracks. `af3_analysis/` is the newer, canonical copy (adds V3 + `reporting/`).
- **B.** `af3_analysis/` is only **partially tracked** (46 V3 files), and
  `af3inputbuilder/` is **not tracked at all** by the outer repo because it
  contains its own embedded `.git`.

---

## 1. `af3inputbuilder/` — matches expectation (with caveats)

Own git repository (HEAD `34a2622`), **69 tracked files**, **0 modified tracked files**.
The surviving scripts are unchanged, as expected.

| Item | State |
|---|---|
| `af3.py` | The changed **main menu** — committed (no pending diff). Menu: 1 Job Builder, 2 MSA Extractor, 3 Ion/Ligand Sweep, 4 JSON Validator, 5 Analysis Pipeline → imports `af3_analysis`. |
| `af3_builder/` | Main package: `condition_manifest/` (11 files), `core/` (6), `utils/` (6), `fixtures/` (3), `ui/` (2), `validation/` (1), plus `cli.py` + `__init__.py`. |
| `scripts/` | Wizards + tools: `af3_wizard.py`, `msa_wizard.py`, `add_ions_wizard.py`, `af3_json_validator.py`, `msa_extractor.py`, `af3_condition_centric_extraction.py`, `add_ions.py`. |
| Input JSONs | `pou_segmented_json/` (8), `priority1_json/` (14), `priority2_json/` (12) — untracked. |
| Extras | `smoke_test/`, `scripts/af3_metadata/` — untracked. |

**Uncommitted changes in its nested repo:**

- **27 deletions** (old analysis code removed when analysis moved out):
  `scripts/analysisscripts/*` (23 files), `scripts/af3_analysis.py`,
  `scripts/tests/test_analysisscripts.py`, `scripts/tests/test_integration.py`,
  `itest.txt`.
- **5 untracked directories:** `pou_segmented_json/`, `priority1_json/`,
  `priority2_json/`, `scripts/af3_metadata/`, `smoke_test/`.

**Caveats:** `af3_builder/tests/` is **empty** (0 files). Because of the
embedded `.git`, the outer repo tracks **0** files under `af3inputbuilder/`.

---

## 2. `af3_analysis/` — three analysis versions present

| Version | Location | Entry point | Output directory |
|---|---|---|---|
| **V1** | `visualization/` (`core_plots.py`, `orchestrator.py`, `utils.py`) | `generate_all_figures()` | `<run_dir>/figures/` |
| **V2** | `visualization/v2/` (figure7, figure8, effects, structural, factors, labels, validation, config) | `generate_all_figures_v2()` | `<run_dir>/figures/v2/` |
| **V3** | `visualization_v3/` (20 figures F01–F20 + analysis/ + structural/ layers) | `run_v3_pipeline()` / `python -m af3_analysis.visualization_v3` | `<run_dir>/v3/` |

Shared algorithmic core: `io/`, `preprocessing/`, `registry/`, `schemas/`,
`statistical/`, `structural/`, `exploratory/`, `workflow/`, `pipeline.py`.

**Health checks (all green):**

- All three versions import cleanly (`af3_analysis.visualization`,
  `.visualization.v2`, `.visualization_v3`).
- `python -m af3_analysis.visualization_v3 --help` works.
- **287 tests pass in 396s (6:36)** — no failures, no regressions.

Output directories on disk show all three generations, e.g.
`test_struct_fig/{figures/, figures/v2/, v3/{figures,tables,cache,report,metadata}/}`.

---

## 3. Findings

| # | Severity | Finding |
|---|---|---|
| 1 | **High** | **Root-level duplicate package.** 91 tracked files at the repo root — `io/`, `structural/`, `statistical/`, `visualization/`, `preprocessing/`, `registry/`, `schemas/`, `exploratory/`, `workflow/`, `tests/`, `cli.py`, `config.py`, `pipeline.py`, `experiment_metadata.py`, `logging_utils.py`, `__init__.py`, `__main__.py`, `py.typed`, `.gitignore` — are **byte-identical to `af3_analysis/` modulo CRLF** (verified file-by-file). This is the pre-V3 location: root `visualization/` has V1+V2 but no V3, and root `tests/` has 11 files with no V3 tests. The outer git repo tracks **this** copy. |
| 2 | **High** | `af3_analysis/` is only **46/… tracked**; the entire core package under it is untracked. Version control therefore does not currently follow the canonical code. |
| 3 | **Med** | **Two embedded git repos** prevent outer-repo tracking: `af3inputbuilder/.git` (active, 69 files) and `af3_analysis/.git.nested-backup/` (leftover from a previous session). |
| 4 | **Med** | `af3inputbuilder` nested repo is **dirty**: 27 unstaged deletions + 5 untracked directories. |
| 5 | **Med** | `af3_analysis/reporting/` is **empty** (no files, no `__init__.py`) yet it appears in `af3_analysis/__init__.py`'s `__all__` and in V3's `ARCHITECTURE.md`. It imports only as an implicit namespace package. |
| 6 | **Low** | `af3_analysis/__init__.py` `__all__` is **stale**: lists `errors` (no such module), `statistics` (the directory is `statistical/`) and `reporting` (empty); omits `visualization_v3`. Cosmetic — `__all__` is not executed. |
| 7 | **Low** | **V3 is not wired into the pipeline orchestrator.** `pipeline.py` supports `visualization_version ∈ {v1, v2, both}`; V3 is standalone-only. The "3 versions" exist but only two are driven by `run_pipeline`. |
| 8 | **Low** | Empty placeholder dirs: `af3_builder/tests/`, `af3_analysis/tests/unit/` (only `__init__.py`). |
| 9 | **Low** | Stray artifacts: `tmprb03f3rv.cif` (891 B) at repo root, `.pytest_cache/`, `__pycache__/`, `debug_struct_test/`. |
| 10 | **Info** | Test suite is slow (~6:36) due to real-data fixtures — `test_visualization.py` (~90 s), `test_wide_to_long.py` (~94 s), `test_analysis.py` (~57 s), `test_extraction.py` (~57 s). Not a hang. |

---

## 4. Generated data present (not source)

`run_20260828_001617/`, `run_20260828_002706/`, `run_20260830_214934/`,
`af3_analysis_output_struct_test/`, `af3_analysis_output_struct_fig_test/`,
`debug_struct_test/`, `testdata/` (contains `pou2/` raw AF3 CIFs).

These are pipeline outputs referenced by tests; they are not part of the two
source components.

---

## 5. Recommended next steps (not executed — audit is read-only)

Ordered by leverage:

1. **Decide the canonical package location** (`af3_analysis/`) and **archive the
   root duplicate** rather than deleting it (per repository rules). This removes
   the two-copy ambiguity and is the single highest-value fix.
2. **Bring `af3inputbuilder/` and the core `af3_analysis/` under outer-repo
   tracking** — resolve the embedded `.git` repos (the nested backup is already
   movable/deletable; the active one needs a decision: subtree, submodule, or
   demote to a plain directory).
3. **Commit the 27 `af3inputbuilder` deletions** so its repo is clean and the
   analysis/builder split is recorded.
4. **Resolve `reporting/`** — either implement it or drop it from `__all__` and
   the V3 architecture doc; fix the stale `errors`/`statistics` entries.
5. **Optionally expose V3 via `pipeline.py`** (`visualization_version` gains
   `"v3"`/`"all"`) so one orchestrator drives all three versions.
6. **Clean stray artifacts** (`tmprb03f3rv.cif`, empty dirs) once classified safe.
