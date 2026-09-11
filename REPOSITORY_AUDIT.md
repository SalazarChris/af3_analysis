# Repository Audit

**Date:** 2026-08-26
**Status:** Read-only audit — no files modified.

---

## 1. Directory Structure

```
.
├── AGENTS.md                          # Agent instructions
├── knowledge.md                       # Project knowledge base
├── REPOSITORY_AUDIT.md                # This file
├── .freebuff/project-id               # Freebuff metadata
│
├── af3inputbuilder/                   # AF3 input-building system
│   ├── af3.py                         # Main menu / entry point (1424 lines)
│   ├── .gitignore
│   ├── .gitattributes
│   ├── af3_builder/                   # Core library
│   │   ├── __init__.py                # Public API re-exports
│   │   ├── cli.py                     # Expert CLI menu
│   │   ├── core/                      # Entities, job builder, seeds, references
│   │   ├── ui/                        # Terminal UI (colors, prompts, interactive wizards)
│   │   ├── utils/                     # I/O, templater, MSA, JSON inline, bonds
│   │   ├── validation/                # AF3 JSON schema validator
│   │   ├── fixtures/data/             # Fixture data (empty)
│   │   └── tests/                     # Empty
│   └── scripts/                       # Standalone tool scripts
│       ├── af3_wizard.py              # Beginner-friendly wizard
│       ├── af3_json_validator.py      # CLI JSON validator
│       ├── af3_condition_centric_extraction.py  # Metric extraction
│       ├── add_ions.py                # Ion/ligand sweep generator
│       ├── add_ions_wizard.py         # Interactive ion sweep wizard
│       ├── msa_extractor.py           # MSA extraction from results
│       ├── msa_wizard.py              # Interactive MSA extraction
│       ├── WINDOWS_PYTHON_SETUP_CHEATSHEET.txt  # User-specific setup notes
│       └── af3_metadata/              # Metadata loader/matcher/registry
│
├── af3_analysis/                      # AF3 output-analysis package
│   ├── __init__.py                    # Package init, version 0.1.0
│   ├── __main__.py                    # Module entry: python -m af3_analysis
│   ├── cli.py                         # CLI argument parser
│   ├── config.py                      # AnalysisConfig dataclass, JSON config loader
│   ├── pipeline.py                    # Pipeline orchestrator (6 stages)
│   ├── logging_utils.py               # Dual-output logger (run.log + run.jsonl)
│   ├── py.typed                       # PEP 561 marker
│   ├── io/                            # Data loading, provenance, Parquet store
│   ├── preprocessing/                 # Metadata, canonicalization, QC, aggregation
│   ├── registry/                      # Metric registry and resolution
│   ├── schemas/                       # Enums, records, table schemas, validation
│   ├── statistics/                    # Descriptive, resampling, comparisons, variance, multiplicity
│   ├── structural/                    # Empty
│   ├── reporting/                     # Empty
│   ├── exploratory/                   # Design inventory, variance, distributions, factors
│   ├── visualization/                 # Plotting orchestrator + 8 core plots
│   ├── workflow/                      # RunContext for pipeline state
│   └── tests/                         # Empty (only __init__.py)
│
└── testdata/                          # Development test data
    └── pou2/
        ├── metrics_example/           # Sample confidence JSONs (6 files)
        ├── pou_baseline/              # Full baseline: 10 seeds × 5 samples + CIF + CSVs
        ├── pou_dna/
        ├── pou_sep102/
        ├── pou_sep102_dna/
        ├── pou_tpo101/
        ├── pou_tpo101_dna/
        ├── pou_tpo101_sep102/
        └── pou_tpo101_sep102_dna/
```

---

## 2. Main Menu / Entry Point

**File:** `af3inputbuilder/af3.py` (1424 lines)

This is the single unified entry point. It presents a 6-option menu:

| Option | Action | Mechanism |
|--------|--------|-----------|
| 1 | Job Builder | Imports `af3_wizard.run_wizard` via `importlib` |
| 2 | MSA Extractor | Imports `msa_wizard.main` via `importlib` |
| 3 | Ion/Ligand Sweep | Imports `add_ions_wizard.run_wizard` via `importlib` |
| 4 | JSON Validator | Imports `af3_json_validator.main` via `importlib` |
| 5 | Analysis Pipeline | Calls `run_analysis()` (inline, spawns subprocess) |
| 6 | Complete Analysis Pipeline | Calls `run_complete_analysis_pipeline()` (inline, imports af3_analysis) |

**Classification: KEEP** — but note fragility (see §10).

---

## 3. How the Main Menu Imports the AF3 Builder

`af3.py` uses `sys.path` manipulation to make `af3_builder` importable:

```python
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)      # adds af3inputbuilder/ itself
sys.path.insert(0, _SCRIPTS)   # adds af3inputbuilder/scripts/
```

It imports UI helpers from `af3_builder` with a fallback if the import fails:

```python
try:
    from af3_builder import RESET, BOLD, DIM, ...
except ImportError:
    # inline fallback definitions
```

Scripts in `scripts/` do the same bootstrap from their own `__file__`:

```python
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)      # adds af3inputbuilder/
```

**Classification: NEEDS-REVIEW** — Works but fragile. No package installation or `__init__.py` at the `af3inputbuilder/` level.

---

## 4. AF3 Builder Structure

`af3inputbuilder/af3_builder/` is a well-organized library:

| Subpackage | Purpose |
|------------|---------|
| `core/` | `JobBuilder`, entity classes (Protein/RNA/DNA/Ligand), ID manager, seeds, reference data |
| `ui/` | Terminal UI helpers (`ui.py`) and interactive wizards (`interactive.py`) |
| `utils/` | I/O (JSON load/save/autosave), templater, MSA handling, CIF slicer, bonds |
| `validation/` | AF3 JSON schema validator |

**Public API** (from `__init__.py`): Clean re-exports of ~40 symbols.

**Classification: KEEP** — Well-structured, modular, good separation of concerns.

---

## 5. AF3 Analysis Structure

`af3_analysis/` is a more mature scientific package:

| Subpackage | Purpose | Status |
|------------|---------|--------|
| `io/` | Raw AF3 reader, Phase1 CSV loader, provenance, Parquet store, artifact inventory | Implemented |
| `preprocessing/` | Metadata, canonicalization, QC, aggregation, duplicate resolution, mappings | Implemented |
| `registry/` | Metric registry loading and resolution | Implemented |
| `schemas/` | Enums, frozen records, table schemas with validation, AF3 metric registry | Implemented |
| `statistics/` | Descriptive, resampling (bootstrap/permutation), comparisons, variance decomposition, multiplicity (Holm) | Implemented |
| `exploratory/` | Design inventory, variance components, distributions, factor marginals | Implemented |
| `visualization/` | Orchestrator + 8 core plots (QC, distributions, factorial, forest, variability, trajectories, ECDF, scatter) | Implemented |
| `workflow/` | RunContext for mutable pipeline state | Implemented |
| `structural/` | — | **Empty** |
| `reporting/` | — | **Empty** |
| `tests/` | — | **Empty** (only `__init__.py`) |

**Classification: KEEP** — Scientifically well-designed with proper separation. `structural/` and `reporting/` are stubs.

---

## 6. Major Analysis Pipeline Entry Points

There are **three** different analysis launch paths:

### Path A: Menu option 6 → `run_complete_analysis_pipeline()`
- In `af3.py`, interactively collects parameters
- Imports `af3_analysis.config.create_config_interactive` and `af3_analysis.pipeline.run_pipeline`
- Runs the full 6-stage pipeline in-process

### Path B: Menu option 5 → `run_analysis()`
- In `af3.py`, collects parameters with "smart defaults"
- Spawns a subprocess: `python scripts/af3_analysis.py --models ... --output ...`
- References a script that **does not exist** in the repository (`scripts/af3_analysis.py`)

### Path C: `_run_cli_command()` (not exposed in main menu)
- In `af3.py`, builds a `python -m af3analysis` command
- References `af3analysis` (different name from `af3_analysis`)
- Spawns as subprocess

**Classification:**
- Path A: **KEEP** — functional
- Path B: **DELETE-CANDIDATE** — references nonexistent script
- Path C: **DELETE-CANDIDATE** — references nonexistent package name, not wired into menu

### Pipeline stages (from `pipeline.py`):
1. Extract metrics (loads script from `af3inputbuilder/scripts/af3_condition_centric_extraction.py` via file path)
2. Validate (placeholder)
3. Build manifest
4. Run analysis (seed aggregation, descriptive stats)
5. Generate figures
6. Generate reports (placeholder)

---

## 7. Important Configuration Files

| File | Purpose | Classification |
|------|---------|---------------|
| `af3inputbuilder/.gitignore` | Git ignore rules | **KEEP** |
| `af3inputbuilder/.gitattributes` | Line ending normalization | **KEEP** |
| `knowledge.md` | Project knowledge base | **KEEP** |
| `AGENTS.md` | Agent instructions | **KEEP** |

**Missing configuration files:**
- No `pyproject.toml` or `setup.py` (neither package is installable)
- No `requirements.txt` (dependencies are only listed inline in `af3.py`)
- No `conftest.py` or `pytest.ini`
- No `.flake8`, `ruff.toml`, or `mypy.ini`

**Classification: NEEDS-REVIEW** — The lack of packaging means imports depend entirely on `sys.path` hacks.

---

## 8. Test Files and Test-Data Usage

### Test files
- `af3inputbuilder/af3_builder/tests/` — **Empty directory**
- `af3_analysis/tests/__init__.py` — Only an init file
- `af3_analysis/tests/unit/` — **Empty directory**
- No `test_*.py` or `*_test.py` files anywhere
- No `conftest.py`

### Test data
- `testdata/pou2/metrics_example/` — 6 JSON files (sample confidences, ranking scores)
- `testdata/pou2/pou_baseline/` — Full baseline: 10 seeds × 5 samples with `_data.json`, `_confidences.json`, `_summary_confidences.json`, `_ranking_scores.csv`, `_model.cif` + `TERMS_OF_USE.md`
- `testdata/pou2/pou_*/` — 8 condition directories (pou_dna, pou_tpo101, pou_sep102, etc.)

**Classification:**
- Empty test dirs: **NEEDS-REVIEW** — Should have at least basic smoke tests
- Test data: **KEEP** — Used by extraction pipeline and visualization; has TERMS_OF_USE.md

---

## 9. Important Dependencies

Dependencies are only listed in `af3.py` (inline pip check):

```python
_REQUIRED_PACKAGES = [
    ("numpy",      "numpy>=1.24.0"),
    ("pandas",     "pandas>=2.0.0"),
    ("scipy",      "scipy>=1.10.0"),
    ("matplotlib", "matplotlib>=3.7.0"),
    ("Bio",        "biopython>=1.81"),
    ("tmtools",    "tmtools"),
    ("sklearn",    "scikit-learn"),
]
```

Additional implicit dependencies found in imports:
- `pyarrow` — used in `af3_analysis/io/parquet_store.py`
- `gemmi` — referenced in `af3_builder/utils/cif_slicer.py`
- `af3_selfreference` — imported in `af3.py` line 818 (not in this repo)
- `af3analysis` / `af3_analysis` — the analysis package itself
- `json`, `hashlib`, `re`, `pathlib`, `dataclasses` — stdlib

**Classification:**
- Core deps (numpy, pandas, scipy, matplotlib): **KEEP**
- BioPython, tmtools, scikit-learn: **KEEP**
- pyarrow, gemmi: **NEEDS-REVIEW** — not listed in _REQUIRED_PACKAGES
- `af3_selfreference`: **NEEDS-REVIEW** — external package not in repo; if not installed, menu option 6's self-referential analysis silently fails

---

## 10. Hard-Coded Paths and Fragile Path Handling

### Critical: User-specific hardcoded path

**File:** `af3inputbuilder/af3.py`, lines 161-167:

```python
_AF3_THESIS_PY = r"C:\Users\Chris\.conda\envs\af3_thesis\python.exe"
_ANALYSIS_PYTHON = _AF3_THESIS_PY if os.path.isfile(_AF3_THESIS_PY) else sys.executable
```

Also extensively in `WINDOWS_PYTHON_SETUP_CHEATSHEET.txt` (18 occurrences of `C:\Users\Chris\`).

**Classification: DELETE-CANDIDATE** (the cheatsheet) / **NEEDS-REVIEW** (the path in af3.py; fallback to sys.executable is reasonable but the constant should be configurable).

### sys.path manipulation (13 occurrences across 8 files)

| File | Lines | What it does |
|------|-------|-------------|
| `af3.py` | 27, 32, 36, 119, 800, 1027, 1290 | Adds af3inputbuilder/, scripts/, af3analysis paths |
| `scripts/af3_wizard.py` | 28 | Adds af3inputbuilder/ |
| `scripts/af3_json_validator.py` | 38 | Adds af3inputbuilder/ |
| `scripts/add_ions_wizard.py` | 17, 19 | Adds af3inputbuilder/ and scripts/ |
| `scripts/add_ions.py` | 44 | Adds af3inputbuilder/ |
| `scripts/msa_wizard.py` | 18 | Adds af3inputbuilder/ |

### Incorrect path resolution

**File:** `af3.py`, line 116:
```python
af3analysis_path = os.path.join(af3_toolkit_root, "af3_analysis", "src")
```
This looks for `af3_analysis/src/` which does not exist — the package is at `af3_analysis/` directly. This path fails silently and falls through.

**File:** `pipeline.py`, line 73:
```python
script_path = (
    Path(__file__).resolve().parents[3]
    / "af3inputbuilder"
    / "scripts"
    / "af3_condition_centric_extraction.py"
)
```
Navigates up 3 parent directories from `af3_analysis/pipeline.py`, which assumes a specific directory nesting that could break if the repo is restructured.

**Classification: NEEDS-REVIEW** — All of these are fragile. The `sys.path` approach works but is not robust.

---

## 11. Duplicated, Obsolete, or Suspicious Files

### Duplicated UI helpers

`af3.py` defines its own inline fallback versions of UI functions (`_rule`, `_banner`, `_ok`, `_err`, `_pause`, `_ask_input`, `_ask_number`, `_ask_yn`) that duplicate `af3_builder.ui.ui`. The wizard scripts also define their own fallback versions. Three independent sets of UI fallbacks exist.

**Classification: NEEDS-REVIEW** — Unify to a single fallback location.

### Duplicate analysis paths

Three different code paths attempt to invoke AF3 analysis (see §6). Options 5 and 6 in the menu overlap significantly. Option 5 references a nonexistent script. The `_run_cli_command()` function (not wired to any menu option) references `af3analysis` (wrong package name).

**Classification: DELETE-CANDIDATE** — options 5 and `_run_cli_command()`; **KEEP** option 6 as the primary pipeline.

### Script duplication with library

`scripts/af3_wizard.py` (531 lines) reimplements much of what `af3_builder.cli.py` already provides. The wizard adds a "bus-line" quick-start flow and some UX polish, but shares substantial duplicated logic with the CLI.

**Classification: NEEDS-REVIEW** — Consider whether the wizard can delegate to the CLI rather than reimplementing.

### Suspicious files

| File | Issue |
|------|-------|
| `WINDOWS_PYTHON_SETUP_CHEATSHEET.txt` | User-specific Windows paths; not useful to other developers |
| `af3_builder/fixtures/data/` | Empty directory |
| `af3_builder/tests/` | Empty directory |
| `af3_analysis/structural/` | Empty directory (stub) |
| `af3_analysis/reporting/` | Empty directory (stub) |
| `af3_analysis/tests/unit/` | Empty directory |
| `af3_analysis/py.typed` | PEP 561 marker present but package is not pip-installable |

**Classifications:**
- `WINDOWS_PYTHON_SETUP_CHEATSHEET.txt`: **ARCHIVE** or **DELETE-CANDIDATE**
- Empty directories (`fixtures/data/`, `tests/`, `structural/`, `reporting/`, `unit/`): **NEEDS-REVIEW** — Either populate or remove
- `py.typed`: **NEEDS-REVIEW** — Meaningless without installable packaging

### External package references not in repo

- `af3_selfreference` (imported in `af3.py` line 818) — Not in this repository
- `af3analysis` (referenced in `_run_cli_command`) — Appears to be a stale name; actual package is `af3_analysis`

**Classification: NEEDS-REVIEW** — `af3_selfreference` must be documented as an external dependency or vendored. The `af3analysis` references are stale.

---

## 12. Main Menu Relocation Feasibility

**Current location:** `af3inputbuilder/af3.py`

**Question:** Would moving it to the repository root (`app.py` or `main.py`) be easy or risky?

### What would break

1. **`sys.path` manipulation** — `af3.py` adds its own directory to `sys.path`. Moving to root means `_HERE` changes; the scripts/ bootstrap would need updating.

2. **`_AF3_THESIS_PY` path** — Uses `_HERE` relative path; would still work if moved.

3. **`_SCRIPTS` path** — Currently `os.path.join(_HERE, "scripts")`. After move, this becomes `os.path.join(_HERE, "af3inputbuilder", "scripts")`.

4. **`_launch()` calls** — Uses `importlib.import_module(module_name)` which depends on `sys.path` containing `scripts/`. Would still work if path is adjusted.

5. **`af3_analysis` imports** — Uses `os.path.dirname(_HERE)` to find sibling directory. Would need updating since root is parent of both packages.

6. **`af3_builder` imports** — Uses `from af3_builder import ...` which depends on `sys.path` containing `af3inputbuilder/`. Would still work if path is adjusted.

### Verdict

**RISK: MODERATE**

The relocation is mechanically straightforward (update 3-4 path references) but requires careful testing because:
- There are 7 `sys.path.insert` calls in `af3.py`
- The `_launch()` function dynamically imports scripts by name
- The analysis pipeline path resolution assumes specific directory nesting
- No tests exist to verify nothing breaks

**Recommendation:** Do NOT move until:
1. A proper `pyproject.toml` exists for at least one package
2. Basic smoke tests exist for the main menu
3. The 3 stale analysis paths are cleaned up

---

## 13. Recommended Target Architecture

The smallest sensible next steps, in priority order:

### Immediate (no breaking changes)

1. **Add `pyproject.toml`** at the repo root for `af3_analysis` (it's already a proper package). This eliminates the need for `sys.path` hacks when running `python -m af3_analysis`.

2. **Add `requirements.txt`** at the repo root listing all dependencies (numpy, pandas, scipy, matplotlib, biopython, tmtools, scikit-learn, pyarrow, gemmi).

3. **Remove or archive dead code:**
   - Delete `_run_cli_command()` from `af3.py` (references nonexistent package)
   - Delete menu option 5 (`run_analysis()`) that references nonexistent `scripts/af3_analysis.py`
   - Archive `WINDOWS_PYTHON_SETUP_CHEATSHEET.txt`

4. **Fix the broken path** in `af3.py` line 116: `af3_analysis/src/` should be `af3_analysis/`.

5. **Add a README.md** at the repo root describing how to run the project.

### Short-term (minimal disruption)

6. **Create a thin `af3inputbuilder/__init__.py`** so the scripts can import via package path rather than `sys.path` hacks.

7. **Consolidate UI fallbacks** — Define fallbacks once in `af3_builder/ui/ui.py` and import from there.

8. **Document `af3_selfreference`** as a required external dependency or remove the self-referential analysis menu option.

### Medium-term (structural improvement)

9. **Move main menu to repo root** — After steps 1-8, this becomes safe. Create `app.py` at root that imports from both `af3inputbuilder` and `af3_analysis` as proper packages.

10. **Add basic tests** — At minimum, smoke tests for: (a) `af3_builder.JobBuilder` round-trip, (b) `af3_analysis.config.AnalysisConfig` creation, (c) main menu import.

---

## Summary Classification Table

| Item | Classification |
|------|---------------|
| `af3inputbuilder/af3_builder/` | **KEEP** |
| `af3inputbuilder/af3.py` (main menu) | **KEEP** (with fixes) |
| `af3inputbuilder/scripts/af3_wizard.py` | **KEEP** (deduplicate with cli.py) |
| `af3inputbuilder/scripts/af3_json_validator.py` | **KEEP** |
| `af3inputbuilder/scripts/af3_condition_centric_extraction.py` | **KEEP** |
| `af3inputbuilder/scripts/add_ions.py` | **KEEP** |
| `af3inputbuilder/scripts/add_ions_wizard.py` | **KEEP** |
| `af3inputbuilder/scripts/msa_extractor.py` | **KEEP** |
| `af3inputbuilder/scripts/msa_wizard.py` | **KEEP** |
| `af3inputbuilder/scripts/af3_metadata/` | **KEEP** |
| `af3inputbuilder/scripts/WINDOWS_PYTHON_SETUP_CHEATSHEET.txt` | **ARCHIVE** |
| `af3_analysis/` (entire package) | **KEEP** |
| `af3_analysis/structural/` (empty) | **NEEDS-REVIEW** |
| `af3_analysis/reporting/` (empty) | **NEEDS-REVIEW** |
| `af3_analysis/tests/` (empty) | **NEEDS-REVIEW** |
| `testdata/` | **KEEP** |
| `knowledge.md` | **KEEP** |
| `AGENTS.md` | **KEEP** |
| `af3.py` menu option 5 (analysis pipeline) | **DELETE-CANDIDATE** |
| `af3.py` `_run_cli_command()` | **DELETE-CANDIDATE** |
| `af3.py` `_AF3_THESIS_PY` hardcoded path | **NEEDS-REVIEW** |
| `af3.py` broken `af3_analysis/src` path | **NEEDS-REVIEW** |
| `af3_builder/fixtures/data/` (empty) | **NEEDS-REVIEW** |
| `af3_builder/tests/` (empty) | **NEEDS-REVIEW** |
| `py.typed` (without installable package) | **NEEDS-REVIEW** |
