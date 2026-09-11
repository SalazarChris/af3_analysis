# AF3 Analysis — Data-Handling Audit

**Date:** 2026-09-11
**Scope:** `af3_analysis/` — how the code ingests, transforms, keys, aggregates, and
reports AF3 data, and where that diverges from the documented scientific contract.
**Method:** read-only source inspection + evidence from real run outputs under
`run_20260830_214934/`. **No code was changed by this audit.**

> **How to use this document.** It is written for a downstream reviewer/LLM that will
> propose improvements. Section 7 is a **defect register** with `file:line` evidence and
> severity — that is the actionable core. Section 8 lists places where *the same concept
> has multiple conflicting implementations*, which is the single biggest source of
> non-reproducibility. Section 11 lists **scientific decisions that must not be guessed**
> (per `AGENTS.md`). Do not "fix" anything in Section 11 without an explicit decision from
> the thesis author.

---

## 1. Executive summary

`af3_analysis` is presented as a "registry-driven scientific analysis package". In reality
it contains **two very different layers**:

1. **A thin, working pipeline** (`pipeline.py` → `af3inputbuilder` extraction script →
   ad-hoc pandas aggregation → V1/V2 figures). This is what actually runs and produces
   the tables in `run_*/tables/`.
2. **A large, carefully-designed analysis layer** (`preprocessing/`, `registry/`,
   `schemas/`, `statistical/`, `exploratory/`, `workflow/`) that implements definedness
   taxonomy, QC gates, resampling, Holm correction, and effect sizes — but which
   **nothing imports**. It is dead code with respect to the pipeline.

The consequences:

- The `canonical_analysis` mode's headline outputs are computed by ~40 lines of ad-hoc
  pandas in `pipeline.py`, **not** by the apparatus built to do it.
- `pairwise_comparisons.csv` — a required output table for effect-size figures — is
  **written empty** (`pipeline.py:233`). Every "effect size" figure is therefore run on
  an empty frame.
- Structural comparisons are computed as a **sample×sample cross-product** rather than
  matched-sample, inflating comparison counts ~5× and violating independence.
- `comparison_id` is **not unique** (69,000 rows → 230 distinct IDs in the real run).
- Two **different, conflicting** implementations of descriptive/variability statistics
  exist, plus a third ad-hoc one in `pipeline.py`.
- At least one module (`statistical/factorial_models.py`) **does not compile**.

---

## 2. Environment hazard: which `af3_analysis` is even executing?

There are **two copies** of the package, and an editable install points at the *older* one.

| Copy | Path | `.py` files | Notes |
|---|---|---|---|
| **Local (authoritative)** | `C:\Users\Chris\Desktop\FINALVERSIONTHESIS\af3_analysis` | 134 | Has `structural/`, `io/structure_reader.py`, `visualization_v3/`, full tests |
| **Editable-installed (stale)** | `C:\Users\Chris\Desktop\Thesis Project\af3-toolkit\af3_analysis\src\af3_analysis` | 68 | Has a `statistics/` dir (old name), no `structural/` |

The venv contains **two** editable installs:

```
site-packages/__editable__.af3_analysis-0.1.0.pth
site-packages/__editable__.af3analysis-0.1.0.pth
site-packages/__editable___af3analysis_0_1_0_finder.py
```

`__editable__.af3_analysis-0.1.0.pth` contains:
`C:\Users\Chris\Desktop\Thesis Project\af3-toolkit\af3_analysis\src`

**Resolution depends on cwd.** Running from the repo root, the local copy wins (cwd
precedes site-packages on `sys.path`). Running from any other directory — including from
inside `af3_analysis/` — silently loads the **stale** copy. This was observed directly:
`python -m py_compile statistical/*.py` from `af3_analysis/` produced a traceback pointing
at `...\Thesis Project\af3-toolkit\af3_analysis\src\af3_analysis\...`.

This also explains the stale `__all__` in `af3_analysis/__init__.py` (§7, D-22): the
installed copy still has `statistics/`.

**Recommendation seed:** uninstall both editable installs, or re-point them at the
canonical copy, and add a test that asserts `af3_analysis.__file__` is inside the repo.

---

## 3. Package map: wired vs. orphaned

Import counts exclude tests and self-imports.

| Subpackage | Importers outside itself | Status |
|---|---|---|
| `io/` | 4 | partially used (`structure_reader` only; see below) |
| `schemas/` | 8 | used as a type/enum vocabulary |
| `structural/` | 16 | **used** (structural stage 4b) |
| `visualization/` | 21 | **used** (V1 + V2) |
| `visualization_v3/` | 0 | used, but only via `python -m` (expected) |
| `experiment_metadata.py` | 3 | **used** |
| `preprocessing/` | **0** | **orphaned** |
| `registry/` | **0** | **orphaned** |
| `statistical/` | **0** | **orphaned** |
| `exploratory/` | **0** | **orphaned** |
| `workflow/` | **0** | **orphaned** |
| `logging_utils.py` | **0** | **orphaned** |

Within `io/`, only `structure_reader` is actually consumed (14 references). `raw_af3_reader`
and `provenance` have 1 reference each (a test); **`phase1_loader`, `parquet_store`, and
`artifact_inventory` have zero references anywhere** — nothing in the repository calls them.

**Implication:** roughly half the package's ~32k lines are unreachable from any entry point
that the thesis actually runs.

Concretely, these public API functions are defined and exported but **called from nowhere**
in the repository (verified by exhaustive grep, excluding their own `def`/`__all__`):

- `statistical.comparisons.compare_two_conditions`
- `statistical.multiplicity.holm_correction`
- `statistical.factorial_models.run_factorial_analysis`
- `preprocessing.aggregation.aggregate_to_seed_level`
- `preprocessing.quality_control.run_quality_control`
- `preprocessing.canonicalize.canonicalize_artifacts`
- `statistical.eligibility.get_analysis_eligibility`
- `MultiConditionResult` is **never constructed** anywhere

So there is currently **no code path that performs a statistical test with multiplicity
control** — the machinery exists and is exported, but nothing invokes it.

---

## 4. End-to-end data flow (as actually executed)

```
raw AF3 output root
      │
      │  af3inputbuilder/scripts/af3_condition_centric_extraction.py
      │  (loaded by importlib from a file path — see D-16)
      ▼
<tables>/metrics_replicates.csv     one row per (condition, replicate); wide metric columns
<tables>/metrics_conditions.csv     per-condition descriptive stats (snake_case names)
<tables>/condition_manifest.csv, condition_registry.csv
      │
      │  pipeline._stage_run_analysis()            ← ALL of this is ad-hoc pandas
      │    · extract_seed() regex on replicate_id
      │    · select_dtypes([np.number])  ← sweeps in non-metric columns (D-01)
      │    · groupby([condition_id, condition_name, seed]).mean()
      ▼
<tables>/seed_aggregated.csv        one row per (condition, seed)
<tables>/descriptive_stats.csv      mean/std of seed_aggregated by condition
<tables>/pairwise_comparisons.csv   EMPTY (header only)  ← D-02
      │
      ├─► visualization/orchestrator.generate_all_figures()      V1 → <figures>/
      ├─► visualization/orchestrator.generate_all_figures_v2()   V2 → <figures>/v2/
      └─► (V3 run separately: python -m af3_analysis.visualization_v3 → <run>/v3/)

if config.coordinate_analysis_enabled:
      │  structural/discovery → io/structure_reader.parse_mmcif
      │  → structural/alignment (Kabsch) → structural/comparison (metrics)
      ▼
<tables>/structural_predictions.csv, structural_metrics.csv,
         structural_comparisons.csv, structural_effects.csv
```

**Structurally important:** the actual metric extraction lives in the **`af3inputbuilder`
repository**, not in `af3_analysis`. `af3_analysis` only consumes its CSVs. So the
authoritative definition of what `pLDDT_mean`, `pae_mean`, `contact_prob_mean`, etc. *mean*
is outside this package. Any audit of metric semantics must include
`af3inputbuilder/scripts/af3_condition_centric_extraction.py`.

---

## 5. Data contracts and naming

### 5.1 Table schemas (observed, `run_20260830_214934`)

| Table | Rows | Grain | Key columns |
|---|---:|---|---|
| `metrics_replicates.csv` | 2,448 | (condition, replicate) | `condition_id, condition_name, replicate_id` |
| `metrics_conditions.csv` | 24 | condition | `condition_id, n_replicates, n_seeds, <metric>_mean/_sd/_median/_min/_max/_cv` |
| `seed_aggregated.csv` | 240 | (condition, seed) | `condition_id, condition_name, seed` |
| `descriptive_stats.csv` | 24 | condition | `condition_id, <metric>_mean, <metric>_std` |
| `pairwise_comparisons.csv` | **0** | — | `reference, condition, metric, diff_mean` (header only) |
| `structural_predictions.csv` | 1,200 | prediction | `prediction_id, condition_id, seed, sample` |
| `structural_metrics.csv` | 18,400 | (prediction, metric) | `prediction_id, metric_id, scope_type, scope_id, value, status, reason` |
| `structural_comparisons.csv` | 69,000 | (comparison, metric) | `comparison_id, prediction_a_id, prediction_b_id, metric_id` |
| `structural_effects.csv` | 207 | (metric, condition pair) | `metric_id, condition_a, condition_b, n, mean, median, sd` |

### 5.2 Naming is inconsistent between tables for the same quantities

The **same physical quantity** is named differently depending on which table you read:

| Quantity | `seed_aggregated.csv` | `metrics_conditions.csv` |
|---|---|---|
| Residues in chain A | `chain_A_residues` | `chain_A_residue_count_mean` |
| Mean pLDDT chain A | `chain_A_plddt` | `chain_A_avg_plddt_mean` |

And within a single table there are **two competing conventions** for per-condition
summaries: `metrics_conditions.csv` uses `<metric>_mean/_sd/_median/_min/_max/_cv`, while
`descriptive_stats.csv` uses `<metric>_mean/_std`. Any code that joins or reuses these
must special-case each. This is a hazard for provenance ("every reported numerical result
must have traceable provenance").

### 5.3 Redundant metric variants

`metrics_replicates.csv` and `seed_aggregated.csv` carry `pLDDT_mean`, `pLDDT_max`,
`pLDDT_min`, `pLDDT_median` (and the same 4-way split for `pae` and `contact_prob`).
These are **not independent measurements** — they are order statistics of one distribution.
Treating them as a metric family inflates any multiplicity correction (e.g. Holm over 15
near-duplicate "metrics") without adding information. `METRICS_DESCRIPTION.md` should be
the arbiter of which variants are scientifically distinct.

---

## 6. Unit of analysis, missingness, and exclusions

### 6.1 Intended semantics (documented)

`schemas/enums.py` defines an explicit **definedness taxonomy**, which is the correct
design for this domain:

| `Definedness` | Meaning |
|---|---|
| `present` | valid numeric value exists |
| `undefined_by_composition` | cannot exist for this condition (e.g. interface metric in a monomer) |
| `missing_technical` | should exist but extraction failed |
| `not_collected` | Phase-1 schema never extracted it |

`AnalysisLevel` defines `sample` vs `seed`, and `Aggregator`'s docstring states:
*"Samples are nested within seeds. Default inferential analyses use one seed-level
observation per condition × metric × scope."* — i.e. **seed is the unit of analysis**.

### 6.2 Actual behaviour

- `pipeline.py` **never emits a `definedness` column**. The taxonomy is not carried into
  any production table.
- `pipeline.py` aggregates samples → seed with a plain `mean()` (D-01), implicitly treating
  `undefined_by_composition` and `missing_technical` as *the same kind of absent* — i.e. it
  **conflates structural missingness with technical missingness**, which `AGENTS.md`
  explicitly forbids.
- The one place `definedness` *is* computed, `preprocessing/aggregation.py:114/137/152`,
  assigns the whole seed the definedness of the **first sample** (`group.iloc[0]`), so a
  seed with mixed sample-level definedness is mislabelled.
- `structural/` comparisons propagate `alignment_status` / `metric_status`
  (`present`, `not_applicable`, `undefined`, `computation_error`,
  `skipped_alignment_failed`) — this is better, and is the only place missingness is
  genuinely preserved end-to-end.

### 6.3 Exclusions are not accounted for

- `preprocessing/quality_control.py` builds an `_exclusion_log`, but only
  `_file_provenance_checks` appends to it (line 119). Duplicate-key, invalid-seed,
  invalid-metric and `PRESENT_WITHOUT_VALUE` findings are recorded as findings with
  `action="exclude"` but **never enter the exclusion log** — yet `excluded_records` is
  computed as `len(self._exclusion_log)` (line 91). So exclusion counts are systematically
  under-reported.
- `preprocessing/canonicalize.py:90` returns `excluded_count=0` with an explicit
  `# TODO: track excluded records`.
- QC is orphaned (§3), so **no exclusion is actually applied** anywhere in the live pipeline.

---

## 7. Defect register

Severity: **S1** = produces wrong/unusable scientific output or blocks execution;
**S2** = correctness/provenance risk; **S3** = robustness/maintainability.

### S1 — Wrong or unusable output

| ID | Location | Defect |
|---|---|---|
| **D-01** | `pipeline.py:211,219,224` | Sample→seed aggregation uses `df.select_dtypes(include=[np.number])` and `groupby(...).mean()`. This sweeps **non-metric numeric columns** into the aggregation. Confirmed in production: `descriptive_stats.csv` contains `chain_A_residues_mean`/`_std` — statistics of a *fixed sequence length*. `io/wide_to_long.py` explicitly excludes `chain_<X>_residues` as "fixed sequence lengths, not measured confidence values", so the pipeline contradicts the package's own contract. |
| **D-02** | `pipeline.py:233` | `pairwise_comparisons.csv` is written as an **empty** DataFrame literal. No pairwise comparison is ever computed by the live pipeline, yet `plot_effect_size_forest` consumes this file. All effect-size output is vacuous. |
| **D-03** | `statistical/factorial_models.py:33` | **SyntaxError** — the file contains literal `\"\"\"` escapes inside the class body. `python -m py_compile` fails; the module cannot be imported. Its content is also non-functional scaffolding (placeholder `Resampler`, `run_model_fit` returns a stub, `test_estimability` always returns `is_estimable: True`, hardcoded `factors = ["ResidueType", "ChainPairing"]`, `grouping_key="SeedID"`). Treat as **unimplemented, not merely buggy**. |
| **D-04** | `structural/comparison.py:84,176` | `_make_comparison_id()` omits **sample**, so IDs collide across samples. Measured: `structural_comparisons.csv` has **69,000 rows but only 230 unique `comparison_id`** (68,770 duplicates). The table also has no `sample_a`/`sample_b` columns — sample identity survives only inside `prediction_*_id` strings. Any `groupby("comparison_id")` aggregates across samples. |
| **D-05** | `structural/comparison.py:170-180` | `matched_seed` mode nests `for target_struct in group.predictions[seed]` × `for ref_struct in ref_group.predictions[seed]` — a **sample×sample cross-product**, not matched-sample pairing. With 5 samples this yields 25 comparisons per (condition, seed) instead of 5 (→ 5,750 instead of 1,150), and reuses each sample many times, violating independence. |
| **D-06** | `statistical/comparisons.py:84-85` | For **unpaired** comparisons the bootstrap CI is computed on `values_a` alone, so `native_ci_lower/upper` describe the **mean of group A**, not the A−B difference. The CI is attached to a difference statistic and is therefore misleading. (Paired path is correct.) |
| **D-07** | `statistical/descriptive.py:109,112` | `between_seed_var` and `between_cond_var` are computed by the **same expression** (`seed_means.var(ddof=1)`), so "between-condition variance" is a duplicate of "between-seed variance" under a different name. |

### S2 — Correctness / provenance risk

| ID | Location | Defect |
|---|---|---|
| **D-08** | `pipeline.py:151` | `_stage_validate` is a **placeholder** that unconditionally returns `status="pass"`. No validation ever runs. |
| **D-09** | `pipeline.py:321` | `_stage_generate_reports` is a **placeholder**; the `reports/` directory is never populated. |
| **D-10** | `pipeline.py:118-151` vs `preprocessing/aggregation.py` vs `statistical/descriptive.py` | **Three independent implementations** of "aggregate samples → seed → condition summaries" with different semantics (see §8). The pipeline uses the weakest one. |
| **D-11** | `statistical/resampling.py:120+`, `statistical/comparisons.py:134` | Permutation/sign-flip p-values use `mean(stat >= observed)` with **no +1 continuity correction**, so a p-value of exactly `0` is reported rather than `< 1/n_perm`. |
| **D-12** | `statistical/comparisons.py:134` | The paired sign-flip test draws from **global** `np.random` (`np.random.choice`), ignoring the `seed` passed to `TwoConditionComparisons`/`Resampler`. Paired results are **not reproducible** even with `analysis_seed` set. |
| **D-13** | `statistical/comparisons.py:111` | `effect_ci_lower` / `effect_ci_upper` are **always `None`** ("Can be computed via bootstrap if needed"). No effect-size CI is ever produced. `MultiConditionResult.omnibus_p_value` and `holm_adjusted_p_values` are declared but never populated — **no omnibus test exists**. |
| **D-14** | `statistical/multiplicity.py:81` | Holm monotonicity uses `results[-1]` — the last element of a **cross-family** accumulator. It happens to work only because families are processed consecutively; any refactor (e.g. sorting globally) silently breaks it. Also `compute_family_wise_error` (line 102) computes the *fraction of families with ≥1 rejection*, which is **not** FWER. |
| **D-15** | `preprocessing/quality_control.py:157,168` | `isinstance(seed, int)` is **False for `numpy.int64`**, which is what pandas yields from CSV. Every replicate would be flagged `INVALID_SEED`/`INVALID_SAMPLE` → `ERROR` → `CRITICAL_ERRORS` → `block_analysis`. QC would block all analysis if it were ever wired in. |
| **D-16** | `pipeline.py:88-104` | The extraction step loads `af3_condition_centric_extraction.py` by **walking parent directories with `importlib.util.spec_from_file_location`**. This hard-couples `af3_analysis` to a sibling repo's internal file path and is fragile to layout changes (it is already the only cross-repo coupling in the pipeline). |
| **D-17** | `io/structure_reader.py:237-272` | `_extract_atoms` picks **model 1 only** by default (`find_values("_atom_site.pdbx_PDB_model_num")[0]`). AF3 writes multiple ranked models per CIF; all but the top-ranked are silently dropped. The preceding block (lines ~230-243) is dead code (`find_loop` + `break` that does nothing). |
| **D-18** | `io/structure_reader.py:188-190` | `_extract_entities` maps entity→chain by taking the **first** matching chain (`break`). For homo-oligomers (one entity, several chains) **all but one chain are lost**, and `_extract_chains` (line 218) builds chains *from entities*, so chain-level metrics are wrong for such systems. |
| **D-19** | `structural/discovery.py:217-219` | Duplicate `(condition, seed, sample)` records are `continue`d (silently dropped, first kept), which makes the later duplicate detection (building `duplicate_predictions`) structurally **incapable of ever firing**. Duplicates are invisible in the QC report. |
| **D-20** | `structural/discovery.py:38` | `StructureRecord.checksum` defaults to `""` and is **never populated**, so structural predictions have no content-based provenance even though the field exists. |
| **D-21** | `preprocessing/canonicalize.py:156` | `_process_reference` hardcodes `sample-0` when keying the best-model/top-level reference to a prediction. If the reference does not correspond to sample 0, it is mis-keyed. |
| **D-22** | `preprocessing/canonicalize.py:52,90` | `self._measurements` is initialised and **never appended to**, so `CanonicalizationResult.measurements` is always an **empty DataFrame**. The canonicalizer produces no measurements. |
| **D-23** | `preprocessing/canonicalize.py:164-168` | `input_signature` hashes only the **source file paths**, not contents. It changes when a directory is moved and does not reflect composition — yet it is the field that backs "compositional comparability" claims. |
| **D-24** | `preprocessing/duplicate_resolution.py:134,147-162` | `apply_resolutions` compares records by **object identity** (`excluded is record`, line 162) and maintains an unused `excluded_ids` set keyed on `id(resolution)` (lines 147-150). It is only correct if the identical dict objects are passed back in — silent no-op otherwise. The checksum rule can never fire either, because callers supply `source_checksum`, not `checksum`, so resolution always degrades to `keep_first_record`. |
| **D-25** | `visualization/v2/factors.py:40` | `_LIGAND_TOKENS = {"dna"}` is a **hardcoded biological token**, contradicting the module docstring ("without hardcoding biological specifics"). It also gates `has_dna` on the metadata attribute being literally named `dna` (line 83). |
| **D-26** | `pipeline.py:219` | `groupby(["condition_id","condition_name","seed"]).mean()` aggregates using condition_name as part of the key; if two `condition_id`s ever share a name (or a name is missing), groups silently merge. |

### S3 — Robustness / maintainability

| ID | Location | Defect |
|---|---|---|
| D-27 | `visualization/orchestrator.py:43` | `except (FileNotFoundError, Exception)` — the second clause subsumes the first, and swallowing all exceptions turns a corrupt CSV into an empty DataFrame. Figures then silently skip. |
| D-28 | `visualization/orchestrator.py:143` | `global_metrics = ['pLDDT_mean','pLDDT_min','pae_mean','contact_prob_mean','ranking_score']` is a **hardcoded whitelist**; any newly extracted metric is silently excluded from all V1 figures. |
| D-29 | `preprocessing/aggregation.py:84`, `quality_control.py` (multiple), `canonicalize.py` | Widespread `.iterrows()` loops (O(n) Python-level iteration) over tables that reach hundreds of thousands of rows. |
| D-30 | `visualization/v2/factors.py:add_factor_columns` | Mutates the caller's DataFrame **in place** (`df[col] = ...`) and returns it, contradicting the "non-destructive" contract used by `io/wide_to_long.py`. Assigns by `.values` (positional) rather than by index. |
| D-31 | `statistical/eligibility.py:64` | `eligible_metrics` is initialised to `[]` and **never populated**, so the gate result always reports zero eligible metrics. `_has_compositional_comparability` (line ~98) is a stub that returns `True` whenever there are ≥2 conditions. |
| D-32 | `statistical/resampling.py` docstring | Documents `stat_func="mean_ratio"` but the implementation raises `ValueError` for anything except `"mean_diff"`. |
| D-33 | `structural/alignment.py:_compute_rmsd_no_alignment` | `rmsd_pre_alignment` is computed **after centroid centring**, not on raw coordinates. It is labelled as a pre-alignment RMSD but is not the raw quantity; it will be misread as "unbound RMSD". |
| D-34 | `structural/alignment.py:_match_residues` | Matching requires **identical `(chain_id, auth_seq_id)`**; there is no sequence alignment. `min_sequence_identity` therefore measures *fraction of shared residue IDs*, not sequence identity — a mutation or numbering shift makes structures "not comparable" rather than aligned. |
| D-35 | `structural/comparison.py:_compare_single` | Metric computation is wrapped in a bare `except Exception` → `"computation_error"` with **no logging**, hiding genuine bugs behind a status string. |
| D-36 | `io/raw_af3_reader.py` | `ArtifactRecord.af3_version` is declared but **never populated**; version provenance is absent. Regex classification uses greedy `(?P<cond>.+)`, which can over-capture on names containing `_seed-`. |
| D-37 | `af3_analysis/__init__.py:24,31,34` | `__all__` lists `"errors"` (no such module), `"statistics"` (directory is `statistical/`) and `"reporting"` (empty directory). Explains/reflects the stale installed copy (§2). `reporting/` is empty and referenced by V3's architecture doc. |
| D-38 | `pipeline.py:279-305` | `visualization_version` accepts only `v1|v2|both`. **V3 is not reachable from `run_pipeline`**, so the "3 versions" are not driven by one orchestrator. |

---

## 8. Divergent duplicate implementations (the reproducibility core)

The same scientific concept is implemented **three times** with different results:

| Concept | Implementation A | Implementation B | Implementation C |
|---|---|---|---|
| Sample→seed aggregation | `pipeline.py:219` — `groupby(...).mean()` over all numeric cols | `preprocessing/aggregation.py` — NaN-dropping mean, `n_samples`, `sample_sd`, `definedness` | — |
| Within-seed variance | `aggregation.py:219` — `mean(sample_sd) ** 2` | `statistical/descriptive.py:113` — `mean(sample_sd ** 2)` | — |
| Between-condition variance | `statistical/descriptive.py:112` — same as between-seed (D-07) | `pipeline.py` — `groupby(...).agg(['mean','std'])` | — |
| Condition summaries | `pipeline.py:224` — pandas `agg(['mean','std'])`, flattened | `aggregation.py:_build_condition_summaries` (mean/sd/median/min/max + seeds) | `statistical/descriptive.py:compute_condition_summaries` |
| Effect sizes | `statistical/comparisons.py` — Hedges' g / d_z (never called) | `structural/tables.py` — descriptive mean/median/sd only | V3 `analysis/effects.py` |
| Condition-name → factors | `visualization/v2/factors.py` regex (hardcoded `dna`) | `experiment_metadata.py` attribute-driven | `visualization_v3/adapter.py` stem→condition map |

`mean(sd)**2 ≠ mean(sd**2)` (Jensen), so A and B **disagree numerically** on
`within_seed_var`. The two variance decompositions in the codebase are not reconcilable as
written. This is the highest-value thing to consolidate.

---

## 9. Provenance and reproducibility

**Present:** `config_hash` in the run manifest; SHA-256 checksums on CIFs
(`structure_reader._compute_checksum`); `source_checksum` in `structural_predictions.csv`;
`analysis_seed` in `ResamplingConfig`.

**Missing or broken:**

- No `definedness` column in production outputs (§6.2).
- No exclusions recorded (§6.3).
- QC never runs, so no QC trail.
- `StructureRecord.checksum` never populated (D-20).
- `input_signature` hashes paths, not contents (D-23).
- Paired permutation tests ignore the seed (D-12).
- `ResamplingConfig` / `bootstrap_iterations` / `permutation_iterations` in `AnalysisConfig`
  are **never used** by `run_pipeline` — the config exposes knobs that do nothing.
- `af3_version` never captured (D-36).
- The **editable install hazard** (§2) means results depend on the caller's cwd.

---

## 10. Evidence appendix (real run `run_20260830_214934`)

```
conditions                     24
metrics_replicates.csv         2,448 rows
seed_aggregated.csv              240 rows (24 conditions x 10 seeds)
descriptive_stats.csv             24 rows, 76 columns  <- incl. chain_A_residues_mean/_std (D-01)
pairwise_comparisons.csv           0 rows (header only)                (D-02)
structural_predictions.csv     1,200 rows (24 conditions x 10 seeds x 5 samples)
structural_metrics.csv        18,400 rows (1,200 x ~15 metrics)
structural_comparisons.csv    69,000 rows, 230 unique comparison_id     (D-04)
                              5,750 comparisons vs 1,150 matched-sample (D-05)
structural_effects.csv           207 rows, descriptive only (mean/median/sd; no p, no CI, no effect size)
```

`descriptive_stats.csv` column excerpt showing the D-01 contamination:

```
condition_id, ranking_score_mean, ranking_score_std, pLDDT_mean_mean, pLDDT_mean_std, ...,
chain_A_residues_mean, chain_A_residues_std, chain_A_plddt_mean, chain_A_plddt_std, ...
```

---

## 11. Scientific decisions required (do not guess)

Per `AGENTS.md`, these must come from the thesis author, not from an implementer:

1. **Unit of analysis for inferential tests.** `Aggregator`'s docstring says one seed-level
   observation per condition × metric × scope. `structural/comparison.py` currently treats
   *predictions* (sample-level) as the unit. Which is canonical for the structural chapter?
2. **Sample→seed reduction.** Mean of samples within a seed? Median? Min-risk (AF3's own
   ranking)? This choice changes every downstream number.
3. **Paired vs unpaired condition comparison.** If seeds are paired across conditions
   (same seed index → same input noise), the paired sign-flip path is required — but its RNG
   is currently unseeded (D-12).
4. **Multiplicity family definition.** What constitutes a family — per metric? per scope?
   per figure? The code's `family_id` is never populated, so no family is defined.
5. **Which metric variants are distinct.** Whether `pLDDT_min/_max/_median` belong in the
   same correction family as `pLDDT_mean` (§5.3).
6. **Reference condition.** Currently defaulted to alphabetically-first
   (`structural/config.py` `first_condition`; V3 `resolve_reference`). Confirm this is the
   intended baseline.
7. **Interface metric definedness.** `DefinednessClassifier` hardcodes a list
   (`iptm, chain_pair_iptm, interface_pae, interface_pde, contact_prob`) as
   composition-undefined for monomers. Confirm this matches the thesis's definition.
8. **Chain identity convention.** `label_asym_id` is used for matching; confirm this is
   stable across all conditions (D-34).

---

## 12. Prioritized improvement opportunities

**Tier 1 — makes the reported numbers real**

1. Implement `pairwise_comparisons.csv` using the existing (but orphaned)
   `statistical/comparisons.py` + `multiplicity.py`, or explicitly document that the thesis
   reports no pairwise tests. *(D-02)*
2. Fix the sample×sample cross-product in `structural/comparison.py` to matched-sample
   pairing, and include `sample` in `comparison_id`. *(D-04, D-05)*
3. Exclude non-metric numeric columns from `pipeline._stage_run_analysis` (reuse
   `io/wide_to_long`'s classification instead of `select_dtypes`). *(D-01)*
4. Decide and enforce a **single** aggregation/descriptive implementation; delete or
   formally retire the other two. *(§8)*
5. Delete or fully implement `statistical/factorial_models.py`. *(D-03)*

**Tier 2 — makes the science defensible**

6. Carry `definedness` through to production tables; stop conflating
   `undefined_by_composition` with `missing_technical`. *(§6.2)*
7. Wire `QualityControl` in (after fixing the `np.int64` check) and record every exclusion.
   *(D-15, §6.3)*
8. Fix the unpaired bootstrap CI to resample the **difference**. *(D-06)*
9. Seed the paired permutation test from `analysis_seed`. *(D-12)*
10. Compute effect-size CIs; populate `n_permutations`-corrected p-values. *(D-11, D-13)*

**Tier 3 — structural correctness**

11. Fix entity→chain mapping for homo-oligomers; reconcile
    `NormalisedStructure.chains` with `_atom_site`. *(D-18)*
12. Decide and document the model-selection rule in `_extract_atoms` (currently top-ranked
    only). *(D-17)*
13. Populate `StructureRecord.checksum` and make `input_signature` content-based.
    *(D-20, D-23)*

**Tier 4 — environment & hygiene**

14. Remove/repoint the two editable installs so exactly one `af3_analysis` is importable;
    add a guard test. *(§2)*
15. Consolidate table naming (`chain_A_residues` vs `chain_A_residue_count_mean`). *(§5.2)*
16. Fix the stale `__all__` and resolve the empty `reporting/` package. *(D-37)*
17. Make V3 reachable from `run_pipeline` so all three versions share one orchestrator.
    *(D-38)*
18. Move the extraction step's cross-repo `importlib` file-path hack behind an explicit,
    documented interface. *(D-16)*

---

## 13. Suggested next artefact

A short **`SCIENCE_CONTRACT.md`** at the package root stating, in one page:

- the unit of analysis per chapter (seed vs prediction),
- the canonical sample→seed reduction,
- the definedness taxonomy and how each category is treated downstream,
- the multiplicity family definition,
- the reference condition and comparison pairing rule,
- the list of metrics that are scientifically distinct.

Most of the defects above become *unambiguous* once that document exists, and most of the
orphaned code in `preprocessing/`/`statistical/` becomes directly reusable.
