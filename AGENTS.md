# Thesis Repository — Agent Instructions

## Mission

This repository contains the code and data supporting a biomedical engineering thesis involving:

1. AlphaFold 3 input generation.
2. AlphaFold 3 output analysis.
3. Statistical analysis of AF3 confidence metrics.
4. Visualization and reporting of results.

The primary goal is to produce a scientifically defensible, reproducible thesis — not to maximize code complexity or automation.

---

## Critical Rule

**Do not invent or silently change scientific methodology.**

Agents may implement explicitly specified methodology, but must not independently change:

* unit of analysis
* inclusion/exclusion criteria
* missing-data handling
* statistical tests
* effect-size definitions
* multiple-comparison correction
* experimental design
* metric definitions
* biological interpretation
* conclusions

If a scientific decision is unclear, STOP and report the ambiguity rather than guessing.

---

## Repository Safety

Before modifying the repository:

1. Inspect the relevant files.
2. Understand imports and dependencies.
3. Identify generated files versus source files.
4. Make the smallest change that solves the task.
5. Do not perform large-scale refactoring unless explicitly requested.
6. Do not delete files during cleanup unless they have been explicitly classified as safe to delete.
7. Prefer archiving obsolete material over deleting it.

---

## Architecture

The project currently contains three main areas:

### `af3inputbuilder/`

AF3 input-building system.

This currently contains the main menu/application entry point as well as input-generation functionality.

### `af3_analysis/`

AF3 output-analysis system.

This contains extraction, metric processing, statistical analysis, and/or visualization functionality.

### `testdata/`

Test/example data used for development and validation.

Do not assume that testdata is disposable until its usage has been audited.

---

## Main Menu

The main menu currently lives inside `af3inputscripts/`.

Do not move it to the repository root until imports, path handling, package structure, and execution behavior have been audited.

A future architecture may use a root-level `app.py`, but this is a refactoring decision and must be validated before implementation.

---

## Scientific Constraints

The AF3 confidence-metric analysis must respect the established thesis specification.

Important principles include:

* Confidence metrics are not equivalent to structural accuracy.
* A change in confidence does not automatically imply a structural change.
* Missing values must not be imputed unless explicitly specified by the scientific specification.
* Structural missingness must be preserved.
* Mandatory row exclusions must be applied consistently.
* The correct unit of analysis must not be changed merely for computational convenience.
* Seed-level and sample-level quantities must not be treated as interchangeable.
* Tokenisation differences must be treated as potential confounding rather than ignored.
* Every reported numerical result must have traceable provenance to the underlying data.

---

## Agent Workflow

For every task:

1. Read this file.
2. Read relevant project documentation.
3. Inspect the existing implementation.
4. State the intended change briefly.
5. Implement the smallest appropriate change.
6. Run relevant tests or validation.
7. Inspect the resulting output.
8. Report:

   * files changed
   * what changed
   * tests/validation performed
   * remaining issues

Do not modify unrelated files.

---

## Code Quality

Prefer:

* existing project conventions
* small functions
* explicit interfaces
* deterministic behavior
* reproducible outputs
* clear error messages
* tests for important behavior
* pathlib for filesystem paths
* configuration over hard-coded paths

Avoid:

* duplicated logic
* hidden global state
* hard-coded user-specific paths
* unnecessary dependencies
* speculative refactoring
* rewriting working components without evidence

---

## Data Safety

Never overwrite raw/source data during analysis.

Generated outputs should be distinguishable from source data.

Never silently modify input datasets to make a pipeline run.

If required input data is missing, report the missing dependency clearly.

---

## Token/Context Efficiency

Do not read the entire repository for every task.

Start with:

1. this file
2. the relevant documentation
3. the files directly related to the task

Search before reading large files.

Do not repeatedly rediscover repository structure if it is documented.

---

## When Unsure

Ask for clarification or report the ambiguity.

A correct STOP is preferable to a scientifically incorrect implementation.

---

## Current Repository Status

The repository is currently undergoing structural cleanup and documentation.

Do not assume the final architecture has been established yet.

The initial cleanup must be read-only/audit-first.
