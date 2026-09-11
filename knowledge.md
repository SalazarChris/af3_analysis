# Thesis Project Knowledge

## Project

This repository contains the software supporting a biomedical engineering thesis involving AlphaFold 3.

The project has two main systems:

1. `af3inputscripts/` — AF3 input generation.
2. `af3_analysis/` — AF3 output analysis.

The current main menu/application entry point is inside `af3inputscripts/`.

`testdata/` contains development/test data.

## Current priority

The repository is being audited and cleaned before further thesis development.

The immediate goal is:

1. Understand the existing repository.
2. Document its architecture.
3. Identify obsolete/duplicated material.
4. Establish a clean development structure.
5. Preserve working functionality.
6. Avoid unnecessary refactoring.

## Scientific principles

Do not invent scientific methodology.

Do not silently change:

* unit of analysis
* inclusion/exclusion rules
* missing-data handling
* statistical methodology
* metric definitions
* experimental design
* biological interpretation

If a scientific decision is ambiguous, report the ambiguity instead of guessing.

## Development principle

Prefer the smallest change necessary.

Do not reorganize or delete files without explicit approval.

Do not modify source code during the initial repository audit.

## Important

The final repository structure has NOT yet been decided.

In particular, do not move the main menu from `af3inputscripts/` until its imports, paths, and execution behavior have been inspected.
