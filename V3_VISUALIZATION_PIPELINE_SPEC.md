# V3 VISUALIZATION PIPELINE — STEP-BY-STEP IMPLEMENTATION PROMPT

You are implementing **Visualization Pipeline V3** inside an existing AF3 analysis project.

Your task is to build a **new, additive visualization pipeline** that substantially improves the structural/geometric analysis and visualization capabilities of the project.

## CRITICAL RULE

**DO NOT MODIFY, REWRITE, REFACTOR, DELETE, REPLACE, OR CHANGE THE BEHAVIOR OF THE EXISTING ANALYSIS OR VISUALIZATION PIPELINE.**

V3 must coexist with the previous pipeline.

The existing pipeline is the baseline and must continue to work exactly as it currently does.

If you discover something in the old pipeline that could be improved, do NOT modify it as part of V3. Instead:

* work around it;
* create an adapter;
* create a compatibility layer;
* or document the issue for later.

The only exception is if V3 literally cannot function without a change. In that case, STOP before modifying the old code and report the dependency.

---

# PHASE 0 — UNDERSTAND THE MISSION

Before touching any code, understand the intended architecture.

The current project already performs AF3 analysis and generates an existing visualization suite.

V3 should add a new layer focused primarily on:

1. structural geometry;
2. structural differences between conditions;
3. local residue movement;
4. contact-map changes;
5. protein–DNA/interface geometry;
6. structural clustering;
7. seed reproducibility;
8. integration of structural geometry with AF3 confidence metrics.

The scientific question should be approximately:

> Do different experimental/modeling conditions produce reproducible differences in predicted protein geometry, interaction patterns, and AF3 confidence, beyond prediction variability?

Do NOT turn this into:

> Does this condition activate/inactivate the protein?

AF3 predictions alone do not justify that claim.

Use neutral terminology:

* structural difference;
* geometric displacement;
* altered contact pattern;
* interface change;
* structural variability;
* predicted structural cluster;
* confidence-associated structural difference.

---

# PHASE 1 — INSPECT THE EXISTING PROJECT

Do not write code yet.

Inspect the repository thoroughly.

Identify:

* main analysis package;
* visualization package;
* configuration files;
* CLI;
* master manifest;
* condition registry;
* condition metadata;
* extracted CSVs;
* structural outputs;
* confidence metrics;
* existing structural analysis;
* statistical modules;
* existing plotting utilities;
* test suite;
* output directory structure.

Specifically determine:

1. How conditions are represented.
2. How seeds are represented.
3. How predictions/samples/models are represented.
4. How AF3 output files are mapped to conditions.
5. How structures are currently located.
6. Whether structural metrics already exist.
7. Whether RMSD already exists.
8. Whether contact maps already exist.
9. Whether DNA/interface data already exist.
10. How existing figures are generated.
11. How existing configuration is loaded.
12. How existing output directories are created.
13. How existing statistical methods are implemented.

Do not infer architecture from filenames alone.

Read the actual code and tests.

---

# PHASE 2 — CREATE AN ARCHITECTURE MAP

Before implementation, produce an internal architecture map.

Document something equivalent to:

```
Existing AF3 outputs
        ↓
Existing analysis pipeline
        ↓
Existing CSV/tables/structures
        ↓
    V3 adapter
        ↓
V3 normalized data model
        ↓
Structural calculations
        ↓
V3 figure generators
        ↓
V3 output directory
```

The important boundary is:

```
EXISTING PIPELINE → V3 ADAPTER
```

The adapter is responsible for translating the existing project representation into the representation required by V3.

V3 should not force the old pipeline to change.

---

# PHASE 3 — IDENTIFY THE PROTECTED SURFACE

Create a list of existing modules/files that V3 must treat as protected.

At minimum protect:

* existing confidence analysis;
* existing statistical calculations;
* existing figure generators;
* existing configuration semantics;
* existing manifest logic;
* existing output structure.

Before continuing, determine how you will verify that these components remain unchanged.

Possible methods:

* git diff;
* regression tests;
* output checksums;
* existing test suite;
* before/after pipeline execution.

---

# PHASE 4 — CREATE THE V3 NAMESPACE

Create a dedicated V3 module/package.

Use an architecture compatible with the existing repository.

Conceptually:

```
visualization_v3/
    __init__.py
    config.py
    adapter.py
    validation.py
    structural/
        alignment.py
        rmsd.py
        displacement.py
        contacts.py
        interfaces.py
        clustering.py
    metrics/
    plots/
        ...
    reporting/
    tests/
```

Do NOT blindly use these exact filenames if the existing architecture suggests a better equivalent.

The important requirement is modular separation.

---

# PHASE 5 — CREATE V3 CONFIGURATION

V3 must have its own configuration.

Do not alter the old configuration format unless absolutely necessary.

The configuration should control things such as:

* V3 enabled/disabled;
* reference condition;
* enabled figures;
* alignment atom;
* contact threshold;
* minimum structural coverage;
* output DPI;
* output format;
* seed handling;
* clustering method;
* clustering parameters;
* region/site definitions;
* plot style.

Example conceptual configuration:

```
v3:
    enabled: true

    reference:
        condition: null

    figures:
        F01: true
        F02: true
        F03: true
        F04: true
        F05: true
        F06: true
        F07: true
        F08: true
        F09: true
        F10: true
        F11: true
        F12: true
        F13: true
        F14: true
        F15: true
        F16: true
        F17: true
        F18: true
        F19: true
        F20: true

    structure:
        alignment_atom: CA
        contact_distance: 8.0
        minimum_coverage: 0.80

    output:
        dpi: 300
        format: png
```

These are examples, not mandatory exact values.

Do not hard-code scientifically meaningful thresholds.

---

# PHASE 6 — BUILD THE V3 DATA ADAPTER

Implement an adapter that reads the existing pipeline outputs.

The adapter must construct a normalized representation containing, conceptually:

```
Dataset
Condition
Prediction
Seed
Structure
Entity
Chain
Residue
Atom
MetricObservation
```

The exact implementation is your choice.

The key requirement is that V3 understands:

```
condition
    → seed
        → prediction
            → structure
                → entity
                    → chain
                        → residue
                            → atom
```

Do not flatten this hierarchy unnecessarily.

---

# PHASE 7 — VALIDATE CONDITION AND SEED MAPPING

Before any structural calculations, validate metadata.

Ensure:

* every prediction belongs to a condition;
* every prediction has a valid seed when seed-level analysis is intended;
* condition-level rows are not interpreted as predictions;
* summary rows are not interpreted as predictions;
* duplicate rows are detected;
* invalid mappings are reported.

Do NOT use a simplistic assumption such as:

```
"If replicate_id contains seed-N, it must be valid."
```

The existing dataset has previously contained condition-level/top-level rows that do not contain seed identifiers.

Those must be explicitly classified.

Never silently discard them.

---

# PHASE 8 — VALIDATE STRUCTURAL INPUTS

For every structure determine:

* whether coordinates exist;
* which entities exist;
* which chains exist;
* sequence identity;
* residue mapping;
* atom mapping;
* missing residues;
* missing atoms;
* modified residues;
* DNA/RNA presence;
* ligand presence;
* ion presence.

Create a structural QC table.

For every failed comparison, record a reason.

Examples:

```
missing_structure
missing_chain
sequence_mismatch
insufficient_common_residues
missing_reference
unsupported_entity
invalid_coordinates
```

Do not silently drop failed comparisons.

---

# PHASE 9 — IMPLEMENT EXPLICIT REFERENCE RESOLUTION

The reference structure/condition must come from configuration.

Do not assume that:

```
WT = reference
```

or:

```
baseline = reference
```

unless the metadata explicitly says so.

Support:

```
condition → reference
```

and:

```
condition A → condition B
```

where appropriate.

For matched-seed comparisons, use:

```
reference(seed N) ↔ condition(seed N)
```

rather than arbitrary seed pairing.

---

# PHASE 10 — IMPLEMENT STRUCTURAL COMPARABILITY

Before calculating RMSD or displacement, determine the common structural space.

For every comparison calculate:

* common chains;
* common residues;
* common atoms;
* common sequence;
* coverage.

Example:

```
coverage =
    valid_common_residues /
    reference_residues
```

If coverage is below the configured threshold:

```
mark comparison LOW_COVERAGE
```

Do not silently produce a misleading RMSD.

---

# PHASE 11 — IMPLEMENT GLOBAL STRUCTURAL RMSD

Create a clean structural RMSD calculation module.

Support at least:

* Cα RMSD;
* configurable alignment atom.

Potential future support:

* backbone RMSD;
* heavy-atom RMSD.

Return both:

```
RMSD
```

and:

```
coverage
```

Do not mix incompatible entities.

For protein-vs-protein comparison, calculate protein structural RMSD.

For protein + DNA vs protein + DNA, allow separate protein geometry and interface analysis.

Do not automatically include DNA, ions, or ligands in protein RMSD.

---

# PHASE 12 — IMPLEMENT MATCHED-SEED STRUCTURAL ANALYSIS

For every condition/reference pair:

```
seed 1 → compare
seed 2 → compare
...
seed N → compare
```

Produce:

* RMSD per seed;
* coverage per seed;
* valid seed count;
* mean;
* median;
* SD;
* IQR;
* direction consistency where applicable.

The number of seeds must be discovered from the dataset.

Do not assume ten.

---

# PHASE 13 — IMPLEMENT PER-RESIDUE DISPLACEMENT

This is one of the highest-priority V3 analyses.

For each common residue:

```
displacement =
    distance(reference atom,
             condition atom)
```

Default atom:

```
Cα
```

but make it configurable.

Generate a table containing:

```
condition
residue_index
residue_name
displacement
n_valid
coverage
seed_consistency
```

Do not hard-code any biological residues.

---

# PHASE 14 — IMPLEMENT STRUCTURAL DISPLACEMENT HEATMAP

Generate:

```
condition × residue
```

with values representing:

```
mean or median displacement
```

The choice must be configurable.

Also consider a separate seed-consistency representation.

The objective is to distinguish:

```
large but inconsistent movement
```

from:

```
moderate but reproducible movement.
```

---

# PHASE 15 — IMPLEMENT CONTACT MAPS

Create residue-level contact maps.

For each structure:

```
contact(i,j) =
    distance(i,j) < configured threshold
```

Use a configurable distance threshold.

Calculate:

* reference contact map;
* condition contact map;
* gained contacts;
* lost contacts;
* unchanged contacts.

Create:

```
Δcontact_map
```

where:

```
+1 = gained
 0 = unchanged
-1 = lost
```

The exact encoding can be adjusted for visualization.

---

# PHASE 16 — IMPLEMENT CONTACT CHANGE SUMMARY

Calculate per condition:

* contacts gained;
* contacts lost;
* percentage changed;
* seed consistency;
* recurring residue-pair changes.

Generate a summary table.

Do not report a contact as a definitive biological mechanism.

Call it:

```
altered structural contact
```

or:

```
condition-associated contact change.
```

---

# PHASE 17 — IMPLEMENT ENTITY-AWARE INTERFACE ANALYSIS

Do not assume:

```
chain A = protein
chain B = DNA
```

Determine entity roles from metadata and/or structure parsing.

Support:

* protein–DNA;
* protein–RNA;
* protein–ligand;
* protein–protein;
* protein–ion.

For each interface calculate where possible:

* contacts;
* distances;
* contact frequency;
* residue identity;
* partner identity;
* seed consistency.

---

# PHASE 18 — IMPLEMENT PROTEIN–DNA INTERFACE ANALYSIS

When DNA exists, calculate:

* contacting protein residues;
* contacting DNA residues/bases;
* contact counts;
* minimum distances;
* contact persistence;
* condition-specific contact changes.

If interface confidence metrics are available, preserve them separately from geometry.

Do not assume every dataset contains DNA.

If DNA is absent:

```
skip the DNA figure
```

and record:

```
SKIPPED — no DNA entity available
```

Do not crash the entire V3 pipeline.

---

# PHASE 19 — IMPLEMENT LOCAL/PTM GEOMETRY

V3 must support configurable local regions.

A configuration might specify:

```
site:
    residue: 123
    modification: SUMO
    radius: 5
```

But this is metadata.

The code must NOT contain:

```
if residue == 123
```

or any OCT4-specific condition.

For each configured site/region calculate:

* local RMSD;
* local displacement;
* neighboring distances;
* contact changes;
* local confidence;
* interface proximity if applicable.

The implementation must work for arbitrary proteins and residues.

---

# PHASE 20 — IMPLEMENT DOMAIN / REGION MOTION

If region/domain annotations exist, calculate:

* region RMSD;
* centroid displacement;
* inter-region distance;
* relative orientation;
* optional angular changes.

Do not hard-code OCT4 domains.

Region definitions must come from metadata/configuration.

---

# PHASE 21 — IMPLEMENT PAIRWISE STRUCTURAL DISTANCES

Calculate a pairwise structural distance matrix.

Prefer:

```
pairwise RMSD
```

or another clearly documented structural metric.

Preserve metadata for every prediction:

* condition;
* seed;
* prediction ID.

Cache this calculation because it can be expensive.

---

# PHASE 22 — IMPLEMENT STRUCTURAL CLUSTERING

Use the pairwise structural distance matrix.

Potential approaches:

* hierarchical clustering;
* agglomerative clustering;
* k-medoids.

Do not arbitrarily claim that a certain number of clusters represents biological states.

If cluster count is configurable, document it.

Use terminology:

```
predicted structural cluster
```

not:

```
biological conformational state.
```

---

# PHASE 23 — IMPLEMENT MDS/PCA STRUCTURAL SPACE

Where sufficient comparable structures exist, generate a 2D structural-space representation.

Possible methods:

* MDS;
* PCA on a suitable coordinate representation.

Color/group using metadata:

* condition;
* factor;
* DNA presence;
* PTM state;
* seed.

Do not interpret axes as biological mechanisms.

The scientific question is simply:

> Do different conditions occupy distinguishable regions of predicted structural space?

---

# PHASE 24 — INTEGRATE CONFIDENCE WITH GEOMETRY

Build a dedicated analysis layer connecting:

```
structural geometry
        +
AF3 confidence
```

Possible combinations:

```
RMSD vs pLDDT
local displacement vs local pLDDT
structural effect vs confidence effect
contact change vs interface confidence
domain displacement vs PAE
```

Keep these metrics conceptually separate.

Do not merge them into one arbitrary "score."

---

# PHASE 25 — CREATE CONFIDENCE × GEOMETRY FIGURE

Create a figure showing categories such as:

1. low structural difference + high confidence;
2. high structural difference + high confidence;
3. low structural difference + low confidence;
4. high structural difference + low confidence.

The purpose is to distinguish:

```
confident structural differences
```

from:

```
uncertain structural differences.
```

Do not interpret high confidence as experimental validation.

---

# PHASE 26 — IMPLEMENT STRUCTURAL EFFECT SIZES

For every condition/reference comparison calculate, where supported:

* RMSD effect;
* local displacement effect;
* contact-change effect;
* interface-change effect;
* domain-motion effect.

Reuse the existing project's statistical framework.

Do NOT create a competing statistical framework unless necessary.

Include:

* effect estimate;
* uncertainty;
* valid seeds;
* direction consistency;
* structural coverage.

---

# PHASE 27 — IMPLEMENT FACTORIAL STRUCTURAL ANALYSIS

If the condition registry defines a factorial experiment, use those metadata fields.

Do not parse condition names.

For example, if factors are:

```
factor_A
factor_B
DNA
```

estimate structural effects corresponding to:

* main effects;
* two-way interactions;
* higher-order interactions where supported.

Reuse existing statistical methodology.

Do not assume the current OCT4 factors will exist in another project.

---

# PHASE 28 — IMPLEMENT SEED REPRODUCIBILITY

For each structural metric calculate:

* mean;
* median;
* SD;
* IQR;
* valid seeds;
* direction consistency.

Treat seeds as prediction-process robustness samples.

Do not describe them as biological replicates.

Do not describe them as a physical conformational ensemble.

---

# PHASE 29 — CREATE THE V3 FIGURE SUITE

Implement the following figure IDs.

## V3-F01

Structural data completeness / QC

## V3-F02

Global structural difference

## V3-F03

Matched-seed structural difference

## V3-F04

Per-residue structural displacement

## V3-F05

Structural displacement heatmap

## V3-F06

Contact-map difference

## V3-F07

Contact-change summary

## V3-F08

Protein–DNA/interface analysis

## V3-F09

Interface change map

## V3-F10

Local/PTM-site geometry

## V3-F11

Domain/region motion

## V3-F12

Structural clustering

## V3-F13

Structural similarity matrix

## V3-F14

Structural MDS/PCA

## V3-F15

Confidence × geometry

## V3-F16

Confidence change vs structural change

## V3-F17

Seed reproducibility

## V3-F18

Structural effect sizes

## V3-F19

Factorial structural effects

## V3-F20

Structure–confidence relationship matrix

Every figure must have:

* unique ID;
* clear title;
* axis labels;
* documented unit of analysis;
* documented metric;
* valid observation count;
* reference information;
* warnings where applicable.

---

# PHASE 30 — DO NOT DUPLICATE THE EXISTING FIGURES

The existing pipeline already contains figures covering:

1. QC completeness;
2. seed distributions;
3. factorial interaction analysis;
4. confidence relationships;
5. between-seed variability;
6. seed × condition heatmaps;
7. ECDF;
8. metric relationship matrix;
9. RMSD vs reference;
10. between-condition structural variability.

Do NOT delete or replace them.

If V3 creates a similar-looking figure, it must have a clearly different scientific purpose.

For example:

Old:

```
RMSD distribution
```

V3:

```
matched-seed structural effect + coverage + reproducibility
```

That is a meaningful extension.

---

# PHASE 31 — CREATE CENTRAL FIGURE MANIFEST

Create:

```
figure_manifest.json
```

or an equivalent machine-readable manifest.

For every figure store:

```
figure_id
name
status
output_path
input_data
unit_of_analysis
reference
metric
n_valid
warnings
skipped_reason
```

Possible statuses:

```
SUCCESS
SKIPPED
FAILED
```

Never hide missing figures.

---

# PHASE 32 — CREATE V3 OUTPUT DIRECTORY

Never write V3 figures into the old figure directory.

Use a separate directory such as:

```
run_YYYYMMDD_HHMMSS/
    v3/
        figures/
        tables/
        metadata/
        validation/
        logs/
        report/
```

Never overwrite old results.

---

# PHASE 33 — CREATE V3 TABLES

Generate at least:

```
structural_summary.csv
per_residue_displacement.csv
contact_changes.csv
interface_contacts.csv
structural_clusters.csv
confidence_geometry.csv
```

Add other tables where useful.

Every table must have a documented unit of analysis.

---

# PHASE 34 — ADD ROBUST LOGGING

The current project has previously experienced a Matplotlib rendering hang.

V3 must therefore log every major stage.

Example:

```
[V3] Starting F04
[V3] Loading structures
[V3] 250 structures discovered
[V3] Building residue mappings
[V3] Calculating displacement
[V3] Generating figure
[V3] Saving figure
[V3] Completed F04
```

If something becomes slow, the log must identify the exact stage.

---

# PHASE 35 — MAKE PLOTTING FAILURE-SAFE

A single failed figure must not terminate the entire V3 run.

Conceptually:

```
for figure:
    try:
        generate()
    except:
        record_failure()
        continue
```

However, do NOT silently swallow errors.

Record:

* exception;
* traceback;
* figure ID;
* input;
* affected condition if available.

---

# PHASE 36 — PROTECT AGAINST PATHOLOGICAL FIGURES

Before rendering:

* estimate number of categories;
* estimate number of labels;
* use sensible figure dimensions;
* avoid thousands of text annotations;
* avoid unnecessary legends;
* avoid giant tick-label objects;
* avoid redundant plotting calls.

Do not silently remove categories just to make a figure smaller.

If a visualization becomes too dense, create a more appropriate representation.

---

# PHASE 37 — ADD CACHING

Cache expensive operations:

* residue mappings;
* structure mappings;
* alignments;
* RMSD matrices;
* contact maps;
* interface maps;
* structural summaries.

Cache keys must depend on relevant inputs/configuration.

Never silently reuse stale cached results.

---

# PHASE 38 — CREATE UNIT TESTS

Test:

### Metadata

* condition mapping;
* seed mapping;
* prediction mapping;
* reference resolution.

### Structural mapping

* common residues;
* common atoms;
* sequence mismatch;
* missing coordinates;
* modified residues.

### Metrics

* RMSD;
* coverage;
* displacement;
* contact maps;
* interface contacts.

### Clustering

* distance matrix;
* clustering;
* metadata preservation.

### Robustness

* missing DNA;
* missing reference;
* monomer;
* multimer;
* incomplete structure;
* invalid seed;
* duplicate prediction;
* zero valid comparisons.

---

# PHASE 39 — CREATE V3 INTEGRATION TEST

Create a small fixture dataset.

Run the complete V3 pipeline.

Verify:

* pipeline finishes;
* figures are generated;
* tables are generated;
* figure manifest exists;
* logs exist;
* report exists;
* skipped figures are documented;
* failures do not terminate unrelated figures.

---

# PHASE 40 — RUN EXISTING PIPELINE REGRESSION

This step is mandatory.

Before V3:

```
run existing pipeline
```

Record:

* tests;
* output state;
* important files/checksums where practical.

Then implement V3.

Run the existing pipeline again.

Compare.

The result must explicitly state:

```
EXISTING PIPELINE PRESERVED: PASS
```

or:

```
EXISTING PIPELINE PRESERVED: FAIL
```

If FAIL:

STOP and diagnose before considering V3 complete.

---

# PHASE 41 — RUN V3 ON THE REAL DATASET

Only after tests pass.

Run V3 against the real dataset.

Do not modify the dataset.

Do not manually clean the dataset just for V3.

All exclusions must be programmatically recorded.

---

# PHASE 42 — AUDIT THE REAL OUTPUT

Inspect every generated figure.

Check for:

* unreadable labels;
* missing conditions;
* inconsistent ordering;
* incorrect legends;
* incorrect reference;
* incorrect seed mapping;
* structural coverage problems;
* pathological axes;
* misleading significance markers;
* missing metadata;
* inconsistent terminology.

Do not assume successful execution means scientific correctness.

---

# PHASE 43 — GENERATE V3 REPORT

Create:

```
V3_VISUALIZATION_REPORT.md
```

Include:

## Dataset

* conditions;
* seeds;
* predictions;
* structural coverage.

## QC

* successful mappings;
* failed mappings;
* low-coverage comparisons.

## Figures

For every figure:

* purpose;
* input;
* unit of analysis;
* valid observations;
* warnings.

## Structural findings

Only report descriptive/statistical observations supported by the data.

## Limitations

Explicitly state:

* confidence ≠ accuracy;
* structural difference ≠ biological function;
* AF3 seeds ≠ physical ensemble;
* static AF3 structures do not represent full dynamics;
* structural comparison depends on coordinate/sequence coverage.

---

# PHASE 44 — GENERICITY AUDIT

Before declaring completion, search the new V3 code for biological hard-coding.

There must be no logic specifically requiring:

```
OCT4
POU5F1
T101
S102
T235
S236
K123
chain A
chain B
DNA = chain B
WT
baseline
exactly 10 seeds
exactly 250 predictions
```

These are dataset values, not software logic.

If such values are required for the current dataset, they must come from:

* metadata;
* manifest;
* configuration;
* input structure;
* condition registry.

---

# PHASE 45 — SCIENTIFIC TERMINOLOGY AUDIT

Search all generated titles, labels, reports, and comments.

Remove unsupported statements such as:

```
activated
inhibited
active state
inactive state
functional state
```

unless an external validated reference explicitly supports them.

Prefer:

```
structural difference
geometric displacement
condition-associated change
predicted structural cluster
altered interface
confidence-associated structural difference
```

---

# PHASE 46 — FINAL PERFORMANCE AUDIT

Check:

* runtime;
* memory;
* repeated structure parsing;
* repeated alignment;
* repeated plotting;
* cache usage;
* figure rendering time.

The V3 pipeline should not unnecessarily recompute expensive calculations.

---

# PHASE 47 — FINAL DELIVERABLE

At the end, provide a concise implementation report containing:

## Files added

List every V3 file.

## Files modified

This should ideally be:

```
NONE
```

or only minimal additive integration files.

Explicitly explain every modification.

## Existing pipeline regression

```
PASS / FAIL
```

## V3 tests

```
PASS / FAIL
```

## V3 figures

List:

```
F01 ...
F02 ...
...
F20 ...
```

with:

```
SUCCESS / SKIPPED / FAILED
```

## Structural QC

Report:

* number of structures;
* valid comparisons;
* low-coverage comparisons;
* failed mappings.

## Known limitations

List remaining limitations honestly.

## Genericity

Explicitly confirm:

```
V3 contains no OCT4-specific biological logic.
```

---

# FINAL ARCHITECTURAL PRINCIPLE

The final architecture should conceptually be:

```
                     EXISTING PIPELINE
                           │
                           │
                           ▼
                 EXISTING OUTPUTS
                           │
                           │
                     ┌─────▼─────┐
                     │ V3 ADAPTER │
                     └─────┬─────┘
                           │
            ┌──────────────┼──────────────┐
            ▼              ▼              ▼
       CONFIDENCE      GEOMETRY       INTERFACES
            │              │              │
            └──────────────┼──────────────┘
                           ▼
                INTEGRATED V3 ANALYSIS
                           │
                           ▼
                 V3 FIGURE GENERATION
                           │
            ┌──────────────┼──────────────┐
            ▼              ▼              ▼
         FIGURES         TABLES        REPORT
                          
```

The central rule is:

> **V3 extends the project; it does not rewrite the project.**

Build incrementally.

Inspect first.

Modify second.

Test continuously.

Run the old pipeline after implementation.

Run V3 separately.

Audit the outputs.

Only then declare V3 complete.
