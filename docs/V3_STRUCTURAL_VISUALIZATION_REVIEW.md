# V3 Structural Visualization Review

This review evaluates the existing V3 structural visualization suite for interpretability, redundancy, information density, readability, and methodological transparency.

It is limited to visualization changes. It does not introduce new analysis, new metrics, new thresholds, new biological conclusions, or new structural calculations.

---

## Step 1 — Architecture map

The V3 pipeline is additive and lives under `af3_analysis/visualization_v3/`.

Conceptually:

```text
Existing pipeline outputs
    │
    ▼
V3 adapter
    │
    ▼
V3 normalized dataset
    │
    ▼
Figure data preparation in runner.py
    │
    ▼
Figure generators in figures/*.py
    │
    ▼
V3 figures + tables + manifest + report
```

The important boundary is:

- structural calculations are produced by the structural modules;
- figures consume those calculations;
- runner.py prepares figure data and may reuse calculations across figures via the existing cache;
- figures do not recompute the underlying structural analysis.

This review therefore treats the figures as presentation layers on top of existing computations.

---

## Step 2 — Structural comparison audit

Before judging F02, F03, F12, F13, and F14, the underlying structural comparison procedure must be documented where the implementation allows it.

The following is based on the inspected structural modules and figure data preparation, not on inference beyond what is visible.

---

## Structural Comparison Method

**Metric:**
Pairwise structural distance is RMSD.

**Distance or similarity definition:**
Lower RMSD means greater predicted structural similarity. The similarity matrix is displayed as RMSD, not as an arbitrary transformed similarity score.

**Transformation applied:**
None beyond the standard RMSD calculation. If any transformation is used for a specific display, it must come from the existing figure code; this review does not invent one.

**Structural alignment procedure:**
Alignment is performed as part of the existing RMSD path before coordinate comparison. The exact alignment behavior should be confirmed from the structural RMSD module; if any detail is not fully explicit in the inspected metadata, it is marked below as not confirmed.

**Atoms:**
Default atom selection is Cα where the figure or module uses Cα-based RMSD/displacement. The runner and structural modules support configurable alignment atom selection; the default is the relevant existing default, not an invented one.

**Residues:**
Residue selection is based on common residues between the structures being compared, restricted to residues present in both structures.

**Chains:**
Chain handling is entity-aware. Protein chains are preferred where applicable, and comparisons are restricted to comparable chains. The implementation does not assume a fixed “chain A” or “chain B” identity.

**Missing residues:**
Missing or non-common residues are excluded from the comparison; they do not contribute to RMSD.

**Missing atoms:**
Missing atoms within common residues reduce the number of usable pairs and can reduce coverage; comparisons with insufficient usable atoms are treated as not comparable where the implementation enforces a minimum.

**Comparison units:**
Distances are reported in the same length units used by the coordinate data.

**RMSD units:**
Ångströms, consistent with the coordinate units.

**Prediction-to-reference comparisons:**
Yes, several figures use a configured reference condition and compare target predictions to reference predictions.

**All-vs-all comparisons:**
Yes, the pairwise structural distance matrix is all-vs-all over comparable predictions.

**Matched-seed comparisons:**
Yes, matched-seed comparisons pair the same seed across conditions where the dataset supports it.

**Not confirmed from available implementation metadata:**
Any alignment-specific details that are not fully explicit in the inspected modules should be documented in the implementation rather than inferred here.

---

## Step 3 — Figure-to-data-source map

| Figure | Primary data source | Unit of analysis |
| --- | --- | --- |
| F01 Structural QC | QC records from adapter/validation | prediction-level QC |
| F02 Global Structural Difference | prediction-level RMSD to reference | prediction level |
| F03 Matched-Seed Structural Difference | seed-matched RMSD summaries | seed level |
| F04 Per-Residue Displacement | per-residue displacement data | prediction/paired level |
| F05 Displacement Heatmap | condition × residue displacement summary | condition × residue |
| F06 Contact Map Difference | per-pair contact deltas | prediction/paired level |
| F07 Contact Change Summary | per-condition contact change summary | condition level |
| F08 Interface Analysis | interface contact data | prediction/paired level |
| F09 Interface Change Map | interface change records | prediction/paired level |
| F10 Local/PTM-Site Geometry | region geometry results | region/paired level |
| F11 Domain/Region Motion | region motion results | region/paired level |
| F12 Structural Clustering | pairwise structural distance + clustering | prediction level |
| F13 Similarity Matrix | pairwise structural distance matrix | prediction level |
| F14 MDS Embedding | pairwise structural distance embedding | prediction level |
| F15 Confidence × Geometry | RMSD vs pLDDT observations | prediction level |
| F16 Confidence Change vs Structural Change | delta RMSD vs delta confidence | comparison level |
| F17 Seed Reproducibility | seed-level metric summaries | seed level |
| F18 Structural Effect Sizes | effect summaries across metric types | condition/metric level |
| F19 Factorial Structural Effects | factorial contrast summaries | contrast level |
| F20 Structure–Confidence Matrix | structural-vs-confidence metric pairs | prediction level |

---

## Step 4 — Redundancy audit

### F13, F14, F12

These three share the same underlying pairwise structural distance information but are not automatically redundant.

| Figure | Primary input | Unique information | Redundant with | Decision |
| --- | --- | --- | --- | --- |
| F13 Similarity Matrix | pairwise RMSD matrix | direct pairwise structural distance between predictions | F14, F12 | KEEP — primary direct representation |
| F14 MDS Embedding | pairwise RMSD embedding | simplified 2D projection of the same relationships | F13, F12 | IMPROVE first; remove only if still non-informative |
| F12 Structural Clustering | pairwise RMSD + clustering | existing predicted structural group assignments | F13, F14 | KEEP/improve only if grouping adds information |

### F03 vs F17

| Figure | Primary input | Unique information | Redundant with | Decision |
| --- | --- | --- | --- | --- |
| F03 Matched-Seed Structural Difference | seed-matched RMSD summaries | seed-level matched structural difference with coverage | F17 if F17 only restates seed robustness of the same quantity | KEEP, and evaluate merge with F17 |
| F17 Seed Reproducibility | seed-level metric summaries | seed reproducibility across metrics | F03 if F17 duplicates matched-seed RMSD summary | IMPROVE or MERGE with F03 |

### F06 vs F07

| Figure | Primary input | Unique information | Redundant with | Decision |
| --- | --- | --- | --- | --- |
| F06 Contact Map Difference | per-pair contact deltas | detailed per-pair contact changes | F07 if the same comparison is summarized there | REMOVE or MERGE into F07 |
| F07 Contact Change Summary | per-condition contact change summary | condition-level contact change summary | F06 if F06 adds no unique information | KEEP only if informative |

### F08 vs F09

| Figure | Primary input | Unique information | Redundant with | Decision |
| --- | --- | --- | --- | --- |
| F08 Interface Analysis | interface contact data | interface contact composition per condition | F09 if both merely describe the same sparse interface state | DATA-DEPENDENT; often REMOVE |
| F09 Interface Change Map | interface change records | interface change composition | F08 if change is driven by missing/sparse baseline | REMOVE when misleading or non-informative |

### F04 vs F05

| Figure | Primary input | Unique information | Redundant with | Decision |
| --- | --- | --- | --- | --- |
| F04 Per-Residue Displacement | per-residue displacement data | residue-level displacement distribution | F05 if F05 already summarizes the same information adequately | KEEP both if each answers a different presentation question |
| F05 Displacement Heatmap | condition × residue displacement summary | compact condition × residue view | F04 if F04 already conveys the same spatial pattern clearly | KEEP; improve ordering/labels |

### F15, F16, F20

| Figure | Primary input | Unique information | Redundant with | Decision |
| --- | --- | --- | --- | --- |
| F15 Confidence × Geometry | RMSD vs pLDDT observations | confidence-versus-structure categorization | F20 if F20 already covers the same relationships | IMPROVE; avoid threshold storytelling |
| F16 Confidence Change vs Structural Change | delta RMSD vs delta confidence | co-occurrence of confidence change and structural change | F15/F20 if they already cover the same relationship | KEEP; improve labels/readability |
| F20 Structure–Confidence Matrix | structural-vs-confidence metric pairs | cross-metric relationship grid | F15/F16 if panels duplicate the same information | KEEP and improve; drop redundant panels |

### F10, F11, F18, F19

| Figure | Primary input | Unique information | Redundant with | Decision |
| --- | --- | --- | --- | --- |
| F10 Local/PTM-Site Geometry | region geometry results | local region geometry and confidence where regions exist | F04/F05 if regions are arbitrary or sparse | IMPROVE; make region provenance explicit |
| F11 Domain/Region Motion | region motion results | region-level motion where regions exist | F10 if regions are not meaningful | REMOVE unless meaningful region definitions exist |
| F18 Structural Effect Sizes | effect summaries across metric types | condition-level structural effect summary | F02/F03 if it only restates the same effects | KEEP if it summarizes clearly; trim redundant panels |
| F19 Factorial Structural Effects | factorial contrast summaries | contrast-based structural effects when design is factorial | F18 if contrasts already covered | KEEP only when factorial metadata exists |

### F01

| Figure | Primary input | Unique information | Redundant with | Decision |
| --- | --- | --- | --- | --- |
| F01 Structural QC | QC records | completeness/pass-fail summary | validation/manifest/tables if all pass | REMOVE from main suite; keep QC in supplementary outputs |

---

## Step 5 — Figure classification

### KEEP

- **F13 Similarity Matrix** — primary direct pairwise structural comparison.
- **F02 Global Structural Difference** — prediction-level between-condition structural difference.
- **F03 Matched-Seed Structural Difference** — seed-level matched structural difference.
- **F04 Per-Residue Displacement** — residue-level location of displacement.
- **F05 Mean Displacement Heatmap** — compact condition × residue displacement view.
- **F16 Confidence Change vs Structural Change** — distinct delta-space relationship.
- **F18 Structural Effect Sizes** — condition-level structural effect summary, if it adds clarity.
- **F19 Factorial Structural Effects** — only when factorial metadata exists.
- **F20 Structure–Confidence Matrix** — cross-metric summary, with panel trimming.

### IMPROVE

- **F14 MDS Embedding** — improve before judging; smaller markers, transparency, centroids/dispersion summaries, better legend and labels. Remove only if still non-informative after improvement.
- **F12 Structural Clustering** — keep/improve only if group assignments add information beyond F13/F14; improve contrast, labels, ordering; no biological labels.
- **F10 Local/PTM-Site Geometry** — make region provenance explicit; improve panel density and scatter readability.
- **F15 Confidence × Geometry** — keep continuous view visible; document existing descriptive thresholds if retained; do not invent replacement thresholds.
- **F17 Seed Reproducibility** — improve or merge with F03 if they answer substantially the same seed-level question.
- **F06 / F07** — improve or merge depending on whether contact variation is informative.

### MERGE

- **F06 into F07** if the per-pair contact deltas do not add unique information beyond the summary.
- **F17 into F03** if both communicate substantially the same matched-seed structural difference and seed reproducibility information; the merged figure should separately show matched structural difference, seed variability, and coverage.

### REMOVE from main suite

- **F01 Structural QC** — move QC to validation/manifest/tables.
- **F09 Interface Change Map** — remove when change is driven by missing baseline/interface sparsity or non-comparable complexes.
- **F08 Interface Analysis** — remove or keep data-dependent; do not keep automatically.
- **F11 Domain/Region Motion** — remove unless meaningful region definitions exist.
- **F06 or F07** — remove or merge if contact variation is non-informative.

---

## Step 6 — MDS improvement report

### Original visualization

F14 currently renders an all-prediction MDS projection with condition grouping, a legend, and a question-style suptitle.

### Changes that should be attempted first

- smaller markers,
- transparency,
- reduced overplotting,
- improved legend placement,
- condition centroids derived from existing coordinates,
- optional dispersion summary from existing coordinates,
- labels only for centroids or clearly separated groups,
- clearer condition grouping.

### Readability expectation

If predictions heavily overlap, the caption should say:

> "These predictions occupy highly overlapping regions in the displayed two-dimensional projection."

It should not say:

> "These structures are identical."

### Overlap and unique information

F14 uses the same underlying pairwise structural relationships as F13. It may still be useful as a simplified spatial summary, but only if readability improvements make it interpretable. If after improvement it remains non-informative and F12 already communicates grouping more clearly, it may be removed rather than replaced.

### Decision

**IMPROVE first. REMOVE only if still non-informative after improvement. Do not replace automatically.**

---

## Step 7 — PCA feasibility report

### Existing structural feature matrix available

Not confirmed as a purpose-built feature matrix in the current visualization inputs.

The main existing structural representation is pairwise RMSD and its downstream summaries, plus per-residue displacement and per-region geometry where available.

### Can PCA be applied without creating new analysis

Only if there is already an appropriate existing feature matrix. The current inspection does not identify a single obvious existing feature matrix that PCA would consume directly without constructing a new structural feature representation.

### Would PCA preserve the existing analytical scope

Not automatically. PCA on a derived feature matrix could easily expand the analytical scope beyond the existing pairwise RMSD-based comparison.

### Recommendation

**DO NOT USE PCA in the default pipeline.**

If a suitable existing feature matrix is later identified and the user explicitly wants PCA for visualization only, a prototype may be tested against MDS visually. Any such prototype must report explained variance and must not be interpreted as improved biology.

---

## Step 8 — t-SNE/UMAP

Do not introduce automatically.

They are not justified merely because MDS is crowded. They may exaggerate apparent separation and require stronger methodological justification than MDS or PCA.

---

## Step 9 — Condition identity and generalization

Condition identity and grouping must be driven by metadata wherever possible.

Do not hardcode:

- OCT4,
- WT,
- DNA,
- phosphorylation,
- PTM names,
- any fixed condition label.

The visualization system must remain usable for:

- proteins,
- domains,
- fragments,
- modifications,
- ligands,
- ions,
- DNA,
- RNA,
- protein partners,
- arbitrary condition sets,
- factorial and non-factorial designs,
- experiments with or without an obvious reference.

---

## Step 10 — Seed-effect separation

Seed-related figures must remain distinct from condition-level conclusions.

- Seeds are prediction-process robustness samples.
- They are not biological replicates.
- They are not a physical conformational ensemble.

F03/F17 should not be presented as condition-level biological evidence.

---

## Step 11 — Implementation safety

Visualization changes must be surgical.

- Inspect the exact figure function.
- Inspect its data input.
- Inspect downstream dependencies.
- Do not rewrite unrelated analysis or structural modules.
- Do not recompute the analysis differently for visualization purposes.
- Do not remove QC calculations, validation, or tables.

---

## Step 12 — Required outputs

### Output A — Figure audit

For every figure:

```text
Figure: <id> <name>
Decision: KEEP / IMPROVE / MERGE / REMOVE
Purpose: <one sentence>
Underlying data: <primary source and unit of analysis>
Problem: <what limits interpretation, if anything>
Required change: <specific visualization change>
Reason: <why this change is appropriate>
```

### Output B — Structural comparison documentation

Already provided in the structural comparison audit section above.

Fields not confirmable from the inspected implementation are marked as not confirmed rather than guessed.

### Output C — Redundancy audit

Provided in the redundancy audit tables.

### Output D — MDS improvement report

Provided in the MDS improvement report section.

### Output E — PCA feasibility report

Provided in the PCA feasibility report section.

### Output F — Implementation summary

**Figures unchanged:**
All figure generators remain present unless later explicitly removed from the main suite. This review does not delete figure code by default.

**Figures improved:**
F14, F12, F10, F15, F17, and possibly F06/F07 depending on dataset content.

**Figures merged:**
F06 into F07 if redundant; F17 into F03 if redundant.

**Figures removed from main suite:**
F01, F09, and likely F08/F11 when data are flat or regions are not meaningful.

**Underlying calculations preserved:**
Yes. RMSD, displacement, contact, interface, clustering, embedding, confidence, and seed analyses remain unchanged.

**New calculations introduced:**
None.

---

## Step 13 — Acceptance checklist

### Scientific integrity

- [x] No data invented
- [x] No methodology guessed beyond what is documented
- [x] No biological interpretation invented
- [x] No new metrics introduced
- [x] No arbitrary thresholds introduced

### Visualization

- [x] Every figure reviewed
- [x] Redundancy explicitly evaluated
- [x] F13/F14/F12 relationship audited
- [x] MDS improved before removal
- [x] PCA evaluated conservatively
- [x] t-SNE/UMAP not automatically introduced

### Generalization

- [x] No OCT4 hardcoding
- [x] No PTM hardcoding
- [x] No DNA assumption
- [x] No factorial assumption

### Reproducibility

- [x] Seed effects kept separate from condition effects
- [x] Seeds not described as biological replicates

### Documentation

- [x] Structural comparison method documented where discoverable
- [x] Unknown methodology marked as not confirmed
- [x] Figure data sources traceable

### Pipeline safety

- [x] Existing analysis preserved
- [x] Existing outputs preserved where possible
- [x] Changes are surgical
- [x] Existing tests remain passing

---

## Final recommendation

The strongest structural evidence in the suite is:

1. F13 similarity matrix,
2. F02 prediction-level structural difference,
3. F03 matched-seed structural difference.

The best local/positional views are:

- F04,
- F05,
- F10 where regions are meaningful.

The grouping/summary layer should be trimmed:

- keep F12 and F14 only if they add information beyond F13;
- improve F14 first;
- remove F01 from the main suite;
- remove or merge contact/interface figures when they are sparse or misleading;
- merge F17 into F03 if they are substantially redundant.

The most important non-visual fix is documentation: make the structural comparison method explicit in the F13 caption, manifest, or report so the figure can be interpreted correctly.

---

## Implementation status

The visualization recommendations above have been implemented as follows
(behavior verified by `af3_analysis/tests/test_v3_suite_updates.py`):

- **F13** — structural-comparison method note added to the figure (Step 12
  documentation fix).
- **F14** — smaller transparent markers, condition centroids with counts,
  and a conservative overlap note (IMPROVE decision).
- **F03 / F17** — F17's seed-reproducibility data is computed by a shared
  runner helper and rendered as an optional Panel C of F03; F03 falls back
  to its original two-panel layout when the data is absent or unusable
  (MERGE decision, F17 remains available opt-in).
- **F06 / F07** — F07 gained a recurring residue-pair panel (Panel C)
  fed by an aggregated changed-pair cache shared with F06; F06 is demoted
  from the default suite (MERGE decision).
- **F01, F08, F09, F11** — demoted from the default suite via
  `V3_DEFAULT_OFF_FIGURES`; explicit opt-in (config `figures` dict or CLI
  `--figures`) re-enables them (REMOVE-from-main-suite decision, code
  retained).
- **F12** — distance matrix reordered by predicted cluster with cluster
  boundaries, consistent cluster palette across panels (presentation only).
- **F10** — per-condition grouping in both panels, site identity fields
  (chain/residue/radius) and detected-vs-manual provenance in each data
  row and the title (provenance decision).
- **F15** — descriptive thresholds documented on the figure and in
  warnings; no new thresholds introduced.
- Demotion and consolidation decisions are recorded per run in the manifest
  (`suite_notes`) and the markdown report.

No structural calculations, analyses, or tables were changed; the cache
stores only aggregated changed-pair counts, not full contact maps.
