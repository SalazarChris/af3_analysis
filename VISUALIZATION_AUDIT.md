# Visualization Audit

**Date:** 2026-08-26
**Status:** Read-only audit — no files modified.

---

## Data Context

All figures consume data produced by the analysis pipeline:

| File | Shape | Key columns |
|------|-------|-------------|
| `seed_aggregated.csv` | 81 rows × 22 cols | 9 conditions × ~10 seeds; 5 global metrics + 3 chain pLDDT + chain residues |
| `descriptive_stats.csv` | 9 rows × 39 cols | Condition-level mean/SD of seed means; no `condition_name` column |
| `pairwise_comparisons.csv` | **0 rows** × 4 cols | Always empty (placeholder) |

**Schema metrics selected by orchestrator:**
Global: `pLDDT_mean`, `pLDDT_min`, `pae_mean`, `contact_prob_mean`, `ranking_score`
Chain: `chain_A_plddt`, `chain_B_plddt`, `chain_C_plddt`

**Conditions:** 8 POU experimental conditions (baseline, ±DNA, ±pSEP102, ±pTPO101)

---

## F1 — QC Completeness Heatmap

**File:** `fig_qc_completeness.png`

### 1. Scientific purpose
Show which metrics are available (non-null) for each condition, identifying missing-data patterns before downstream analysis.

### 2. Exact input data
`seed_aggregated.csv` — one row per (condition × seed). For each metric, computes `cond_df[m].notna().mean()` = fraction of seeds with non-null values.

### 3. Unit of observation
One row per condition (aggregated over seeds). Heatmap cell = fraction of seed-level rows where the metric is non-null.

### 4. Visual elements
- **Heatmap cells:** Fraction non-null (0–1) via viridis colormap
- **Annotations:** Numeric fraction inside each cell
- **Right margin:** `n=XX` annotation showing total seed count per condition
- **X-axis:** Metric names (rotated 30°)
- **Y-axis:** Condition names

### 5. Current strengths
- Correctly identifies that chain_B and chain_C pLDDT are NaN for monomeric conditions (pou_baseline, pou_dna, etc.)
- Clean layout with sample counts annotated outside the heatmap
- Appropriate use of `notna().mean()` for missingness fraction

### 6. Current visual problems
- **Viridis colormap is misleading.** Viridis implies a continuous gradient of "quality," but missingness is essentially binary for this dataset (either 100% present or ~0%). The color scale suggests nuanced differences that don't exist.
- **Most cells are identical (1.0).** With 10 seeds per condition and no extraction failures, all global metrics show 1.0. The plot has low information density — it confirms everything is present but doesn't reveal anything unexpected.
- **Chain residue columns are excluded from schema but not from the heatmap.** The `_get_primary_metrics` heuristic may or may not include them depending on the fallback path.

### 7. Potentially misleading elements
- A reader unfamiliar with the data might interpret the viridis gradient as representing data quality rather than mere presence/absence.
- The `n=XX` annotation could be confused with the sample count rather than the number of seed-level observations.

### 8. Hard-coded assumptions
- Metrics are selected via `_get_primary_metrics()` heuristic (contains "mean", "plddt", "score", etc.)
- Conditions are ordered by `get_condition_order()` which uses a hard-coded `ptm_rank` dict for POU conditions
- Figure height scales with `len(conditions) * 0.6` inches

### 9. Currently valid?
**Partially.** The plot renders correctly and accurately reports missingness. However, for this specific dataset it has low information value because nearly all cells are 1.0. It becomes more useful when extraction failures or condition-dependent missingness exist.

### 10. Recommendation: **KEEP with minor styling adjustment**
The concept is sound. The viridis colormap should be replaced with a sequential blue palette (e.g., `Blues_r`) or a binary palette to avoid implying false precision. The figure is most valuable as a sanity check before analysis.

---

## F2 — Seed Distributions (Box + Strip)

**File:** `fig_seed_distributions.png`

### 1. Scientific purpose
Show the spread and central tendency of seed-level metric values across conditions, allowing visual comparison of distributions.

### 2. Exact input data
`seed_aggregated.csv` — one row per (condition × seed). Each point is a seed-level mean.

### 3. Unit of observation
Each point = one seed's mean value for a given condition and metric. 10 points per condition (10 seeds).

### 4. Visual elements
- **Box plots:** Quartiles of seed-level values per condition (white fill, `#f5f5f5`)
- **Strip plots:** Individual seed-level points, jittered horizontally, colored by condition
- **Panels:** One panel per metric (A through E for 5 global metrics, plus chain metrics)
- **Y-axis:** Metric value
- **X-axis:** Condition labels (display names, rotated 30°)

### 5. Current strengths
- Box + strip combination is appropriate for small n (10 seeds)
- `showfliers=False` prevents outlier dots from dominating
- Conditions are consistently ordered and colored across all panels
- Panel labeling (A, B, C...) is clear

### 6. Current visual problems
- **Box plots are misleading with n=10.** The box shows quartiles, but with only 10 points the quartile boundaries are unstable and the "whiskers" (1.5×IQR) may not be meaningful.
- **Jitter overlap.** With 10 points per condition and `jitter=0.25`, some points overlap, especially when values are close. The `alpha=0.6` helps but doesn't fully resolve it.
- **Chain metrics with mostly NaN values.** chain_B and chain_C pLDDT panels will show empty or near-empty plots for most conditions (only 2-3 conditions have data), creating wasted space.
- **No indication of sample size per condition.** The reader cannot tell whether a condition has 10 seeds or fewer.
- **X-axis label is empty** (`ax.set_xlabel("")`) — the x-axis represents conditions, which should be labeled.

### 7. Potentially misleading elements
- The box plot median line with n=10 can be mistaken for a robust central tendency estimate when it's actually just the 5th and 6th sorted values.
- The `showfliers=False` hides legitimate extreme values that might be scientifically important.

### 8. Hard-coded assumptions
- `figsize=(11, 4.5 * len(metrics))` — height per panel is fixed
- `alpha=0.6`, `size=5`, `jitter=0.25` — overlap parameters are fixed
- Box plot color `#f5f5f5` (near-white) is hard-coded
- Condition ordering via `get_condition_order()` uses POU-specific rank dict

### 9. Currently valid?
**Yes.** The plot renders correctly for the available data. The chain_B and chain_C panels will be sparse but not incorrect.

### 10. Recommendation: **MODIFY**
- Add `n=XX` annotations above each box to show seed count
- Consider replacing box plots with violin plots or strip-only plots for n=10
- Add x-axis label ("Condition")
- Remove or collapse chain metric panels that are mostly empty
- Increase jitter or use swarm plot for better point separation

---

## F3 — Factorial Interaction Plot

**File:** `fig_factorial_interaction.png`

### 1. Scientific purpose
Show how the effect of DNA substrate (factor 1) varies across PTM modifications (factor 2), revealing interaction effects in the factorial design.

### 2. Exact input data
`descriptive_stats.csv` — one row per condition. Uses `{metric}_mean` columns (condition-level mean of seed means).

### 3. Unit of observation
Each point = one condition's mean value. Lines connect conditions within the same PTM level across DNA states.

### 4. Visual elements
- **Lines:** Connect "No DNA" (0) to "DNA Present" (1) for each PTM modification
- **Points (markers):** Condition-level means at each factor level
- **Hue:** PTM modification type (Baseline, pSEP102, pTPO101, pTPO101+pSEP102)
- **Panels:** One per metric
- **X-axis:** DNA substrate (0/1)
- **Y-axis:** Mean metric value

### 5. Current strengths
- The interaction concept is scientifically appropriate for this factorial design
- `ptm_pretty` labels provide clear scientific names for PTM states
- The line-per-PTM layout makes interactions visually intuitive (non-parallel lines = interaction)
- Factor inference from condition names is clever and works for the POU naming convention

### 6. Current visual problems
- **`condition_name` is missing from `descriptive_stats.csv`.** The orchestrator patches this by joining from `seed_aggregated.csv`, but the factor inference code in the plot function itself also reconstructs `factor_PTM` and `factor_DNA` from condition names. This double reconstruction is fragile.
- **Points are not visible.** `sns.lineplot` with `marker='o'` draws markers, but with only 2 x-values (0 and 1) per line, the markers may overlap or be hidden by the line.
- **No error bars or confidence intervals.** The plot shows point estimates only, with no indication of uncertainty.
- **Factor levels are categorical but plotted on a continuous 0–1 axis.** The x-axis shows numeric 0/1, which could imply a continuous relationship between "No DNA" and "DNA Present."

### 7. Potentially misleading elements
- The connecting lines imply a continuous transition between DNA=0 and DNA=1, which is categorical. A reader might interpret the slope as a rate of change.
- The absence of error bars means the reader cannot assess whether the DNA effect is meaningful relative to seed-level variability.
- The `sns.lineplot` aggregation (default `estimator=mean`) is redundant since the input is already a single mean per condition.

### 8. Hard-coded assumptions
- `factor_DNA` is derived from whether condition name ends with `_dna`
- `factor_PTM` is derived by stripping `pou_` prefix and `_dna` suffix
- `ptm_pretty` dict maps 5 specific PTM states to display names
- `main_factor = factor_cols[0]` always uses the first factor column as x-axis
- X-axis is always `factor_DNA` (binary 0/1)
- Legend title is always "PTM Modification"

### 9. Currently valid?
**Conditionally.** The plot renders correctly for the POU dataset. However, the factor inference is entirely POU-specific and would break for any other experimental design. The plot is scientifically meaningful for this dataset but would need generalization for other experiments.

### 10. Recommendation: **MODIFY**
- Add error bars (SE or CI) to each point
- Use categorical x-axis (strings, not 0/1 numbers)
- Make the factor inference configurable rather than hard-coded to POU naming
- Consider adding individual seed-level points as jittered background

---

## F4 — Effect Size Forest Plot

**File:** `fig_effect_size_forest.png`

### 1. Scientific purpose
Show pairwise effect sizes (Hedges' g or mean difference) between conditions with confidence intervals, enabling visual assessment of which comparisons are statistically meaningful.

### 2. Exact input data
`pairwise_comparisons.csv` — currently **always empty** (0 rows).

### 3. Unit of observation
Each row = one pairwise comparison (condition vs reference) for one metric.

### 4. Visual elements
- **Points:** Effect size estimate
- **Error bars:** Confidence intervals (if `hedges_g_ci_lower/upper` columns exist)
- **Vertical dashed line:** Zero effect (null hypothesis)
- **Y-axis:** Contrast labels ("Condition A vs Reference")
- **Panels:** One per metric

### 5. Current strengths
- The forest plot format is standard and appropriate for effect size visualization
- The zero-reference line aids interpretation
- Sorting by effect size aids comparison

### 6. Current visual problems
- **The figure is never generated.** `pairwise_comparisons.csv` is always empty because the pipeline writes a dummy empty DataFrame. The function returns immediately on `len(pairwise_comparisons_df) == 0`.
- **Dead code.** This entire function is effectively dead code until real pairwise comparisons are implemented.

### 7. Potentially misleading elements
- None currently (figure is never produced). If it were produced with the dummy empty data, it would show nothing.

### 8. Hard-coded assumptions
- Falls back to `diff_mean` if `hedges_g` column doesn't exist
- Uses `get_display_label()` for condition names
- Contrast label format: `"Condition vs Reference"`

### 9. Currently valid?
**No.** The figure is never produced because the input data is always empty.

### 10. Recommendation: **KEEP but document as pending**
The plot design is sound. It needs real pairwise comparison data to function. This is a pipeline issue (P4 from DATA_PIPELINE_AUDIT.md), not a visualization issue.

---

## F5 — Variability (Between-Seed SD)

**File:** `fig_variability.png`

### 1. Scientific purpose
Show how much each condition varies across seeds for each metric, identifying conditions with high or low prediction consistency.

### 2. Exact input data
`descriptive_stats.csv` — uses `{metric}_std` columns (SD of seed-level means within each condition).

### 3. Unit of observation
Each bar = one condition's between-seed SD for one metric.

### 4. Visual elements
- **Bars:** Height = SD of seed-level means
- **Color:** Condition-colored bars (consistent with other figures)
- **Panels:** One per metric
- **Y-axis:** SD value
- **X-axis:** Condition names

### 5. Current strengths
- The concept is scientifically valuable — between-seed variability is a key indicator of prediction reliability
- Consistent condition ordering and coloring with other figures
- Clear panel labeling

### 6. Current visual problems
- **No error bars on the SD bars themselves.** The SD is a point estimate; with only 10 seeds, the SD has substantial uncertainty that is not communicated.
- **Y-axis always starts at 0.** For metrics where all conditions have similar low variability, the bars are all tiny and indistinguishable.
- **Condition labels are repeated for every panel.** With 5+ panels, the x-axis labels consume significant vertical space.
- **`condition_name` is missing from `descriptive_stats.csv`.** The orchestrator patches this, but the plot function itself requires it.

### 7. Potentially misleading elements
- Low SD bars might be interpreted as "precise" when they could simply reflect that the metric has low dynamic range (e.g., if all conditions produce nearly identical pLDDT values, the SD will be low regardless of prediction quality).
- The bar chart format implies that SD is a "bigger is worse" quantity, which is not always the case.

### 8. Hard-coded assumptions
- Uses `{metric}_std` column naming convention from `descriptive_stats.csv`
- Condition ordering via `get_condition_order()` (POU-specific)
- Bar colors from `get_condition_style_map()` (husl palette)

### 9. Currently valid?
**Yes.** The plot renders correctly for the available data, though the interpretation caveats above apply.

### 10. Recommendation: **MODIFY**
- Add confidence intervals on the SD estimates (e.g., via bootstrap)
- Consider log scale for y-axis if SD values span orders of magnitude
- Add `n=XX` annotations to show seed count per condition
- Consider a heatmap alternative for many metrics (rows=conditions, columns=metrics, color=SD)

---

## F6 — Seed Trajectories

**File:** `fig_seed_trajectories.png`

### 1. Scientific purpose
Show how individual seeds behave across conditions, revealing seed-level consistency or instability in the ranking of conditions.

### 2. Exact input data
`seed_aggregated.csv` — one row per (condition × seed). Each line connects one seed's values across conditions.

### 3. Unit of observation
Each line = one seed's metric values across all conditions.

### 4. Visual elements
- **Lines:** Connect one seed's values across conditions (hue = seed ID)
- **Color:** Discrete palette (tab10 or tab20) for seed identification
- **Alpha:** Transparency scaled by seed count (`5.0 / n_seeds`)
- **Panels:** One per metric
- **Legend:** Seed IDs (only on first panel)

### 5. Current strengths
- The trajectory concept is valuable — it shows whether the same seed consistently produces high or low values across conditions
- Alpha transparency helps with overplotting
- Legend is placed outside the plot area to avoid data overlap

### 6. Current visual problems
- **Severe overplotting.** With 10 seeds × 8 conditions, there are 10 lines crossing 8 x-positions. The lines frequently overlap, especially in the middle of the distribution.
- **Lines between conditions are meaningless.** Connecting seed-1's value in "Baseline" to seed-1's value in "POU + DNA" implies a trajectory, but the conditions are categorical — there is no meaningful ordering that makes a "trajectory" interpretable.
- **Legend is cluttered.** 10 seed entries in the legend on the first panel are hard to distinguish.
- **Color assignment is arbitrary.** Seeds are assigned colors by tab10/tab20, but the seed numbers don't have inherent ordering significance.
- **Same information as F2 but less readable.** The box+strip plot (F2) shows the same distribution information more clearly.

### 7. Potentially misleading elements
- The connecting lines strongly imply continuity or ordering between conditions, which is categorical. A reader might interpret the slope between conditions as a meaningful rate of change.
- The trajectory format might suggest a time-series or sequential relationship between conditions.

### 8. Hard-coded assumptions
- `alpha_val = max(0.2, min(0.8, 5.0 / n_seeds))` — transparency formula
- `linewidth=1.8` — fixed line width
- `sort=False` in `sns.lineplot` — lines connect data points in the order they appear in the DataFrame, not sorted by x-value
- `sns.color_palette("tab10"/"tab20")` — fixed palette choice
- Legend only on first subplot

### 9. Currently valid?
**Yes.** The plot renders correctly. However, the scientific interpretation is questionable due to the categorical x-axis.

### 10. Recommendation: **REDESIGN or REMOVE**
This figure has the weakest scientific justification of the 8. The trajectory format is inappropriate for categorical conditions. Consider:
- **Replace with a heatmap** (rows=seeds, columns=conditions, color=metric value) — this would show the same seed×condition pattern without implying continuity
- **Remove entirely** if F2 (box+strip) adequately covers the distribution information
- **Keep only if** the seed-level consistency story is explicitly part of the thesis narrative

---

## F7 — ECDF Overlay

**File:** `fig_ecdf_overlay.png`

### 1. Scientific purpose
Compare the cumulative distributions of seed-level metric values across conditions, revealing differences in location, spread, and shape without binning assumptions.

### 2. Exact input data
`seed_aggregated.csv` — one row per (condition × seed).

### 3. Unit of observation
Each ECDF curve = cumulative distribution of seed-level values for one condition.

### 4. Visual elements
- **ECDF curves:** One per condition, showing proportion of seed values ≤ x
- **Color:** Condition-colored curves (consistent palette)
- **Panels:** One per metric
- **X-axis:** Metric value
- **Y-axis:** Proportion (0–1)

### 5. Current strengths
- ECDF is a powerful non-parametric visualization that shows the full distribution shape
- No binning or smoothing assumptions
- Condition coloring is consistent with other figures
- Overlaid curves allow direct comparison

### 6. Current visual problems
- **With only 10 points per ECDF, the curves are jagged.** The step function has large jumps, making fine-grained comparison difficult.
- **8 overlapping curves can be hard to distinguish.** Especially in the middle of the distribution where curves cross.
- **The `husl` palette may not provide sufficient contrast** for 8 conditions, especially for colorblind readers.
- **No sample size annotation.** The reader cannot tell how many seeds contributed to each curve.

### 7. Potentially misleading elements
- With n=10, the ECDF is a very rough estimate of the true distribution. The curves might suggest more precision than the data supports.
- Crossing ECDF curves can be difficult to interpret — the reader must mentally track which line is which.

### 8. Hard-coded assumptions
- `linewidth=2.2` — fixed line width
- `husl` palette — fixed color scheme
- Condition ordering via `get_condition_order()` (POU-specific)
- Legend placed outside plot via `adjust_legend()`

### 9. Currently valid?
**Yes.** The plot renders correctly. The ECDF format is scientifically sound for comparing distributions.

### 10. Recommendation: **KEEP**
The ECDF overlay is a strong choice for distribution comparison. Minor improvements could include:
- Adding `n=XX` to the legend entries
- Using a colorblind-friendly palette (e.g., `colorblind` or `Set2`)
- Adding vertical reference lines for overall median or baseline median

---

## F8 — Metric Relationship Matrix (Scatter Matrix)

**File:** `fig_scatter_matrix.png`

### 1. Scientific purpose
Show pairwise relationships between key metrics, revealing correlations, clusters, and outliers across conditions.

### 2. Exact input data
`seed_aggregated.csv` — one row per (condition × seed). Uses a subset of metrics.

### 3. Unit of observation
Each point = one seed-level observation for one condition.

### 4. Visual elements
- **Scatter plots:** Pairwise metric relationships (lower triangle only via `corner=True`)
- **Color:** Condition-colored points
- **Diagonal:** Histograms/KDEs of each metric's marginal distribution
- **Metrics included:** `pLDDT_mean`, `pae_mean`, `contact_prob_mean`, `ranking_score`, plus first chain pLDDT

### 5. Current strengths
- `corner=True` eliminates redundant upper triangle
- Small point size (`s=20`) and alpha (`0.6`) help with overplotting
- Consistent condition coloring
- `sns.pairplot` is a well-tested, reliable function

### 6. Current visual problems
- **Hard-coded metric list doesn't match the schema.** The function uses its own `f8_metrics` list (`['pLDDT_mean', 'pae_mean', 'contact_prob_mean', 'ranking_score']`) instead of the schema's `metrics`. If the schema includes chain metrics, the first chain metric is appended, but this logic is separate from the schema.
- **chain_B and chain_C are mostly NaN.** Including them creates sparse panels.
- **With 8 conditions and 81 points total, dense regions overlap significantly.** The point colors blend together.
- **No correlation annotations.** The scatter matrix shows relationships but doesn't quantify them.
- **Metric labels on axes are raw column names** (e.g., "pLDDT_mean") rather than display labels, though `get_display_label` is applied via `rename_map`.

### 7. Potentially misleading elements
- Overlapping points in dense regions can make it appear that two conditions overlap when they actually occupy different regions of the metric space.
- The scatter matrix doesn't account for the hierarchical structure (seeds nested within conditions).

### 8. Hard-coded assumptions
- `f8_metrics = ['pLDDT_mean', 'pae_mean', 'contact_prob_mean', 'ranking_score']` — fixed metric list
- `height=2.2` — fixed panel size
- `plot_kws={'alpha': 0.6, 's': 20}` — fixed point appearance
- Only the first chain metric is added (if any)
- Title is hard-coded: "Figure 8: Metric Relationship Matrix"

### 9. Currently valid?
**Yes.** The plot renders correctly. The scatter matrix is a standard exploratory tool.

### 10. Recommendation: **MODIFY**
- Use the schema's metric list instead of the hard-coded `f8_metrics`
- Filter out chain metrics that are mostly NaN (>50% missing)
- Add Pearson correlation coefficients to each panel
- Consider using `sns.pairplot` with `diag_kind='kde'` for smoother marginal distributions
- Use a colorblind-friendly palette

---

## Cross-Cutting Issues

### Issue 1: Condition ordering is POU-specific
All 8 figures use `get_condition_order()` from `utils.py`, which contains a hard-coded `ptm_rank` dict for the 8 POU conditions. Any condition not in this dict gets rank 99 (sorted last). This works for the current dataset but breaks for any other experimental design.

**Recommendation:** Make condition ordering configurable via the schema, with the hard-coded dict as a fallback.

### Issue 2: Color palette is not colorblind-friendly
`get_condition_style_map()` uses `sns.color_palette("husl")`, which is not designed for colorblind accessibility. With 8 conditions, the palette may be difficult to distinguish for readers with deuteranopia or protanopia.

**Recommendation:** Use `sns.color_palette("colorblind")` or `sns.color_palette("Set2")` for better accessibility.

### Issue 3: No uncertainty visualization
None of the figures (except F4, which is dead) show uncertainty estimates. F2 shows distributions but not CIs. F3 shows point estimates without error bars. F5 shows SD without CI on the SD.

**Recommendation:** Add bootstrap CIs or SE bars where appropriate, especially in F3 (interaction) and F5 (variability).

### Issue 4: Figure naming convention
Figures are named `fig_*.png` with a hard-coded "Figure N:" title prefix. The numbering (1-8) is fixed regardless of which figures are actually generated.

**Recommendation:** Use a consistent naming scheme that reflects the content, not a sequential number.

### Issue 5: `descriptive_stats.csv` lacks `condition_name`
F3 and F5 depend on `condition_name` being present in `descriptive_stats.csv`, but the pipeline's `_stage_run_analysis` drops it during the groupby. The orchestrator patches this by joining from `seed_aggregated.csv`.

**Recommendation:** Include `condition_name` in the groupby output, or document the patch as a known dependency.

---

## Summary Table

| Figure | Purpose | Currently valid? | Recommendation | Priority |
|--------|---------|-----------------|----------------|----------|
| F1 | QC completeness | Partially (low info density) | KEEP (minor styling) | Low |
| F2 | Seed distributions | Yes | MODIFY (add n, fix chain panels) | Medium |
| F3 | Factorial interaction | Conditionally (POU-specific) | MODIFY (add error bars, categorical x) | High |
| F4 | Effect size forest | No (always empty) | KEEP (document as pending) | Low |
| F5 | Variability | Yes | MODIFY (add CI, consider heatmap) | Medium |
| F6 | Seed trajectories | Yes (but weak concept) | REDESIGN or REMOVE | High |
| F7 | ECDF overlay | Yes | KEEP (minor improvements) | Low |
| F8 | Scatter matrix | Yes | MODIFY (use schema, filter NaNs) | Medium |
