# Structural/Geometric Analysis Subsystem — Architecture

> **Status:** Design document — not yet implemented.
>
> This document specifies the architecture for a generic structural/geometric
> analysis layer inside the `af3_analysis` package. It is designed for arbitrary
> AF3 projects and must not contain hard-coded biological assumptions.

---

## 1. Architecture Overview

The structural subsystem adds 3D geometry analysis to the existing confidence-metric
pipeline. It reads AF3 predicted structure files (mmCIF), extracts coordinates,
computes geometric metrics, compares structures across conditions, and produces
normalised output tables that the existing statistical and visualisation layers
consume.

```
AF3 mmCIF files
│
▼
┌─────────────────────────┐
│  io/structure_reader.py  │  Parse mmCIF → NormalisedStructure
└────────────┬────────────┘
             ▼
┌─────────────────────────┐
│  structural/representation.py  │  Canonical internal model
└────────────┬────────────┘
             ▼
┌─────────────────────────┐
│  structural/alignment.py     │  Pairwise / reference alignment
└────────────┬────────────┘
             ▼
┌─────────────────────────────────────┐
│  structural/metrics/                │  Pluggable metric implementations
│    rmsd.py                          │
│    centroid.py                      │
│    contacts.py                      │
│    distance_matrix.py               │
│    interface.py                     │
│    registry.py                      │
└────────────┬────────────────────────┘
             ▼
┌─────────────────────────┐
│  structural/comparison.py   │  Cross-condition structural comparisons
└────────────┬────────────┘
             ▼
┌─────────────────────────┐
│  structural/tables.py       │  Normalised output CSV/Parquet schemas
└────────────┬────────────┘
             ▼
     Existing statistical + visualisation layers
```

**Principle:** each layer depends only on layers above it, never below. The
output of the entire subsystem is normalised CSV/DataFrame tables that the
existing pipeline consumes through the same `tables/` directory interface used
by confidence metrics.

---

## 2. Module Structure

All new code lives under `af3_analysis/structural/`. The `af3_analysis/io/`
package gains one new file for CIF parsing.

```
af3_analysis/
├── io/
│   ├── raw_af3_reader.py          # existing — already discovers *.cif
│   ├── structure_reader.py        # NEW — parse mmCIF → NormalisedStructure
│   └── ...
├── structural/                    # NEW package
│   ├── __init__.py
│   ├── representation.py          # NormalisedStructure, Chain, Residue, Atom
│   ├── alignment.py               # Alignment engine + ComparisonResult
│   ├── comparison.py              # Cross-condition comparison orchestrator
│   ├── tables.py                  # Output table schemas
│   ├── config.py                  # StructuralConfig dataclass
│   ├── enums.py                   # Structural status enums
│   ├── metrics/                   # Pluggable metric sub-package
│   │   ├── __init__.py
│   │   ├── base.py                # StructuralMetric ABC
│   │   ├── rmsd.py                # Global + chain RMSD
│   │   ├── centroid.py            # Centroid displacement
│   │   ├── contacts.py            # Contact-map differences
│   │   ├── distance_matrix.py     # Pairwise Cα distance changes
│   │   ├── interface.py           # Inter-chain interface metrics
│   │   └── registry.py            # MetricRegistry for structural metrics
│   └── regions.py                 # Region definition + selection
└── ...
```

### 2.1 Responsibility Summary

| Module | Responsibility | No knowledge of |
|---|---|---|
| `io/structure_reader.py` | Parse mmCIF → `NormalisedStructure` | Conditions, biology |
| `representation.py` | Canonical data model | File format, conditions |
| `alignment.py` | Pairwise atom matching + superposition | Metric calculation |
| `metrics/base.py` | Metric ABC + applicability declaration | Other metrics |
| `metrics/*.py` | Individual metric implementations | Each other |
| `comparison.py` | Orchestrator: which pairs to compare | Individual metrics |
| `tables.py` | Output DataFrame contracts | Statistics |
| `config.py` | Structural analysis configuration | Any particular project |
| `regions.py` | User-supplied region definitions | Specific domain names |
| `enums.py` | Status codes | Everything else |

---

## 3. Data-Flow Diagram

```
                     ┌────────────────────────┐
                     │  Raw AF3 Root Directory │
                     │  <run>/cond/seed_s/*_model.cif │
                     └───────────┬────────────┘
                                 │
            ┌────────────────────┼────────────────────┐
            ▼                    ▼                    ▼
   ┌────────────────┐  ┌────────────────┐  ┌────────────────┐
   │ Replicate A,s1 │  │ Replicate A,s2 │  │ Replicate B,s1 │ ...
   └───────┬────────┘  └───────┬────────┘  └───────┬────────┘
           │                   │                   │
           ▼                   ▼                   ▼
   ┌──────────────────────────────────────────────────────────┐
   │  io/structure_reader.py  (mmCIF → NormalisedStructure)   │
   └───────────────────────────┬──────────────────────────────┘
                               │
                               ▼
   ┌──────────────────────────────────────────────────────────┐
   │  structural/representation.py                            │
   │  (one NormalisedStructure per CIF file)                  │
   └───────────────────────────┬──────────────────────────────┘
                               │
            ┌──────────────────┼──────────────────┐
            ▼                  ▼                  ▼
   ┌──────────────┐  ┌────────────────┐  ┌────────────────┐
   │   Within-    │  │  Within-cond   │  │  Between-cond  │
   │   condition  │  │  vs reference  │  │  vs reference  │
   │  variability │  │  (matched seed)│  │  (all pairs)   │
   └──────┬───────┘  └───────┬────────┘  └───────┬────────┘
          │                  │                   │
          ▼                  ▼                   ▼
   ┌──────────────────────────────────────────────────────────┐
   │  structural/tables.py                                     │
   │  (normalised CSVs → tables/ directory)                   │
   └───────────────────────────┬──────────────────────────────┘
                               │
                               ▼
   ┌──────────────────────────────────────────────────────────┐
   │  Existing statistical/ and visualization/ layers          │
   │  (consume structural tables like any other metric)       │
   └──────────────────────────────────────────────────────────┘
```

---

## 4. Normalised Structural Data Model

### 4.1 Core Records

All records follow the existing frozen-dataclass convention from
`schemas/records.py`.

```python
@dataclass(frozen=True)
class AtomRecord:
    """Single atom in a normalised structure."""
    atom_id: int                    # Original mmCIF _atom_site.id
    atom_name: str                  # e.g. "CA", "N", "O", "P"
    atom_type: str                  # e.g. "C", "N", "O", "S", "P"
    comp_id: str                    # Residue/monomer identifier (e.g. "ALA", "DA")
    chain_id: str                   # Original label_asym_id (e.g. "A", "B")
    entity_id: int                  # Original label_entity_id
    seq_id: Optional[int]           # label_seq_id (None = insertion code / missing)
    auth_seq_id: Optional[int]      # auth_seq_id
    auth_asym_id: str               # auth_asym_id
    coords: tuple[float, float, float]  # (x, y, z) Ångströms
    b_factor: Optional[float]       # B_iso_or_equiv (proxy for local confidence)
    occupancy: Optional[float]      # occupancy
    model_num: int                  # pdbx_PDB_model_num


@dataclass(frozen=True)
class ResidueRecord:
    """Residue-level summary (aggregated from atoms)."""
    comp_id: str
    chain_id: str
    entity_id: int
    seq_id: Optional[int]
    auth_seq_id: Optional[int]
    atom_names: tuple[str, ...]     # sorted tuple of atom names present
    n_atoms: int
    is_complete: bool               # True if all expected backbone atoms present


@dataclass(frozen=True)
class EntityInfo:
    """Entity metadata extracted from mmCIF _entity + _entity_poly."""
    entity_id: int
    chain_id: str                   # label_asym_id
    auth_asym_id: str
    entity_type: str                # "polymer", "non-polymer", "water", "branched"
    polymer_type: Optional[str]     # e.g. "polypeptide(L)", "polydeoxyribonucleotide", None
    description: Optional[str]
    n_residues: int
    n_atoms: int


@dataclass(frozen=True)
class ChainInfo:
    """Per-chain metadata."""
    chain_id: str
    auth_asym_id: str
    entity_id: int
    entity_type: str
    polymer_type: Optional[str]
    n_residues: int
    n_atoms: int
    first_seq_id: Optional[int]
    last_seq_id: Optional[int]
```

### 4.2 NormalisedStructure

```python
@dataclass
class NormalisedStructure:
    """
    Complete normalised representation of one AF3 prediction.

    One NormalisedStructure per (condition, seed, sample) CIF file.
    """
    # Identity
    source_path: Path
    source_checksum: str
    condition_id: str
    seed: int
    sample: int

    # Contents
    atoms: list[AtomRecord]
    residues: list[ResidueRecord]
    entities: list[EntityInfo]
    chains: list[ChainInfo]

    # Derived indices (built on construction)
    _by_chain: dict[str, list[AtomRecord]]     # chain_id → atoms
    _by_entity: dict[int, list[AtomRecord]]    # entity_id → atoms
    _by_seq: dict[tuple[str, int], list[AtomRecord]]  # (chain_id, seq_id) → atoms

    # Convenience queries
    @property
    def chain_ids(self) -> list[str]: ...
    @property
    def entity_ids(self) -> list[int]: ...
    @property
    def n_atoms(self) -> int: ...
    @property
    def n_chains(self) -> int: ...

    def get_chain(self, chain_id: str) -> list[AtomRecord]: ...
    def get_entity(self, entity_id: int) -> list[AtomRecord]: ...
    def get_ca_atoms(self, chain_ids: Optional[list[str]] = None) -> list[AtomRecord]: ...
    def get_backbone_atoms(self, chain_ids: Optional[list[str]] = None) -> list[AtomRecord]: ...
    def get_residue(self, chain_id: str, seq_id: int) -> list[AtomRecord]: ...
```

### 4.3 Design Rules

- **No chain-letter assumptions.** Chain IDs are opaque strings discovered from
  the mmCIF file.
- **No entity-type assumptions.** The model records the entity type as found
  (protein, DNA, RNA, ligand, water, etc.) but does not assume any particular
  type is present.
- **No numbering assumptions.** `seq_id` and `auth_seq_id` are preserved as
  nullable integers; gaps are represented explicitly.
- **Deterministic.** For the same CIF file the NormalisedStructure is identical.
- **Lazy derivatives.** The `_by_chain` etc. indices are built once and reused.

---

## 5. Comparison Model

### 5.1 Pairwise Alignment Result

```python
@dataclass(frozen=True)
class AlignmentResult:
    """Result of aligning two NormalisedStructures."""
    status: ComparisonStatus        # enum (see §11)
    reason: str                     # Human-readable explanation if not "comparable"

    # Atom pairing
    n_common_atoms: int
    paired_atom_ids: list[tuple[int, int]]  # (atom_id_A, atom_id_B)

    # Common residue info
    n_common_residues: int
    common_residues: list[tuple[str, int, int]]  # (chain_id, seq_A, seq_B)

    # Alignment transform (if applicable)
    rotation: Optional[np.ndarray]  # 3×3 rotation matrix
    translation: Optional[np.ndarray]  # 3-element translation vector
    rmsd_pre_alignment: Optional[float]
    rmsd_post_alignment: Optional[float]
```

### 5.2 Cross-Condition Comparison

```python
@dataclass(frozen=True)
class StructuralComparison:
    """One pairwise comparison between two predictions."""
    comparison_id: str
    prediction_a_id: str            # e.g. "condition_seed_sample"
    prediction_b_id: str
    condition_a: str
    condition_b: str
    seed_a: int
    seed_b: int

    alignment: AlignmentResult

    # Metric results (populated by the metrics layer)
    metrics: dict[str, Optional[float]]  # metric_id → value

    status: ComparisonStatus
    reason: str
```

### 5.3 Comparison Strategy

The comparison orchestrator builds a list of `(A, B)` pairs based on:

1. **Explicit reference** — user supplies `reference_condition`; all others
   compared against it.
2. **Matched-seed pairwise** — for each seed, compare all condition pairs at
   that seed. This respects the seed-matching design.
3. **All-pairs** — every prediction compared to every other prediction (for
   within-condition variability analysis).
4. **User-supplied pairs** — explicit list from configuration.

The default strategy is matched-seed comparison against an explicit reference.
The user configures which strategy is active.

---

## 6. Metric Interface

### 6.1 Base Class

```python
from abc import ABC, abstractmethod
from typing import Optional

class StructuralMetric(ABC):
    """
    Abstract base for a pluggable structural metric.

    Every metric declares what it requires (applicability) and computes
    a value for a single pair of aligned structures.
    """

    @property
    @abstractmethod
    def metric_id(self) -> str:
        """Unique identifier, e.g. 'rmsd_global_ca'."""
        ...

    @property
    @abstractmethod
    def display_name(self) -> str:
        """Human-readable name."""
        ...

    @property
    @abstractmethod
    def units(self) -> str:
        """Measurement units, e.g. 'Å', 'Å²', 'count'."""
        ...

    @property
    @abstractmethod
    def scope(self) -> str:
        """'global', 'chain', 'residue', 'interface', 'matrix_summary'."""
        ...

    @property
    @abstractmethod
    def direction(self) -> str:
        """'lower_better' or 'higher_better'."""
        ...

    @abstractmethod
    def requires(self) -> list[str]:
        """
        Declares what entity types or structural features are needed.

        Examples:
            []                           — always applicable
            ["protein"]                  — requires at least one protein chain
            ["protein", "protein"]       — requires at least two protein chains
            ["protein", "dna"]           — requires protein + nucleic acid
            ["protein", "ligand"]        — requires protein + ligand

        The comparison layer checks applicability before invoking the metric.
        """
        ...

    @abstractmethod
    def compute(
        self,
        structure_a: NormalisedStructure,
        structure_b: NormalisedStructure,
        alignment: AlignmentResult,
        *,
        region: Optional["RegionSpec"] = None,
        **kwargs,
    ) -> Optional[float]:
        """
        Compute the metric for one pair of aligned structures.

        Returns None if the metric is undefined for this pair
        (e.g. no common atoms in the specified region).
        """
        ...
```

### 6.2 Applicability Checking

The comparison layer evaluates applicability before running a metric:

```python
def is_applicable(
    metric: StructuralMetric,
    structure_a: NormalisedStructure,
    structure_b: NormalisedStructure,
) -> tuple[bool, str]:
    """Check whether a metric's requirements are met by two structures."""
    required = metric.requires()
    available_a = {e.entity_type for e in structure_a.entities}
    available_b = {e.entity_type for e in structure_b.entities}

    # For interface metrics: need at least N distinct chains
    # For entity-specific metrics: need the entity type in both structures
    ...
```

### 6.3 Initial Metric Set

| `metric_id` | Description | `requires()` | `scope` | `units` |
|---|---|---|---|---|
| `rmsd_global_ca` | Whole-structure Cα RMSD | `["protein"]` | `global` | `Å` |
| `rmsd_backbone` | Whole-structure backbone RMSD | `["protein"]` | `global` | `Å` |
| `rmsd_chain_{id}` | Per-chain Cα RMSD | `[]` | `chain` | `Å` |
| `centroid_displacement` | Cα centroid distance | `["protein"]` | `global` | `Å` |
| `centroid_displacement_{id}` | Per-chain centroid displacement | `[]` | `chain` | `Å` |
| `max_atom_displacement` | Maximum per-atom displacement | `[]` | `global` | `Å` |
| `mean_atom_displacement` | Mean per-atom displacement | `[]` | `global` | `Å` |
| `contact_map_diff` | Contact-map difference (Jaccard) | `["protein"]` | `global` | `dimensionless` |
| `interface_contacts` | Inter-chain contact count change | `["protein", "protein"]` | `interface` | `count` |
| `interface_distance` | Minimum inter-chain heavy-atom distance | `[]` | `interface` | `Å` |
| `pairwise_distance_change` | Mean Cα pairwise distance change | `["protein"]` | `global` | `Å` |
| `region_rmsd` | RMSD over user-specified region | `[]` | `chain` | `Å` |

Metrics are registered in `structural/metrics/registry.py` and loaded at
pipeline start. New metrics are added by subclassing `StructuralMetric` and
registering the instance.

---

## 7. Region Definitions

### 7.1 RegionSpec

```python
@dataclass(frozen=True)
class RegionSpec:
    """
    Generic structural region definition.

    A region selects a subset of atoms from a NormalisedStructure.
    All fields are optional; the intersection of specified filters is used.
    """
    label: Optional[str] = None         # Human-readable name (e.g. "DNA-binding domain")

    chain_ids: Optional[list[str]] = None       # Restrict to these chains
    entity_ids: Optional[list[int]] = None      # Restrict to these entities
    entity_types: Optional[list[str]] = None    # Restrict to these entity types

    seq_range: Optional[tuple[int, int]] = None # Inclusive residue range (auth_seq_id)
    seq_list: Optional[list[int]] = None        # Explicit list of auth_seq_ids

    atom_names: Optional[list[str]] = None      # Restrict to these atom names

    sequence_pattern: Optional[str] = None      # future: sequence-based selection

    def select(self, structure: NormalisedStructure) -> list[AtomRecord]:
        """Return atoms matching this region in the given structure."""
        ...
```

### 7.2 Region Sources

Regions may be supplied from:

1. **Configuration** — static region definitions in the analysis config JSON.
2. **Metadata** — entity-type-based regions auto-generated from the mmCIF
   (e.g. "all protein chains", "all DNA chains").
3. **Computation** — detected regions (e.g. interface residues within a cutoff).
4. **Future: domain annotations** — user-supplied domain boundaries.

The system never assumes that regions exist. If no regions are supplied, metrics
operate on all atoms or all atoms of the required type.

### 7.3 Built-in Region Generators

```python
def protein_chains(structure: NormalisedStructure) -> RegionSpec:
    """Select all protein-type chains."""

def nucleic_acid_chains(structure: NormalisedStructure) -> RegionSpec:
    """Select all nucleic-acid-type chains."""

def interface_region(
    structure_a: NormalisedStructure,
    structure_b: NormalisedStructure,
    chain_a: str,
    chain_b: str,
    cutoff_angstrom: float = 8.0,
) -> RegionSpec:
    """Detect interface residues within cutoff."""
```

---

## 8. Reference Strategy

### 8.1 ReferenceConfig

```python
@dataclass(frozen=True)
class ReferenceConfig:
    """Configuration for how reference structures are selected."""
    strategy: str  # one of the strategies below

    # For explicit_reference:
    reference_condition: Optional[str] = None

    # For explicit_structure:
    reference_path: Optional[Path] = None

    # For matched_seed:
    # (compared to the reference condition at the same seed)
    seed_matching: bool = True
```

### 8.2 Strategy Catalogue

| `strategy` | Description | Requirements |
|---|---|---|
| `explicit_reference` | All conditions compared to one named reference condition | `reference_condition` set |
| `explicit_structure` | All predictions compared to one external CIF file | `reference_path` set |
| `first_condition` | First condition alphabetically used as reference | None |
| `pairwise_all` | Every condition pair compared | None |
| `pooled_reference` | Mean structure of the reference condition used | Reference condition has ≥ 1 prediction |

### 8.3 Rules

- The reference must be explicitly configured or the system must raise an error.
- The default is `explicit_reference` with `reference_condition=None`, which
  means "first condition alphabetically" if none is supplied (and a warning is
  logged).
- The reference structure is *not* privileged geometrically; it is simply the
  comparison anchor.
- Confidence metrics of the reference are not treated as "truth."

---

## 9. Seed Strategy

### 9.1 Principles

AF3 seeds are repeated stochastic predictions, not independent biological
replicates. The structural subsystem must:

- Preserve seed identity through all stages.
- Distinguish within-condition structural variability (across seeds) from
  between-condition structural differences.
- Support matched-seed comparisons (same seed, different conditions).
- Report per-seed metric values before any aggregation.

### 9.2 Comparison Modes

| Mode | Description | Pair count |
|---|---|---|
| `matched_seed` | For each seed, compare condition X vs reference at that seed | `n_seeds × (n_conditions - 1)` |
| `all_vs_all` | Every prediction compared to every other prediction | `n_predictions² / 2` |
| `seed_paired` | Only compare predictions sharing the same seed | Filtered all_vs_all |

### 9.3 Within-Condition Variability

For each condition, the subsystem computes structural variability metrics:

- **Pairwise RMSD between seeds** — how much does the predicted structure
  vary across seeds?
- **Mean pairwise RMSD** — summary statistic.
- **Structural ensemble spread** — range of RMSD values.

These are reported as structural-quality metrics, not biological measurements.

### 9.4 Seed Preservation

Every output table row includes `seed` as an identifier column. Aggregation
to condition-level summaries follows the same seed-first pattern used by the
confidence-metric pipeline (`preprocessing/aggregation.py`).

---

## 10. Statistical Layer

The structural subsystem produces metric values. Statistical comparison of
those values follows the existing pipeline philosophy:

### 10.1 Per-Metric Report

For each structural metric and comparison pair:

| Field | Description |
|---|---|
| `n` | Number of observations (seed-level pairs) |
| `mean` | Mean metric value across seed pairs |
| `median` | Median metric value |
| `sd` | Standard deviation |
| `effect_size` | Hedges' g or Cohen's d (using existing `statistical/comparisons.py`) |
| `ci_lower` | Lower bound of confidence interval |
| `ci_upper` | Upper bound of confidence interval |
| `direction_consistency` | Fraction of seed pairs where difference is in the reported direction |
| `n_comparable` | Number of seed pairs that were successfully compared |
| `n_incomparable` | Number of seed pairs that could not be compared |

### 10.2 Integration with Existing Statistics

- The structural subsystem writes metric values to `structural_metrics.csv`.
- The existing `preprocessing/aggregation.py` aggregates seed-level structural
  metrics identically to confidence metrics.
- The existing `statistical/comparisons.py` compares conditions on structural
  metrics using the same permutation/bootstrap framework.
- No new statistical methods are introduced.

### 10.3 Important Constraints

- Structural metric values are geometric measurements, not confidence scores.
- An RMSD difference does not imply biological activation.
- An effect size for RMSD is a geometric observation, not a biological claim.
- Statistical significance of structural differences is a reportable quantity,
  not a conclusion.

---

## 11. Missingness and Invalid Comparisons

### 11.1 ComparisonStatus Enum

```python
class ComparisonStatus(Enum):
    COMPARABLE = "comparable"
    PARTIALLY_COMPARABLE = "partially_comparable"
    NOT_COMPARABLE = "not_comparable"
```

### 11.2 Detailed Status Codes

| Status | Meaning | Action |
|---|---|---|
| `valid` | Both structures parsed and aligned successfully | Proceed |
| `insufficient_atoms` | Too few atoms for alignment or metric | Skip, log |
| `missing_chain` | Expected chain not found in one structure | Skip, log |
| `sequence_mismatch` | Chains have different sequences at aligned positions | Warn, proceed with subset |
| `no_common_residues` | Zero common residues between two structures | Skip, log |
| `missing_reference` | Reference structure not available | Abort comparison |
| `unsupported_entity_type` | Required entity type not present | Skip metric, log |
| `insufficient_replicates` | Fewer than minimum seed pairs | Report, flag |
| `parse_error` | CIF file could not be parsed | Skip structure, log |
| `model_number_mismatch` | Multiple models in CIF; selection ambiguous | Use first model, warn |
| `empty_structure` | CIF parsed but zero atoms extracted | Skip, log |

### 11.3 Rules

- **No silent imputation.** If a structure is missing, the comparison status is
  `not_comparable` with reason `missing_reference` or `parse_error`.
- **No silent dropping.** Every skipped comparison appears in
  `structural_comparisons.csv` with its status and reason.
- **Machine-readable.** Every output row has a `status` and `reason` column.
- **Auditable.** The count of incomparable pairs is reported in the run manifest.

---

## 12. Confidence vs Geometry

### 12.1 Separation

| Concept | Source | Meaning |
|---|---|---|
| pLDDT | `_confidences.json` | Per-atom predicted confidence |
| PAE | `_confidences.json` | Predicted aligned error between residue pairs |
| RMSD | mmCIF coordinates | Geometric distance between two predicted structures |
| Centroid displacement | mmCIF coordinates | Distance between structural centres |
| Contact map difference | mmCIF coordinates | Difference in residue-residue contacts |

The structural subsystem computes only the geometric columns. Confidence
metrics are already computed by the existing pipeline.

### 12.2 Optional Combination

At the reporting stage (not within the metrics layer), it may be useful to
correlate:

- pLDDT with per-residue displacement
- PAE with interface distance changes
- pLDDT with RMSD

These are **reporting-stage observations**, not metric definitions. They belong
in the visualisation or reporting layer, not in the metrics computation.

### 12.3 Cross-Reference in Output Tables

The `structural_metrics.csv` may include a `confidence_available` flag per
prediction, linking to the existing confidence data. But structural metric
values must never be filled from or modified by confidence values.

---

## 13. Output Contracts

All output tables are written to `<run>/tables/` using the same conventions as
the existing pipeline.

### 13.1 `structural_predictions.csv`

One row per successfully parsed CIF file.

| Column | Type | Description |
|---|---|---|
| `prediction_id` | str | Unique prediction identifier |
| `condition_id` | str | Condition identifier |
| `condition_name` | str | Human-readable condition name |
| `seed` | int | AF3 seed |
| `sample` | int | AF3 sample/model |
| `source_path` | str | Path to CIF file |
| `source_checksum` | str | SHA-256 of CIF file |
| `n_atoms` | int | Total atom count |
| `n_chains` | int | Number of chains |
| `chain_ids` | str | Semicolon-separated chain IDs |
| `entity_types` | str | Semicolon-separated entity types |
| `parse_status` | str | `success` or error code |
| `parse_reason` | str | Explanation if parse failed |

### 13.2 `structural_metrics.csv`

One row per prediction × metric × scope.

| Column | Type | Description |
|---|---|---|
| `prediction_id` | str | Prediction identifier |
| `condition_id` | str | Condition identifier |
| `condition_name` | str | Condition name |
| `seed` | int | Seed |
| `metric_id` | str | Metric identifier |
| `scope_type` | str | `global`, `chain`, `interface` |
| `scope_id` | str | Chain ID or empty for global |
| `value` | float | Metric value (NaN if undefined) |
| `status` | str | `present`, `undefined`, `not_applicable` |
| `reason` | str | Explanation if undefined |

### 13.3 `structural_comparisons.csv`

One row per comparison pair × metric.

| Column | Type | Description |
|---|---|---|
| `comparison_id` | str | Unique comparison identifier |
| `prediction_a_id` | str | Prediction A identifier |
| `prediction_b_id` | str | Prediction B identifier |
| `condition_a` | str | Condition A |
| `condition_b` | str | Condition B |
| `seed_a` | int | Seed of prediction A |
| `seed_b` | int | Seed of prediction B |
| `metric_id` | str | Metric identifier |
| `scope_type` | str | Scope |
| `scope_id` | str | Scope identifier |
| `value` | float | Metric value |
| `alignment_status` | str | Alignment status code |
| `alignment_reason` | str | Explanation |
| `n_common_atoms` | int | Number of common atoms used |

### 13.4 `structural_effects.csv`

One row per metric × comparison-condition-pair (seed-aggregated).

| Column | Type | Description |
|---|---|---|
| `metric_id` | str | Metric identifier |
| `scope_type` | str | Scope |
| `scope_id` | str | Scope identifier |
| `condition_a` | str | Test condition |
| `condition_b` | str | Reference condition |
| `n` | int | Number of seed-pair observations |
| `mean` | float | Mean metric difference |
| `median` | float | Median metric difference |
| `sd` | float | Standard deviation |
| `effect_size` | float | Standardised effect (Hedges' g) |
| `ci_lower` | float | Confidence interval lower bound |
| `ci_upper` | float | Confidence interval upper bound |
| `direction_consistency` | float | Fraction in reported direction |
| `n_comparable` | int | Number of comparable seed pairs |
| `n_incomparable` | int | Number of incomparable seed pairs |

### 13.5 `structural_quality.csv`

One row per condition × structural-quality metric.

| Column | Type | Description |
|---|---|---|
| `condition_id` | str | Condition |
| `quality_metric` | str | e.g. `mean_pairwise_seed_rmsd` |
| `value` | float | Quality metric value |
| `n_predictions` | int | Number of parsed predictions |
| `n_failed` | int | Number of failed parses |

### 13.6 Table Conventions

- All tables use UTF-8 encoding.
- `NaN` is used for missing numeric values.
- Empty strings are used for missing string values.
- Primary keys are composite (see column descriptions).
- All tables can be imported into Parquet via the existing `ParquetStore`.

---

## 14. Configuration Design

### 14.1 StructuralConfig

```python
@dataclass(frozen=True)
class StructuralConfig:
    """
    Configuration for structural analysis.

    All fields have sensible defaults. No project-specific defaults.
    """
    # Enable/disable
    enabled: bool = False

    # Reference
    reference_strategy: str = "explicit_reference"
    reference_condition: Optional[str] = None
    reference_path: Optional[Path] = None

    # Atom selection
    atom_selection: str = "ca"           # "ca", "backbone", "all_heavy", "all"
    min_common_atoms: int = 10           # Minimum atoms for valid alignment

    # Alignment
    alignment_method: str = "kabsch"     # "kabsch", "iterative"
    min_sequence_identity: float = 0.5   # Minimum fraction of common residues

    # Contact detection
    contact_cutoff_angstrom: float = 8.0

    # Region definitions (optional)
    regions: list[RegionSpec] = field(default_factory=list)

    # Metrics to compute (empty = all applicable)
    enabled_metrics: list[str] = field(default_factory=list)

    # Missing-data policy
    drop_incomparable: bool = False      # False = include with status

    # Output settings
    output_csv: bool = True
    output_parquet: bool = True

    # Seed settings
    comparison_mode: str = "matched_seed"  # "matched_seed", "all_vs_all", "seed_paired"
    min_seeds: int = 2                     # Minimum seeds for condition-level summary
```

### 14.2 Integration with AnalysisConfig

The existing `AnalysisConfig` already has `coordinate_analysis_enabled: bool`.
This flag would gate the structural subsystem. The `StructuralConfig` would be
an optional nested config within `AnalysisConfig`:

```python
@dataclass(frozen=True)
class AnalysisConfig:
    ...
    coordinate_analysis_enabled: bool = False
    structural_config: Optional[StructuralConfig] = None
    ...
```

### 14.3 Defaults

All defaults are project-agnostic:

- `atom_selection = "ca"` — standard for protein structural comparison.
- `contact_cutoff_angstrom = 8.0` — standard heavy-atom contact cutoff.
- `reference_strategy = "explicit_reference"` — user must specify.
- `enabled_metrics = []` — compute all applicable metrics.
- No default domain boundaries, no default residue ranges, no default chains.

---

## 15. Integration Points with Existing Pipeline

### 15.1 New Files

| File | Purpose |
|---|---|
| `af3_analysis/io/structure_reader.py` | mmCIF parsing |
| `af3_analysis/structural/__init__.py` | Package init |
| `af3_analysis/structural/representation.py` | Data model |
| `af3_analysis/structural/alignment.py` | Alignment engine |
| `af3_analysis/structural/comparison.py` | Comparison orchestrator |
| `af3_analysis/structural/tables.py` | Output schemas |
| `af3_analysis/structural/config.py` | StructuralConfig |
| `af3_analysis/structural/enums.py` | Status enums |
| `af3_analysis/structural/regions.py` | Region definitions |
| `af3_analysis/structural/metrics/__init__.py` | Metrics sub-package |
| `af3_analysis/structural/metrics/base.py` | Metric ABC |
| `af3_analysis/structural/metrics/rmsd.py` | RMSD metrics |
| `af3_analysis/structural/metrics/centroid.py` | Centroid metrics |
| `af3_analysis/structural/metrics/contacts.py` | Contact-map metrics |
| `af3_analysis/structural/metrics/distance_matrix.py` | Distance metrics |
| `af3_analysis/structural/metrics/interface.py` | Interface metrics |
| `af3_analysis/structural/metrics/registry.py` | Metric registry |

### 15.2 Existing Files Requiring Modification

| File | Change | Risk |
|---|---|---|
| `af3_analysis/config.py` | Add `structural_config` field to `AnalysisConfig`; extend `from_dict` and `to_dict` | Low — additive field |
| `af3_analysis/pipeline.py` | Add structural stage between stage 4 (analysis) and stage 5 (figures) | Low — new stage, no existing stage modified |
| `af3_analysis/io/raw_af3_reader.py` | No change — already discovers `*.cif` files and records `source_model_path` on `Replicate` | None |

### 15.3 Files That Must Remain Untouched

| File | Reason |
|---|---|
| `af3_analysis/schemas/records.py` | Existing confidence-metric records are independent |
| `af3_analysis/schemas/enums.py` | Confidence enums are independent |
| `af3_analysis/schemas/tables.py` | Confidence table schemas are independent |
| `af3_analysis/preprocessing/*` | Preprocessing pipeline for confidence metrics |
| `af3_analysis/statistical/*` | Statistical methods are reused, not modified |
| `af3_analysis/exploratory/*` | Exploratory analysis for confidence metrics |
| `af3_analysis/visualization/*` | Visualisation is consumed at the end, not modified |
| `af3_analysis/experiment_metadata.py` | Experiment metadata is independent |
| `af3inputbuilder/*` | Input builder is upstream, not downstream |

---

## 16. Testing Strategy

### 16.1 Unit Test Categories

All tests go in `af3_analysis/tests/unit/` following existing conventions.

#### Coordinate Parsing

| Test | Description |
|---|---|
| `test_parse_single_chain_protein` | Parse a single-chain protein CIF |
| `test_parse_multi_chain_protein_dna` | Parse a protein + DNA CIF |
| `test_parse_modified_residues` | Parse CIF with phosphorylated residues |
| `test_parse_ligand_atoms` | Parse CIF with non-polymer entities |
| `test_parse_ions` | Parse CIF with ion entities |
| `test_parse_missing_residues` | Handle insertion codes / missing seq_ids |
| `test_parse_multiple_models` | Handle multi-model CIF (use first model) |
| `test_parse_malformed_cif` | Graceful error on corrupted CIF |
| `test_parse_empty_cif` | Graceful error on zero-atom CIF |
| `test_source_checksum_recorded` | Verify provenance checksum |

#### Chain and Entity Extraction

| Test | Description |
|---|---|
| `test_entity_type_detection` | Correctly identify protein/DNA/RNA/ligand |
| `test_chain_id_preservation` | Chain IDs preserved exactly as in CIF |
| `test_polymer_type_recorded` | Polymer type stored from `_entity_poly` |
| `test_residue_count_per_chain` | Correct residue counts per chain |

#### Alignment

| Test | Description |
|---|---|
| `test_identical_structures_zero_rmsd` | Aligning identical coordinates → RMSD = 0 |
| `test_known_translation_rmsd` | Apply known translation, verify RMSD |
| `test_known_rotation_rmsd` | Apply known rotation, verify RMSD |
| `test_partial_overlap_alignment` | Different sequence lengths, partial alignment |
| `test_no_common_residues` | Returns `not_comparable` with reason |
| `test_missing_chain_in_one_structure` | Returns `partially_comparable` |
| `test_different_entity_composition` | Protein-only vs protein+DNA |
| `test_calpha_only_alignment` | Cα-only alignment option |
| `test_backbone_alignment` | Backbone atom alignment option |

#### Metrics

| Test | Description |
|---|---|
| `test_rmsd_identical` | RMSD of identical structures = 0 |
| `test_rmsd_known_value` | RMSD matches known analytical value |
| `test_centroid_zero_displacement` | Centroid displacement of identical structures = 0 |
| `test_centroid_known_displacement` | Known translation → known centroid distance |
| `test_contact_map_identical` | Identical structures → contact map diff = 0 |
| `test_contact_map_single_residue_change` | Single residue move → correct contact changes |
| `test_distance_matrix_known` | Pairwise distances match analytical values |
| `test_metric_applicability_protein_only` | Protein metric returns not_applicable for DNA-only |
| `test_metric_applicability_interface` | Interface metric requires ≥ 2 chains |

#### Comparison

| Test | Description |
|---|---|
| `test_matched_seed_comparison` | Correct pairing for matched-seed mode |
| `test_explicit_reference` | All conditions compared to named reference |
| `test_pairwise_all` | All condition pairs generated |
| `test_incomparable_pair_logged` | Incomparable pairs appear in output with status |
| `test_cross_condition_same_seed` | Compare two conditions at same seed |

#### Seed Handling

| Test | Description |
|---|---|
| `test_seed_identity_preserved` | Seed column present in all outputs |
| `test_within_condition_variability` | Pairwise seed RMSD computed correctly |
| `test_between_condition_difference` | Cross-condition metric differences computed |

#### Reference Selection

| Test | Description |
|---|---|
| `test_first_condition_reference` | Default reference when none specified |
| `test_explicit_reference_condition` | Named reference used correctly |
| `test_external_reference_structure` | External CIF used as reference |
| `test_missing_reference_error` | Error raised when reference not found |

### 16.2 Test Data Strategy

- Use **synthetic CIF files** generated programmatically for unit tests.
- Use **existing testdata CIF files** for integration tests only.
- Never hard-code specific biological assertions in structural tests.
- Test fixtures create NormalisedStructure objects with known coordinates.

---

## 17. Migration and Backward Compatibility

### 17.1 Backward Compatibility

- The structural subsystem is **additive**. It adds new stages and new tables.
- Existing tables (`metrics_replicates.csv`, `seed_aggregated.csv`, etc.) are
  not modified.
- Existing pipeline stages continue to work unchanged.
- The `coordinate_analysis_enabled` flag defaults to `False`, so existing
  analyses are unaffected.

### 17.2 Migration Path

1. **Phase 1 (current design):** Architecture document. No code.
2. **Phase 2:** Implement `io/structure_reader.py` + `representation.py` +
   `alignment.py` with unit tests.
3. **Phase 3:** Implement initial metrics (RMSD, centroid, contacts) with tests.
4. **Phase 4:** Implement `comparison.py` + `tables.py` + pipeline integration.
5. **Phase 5:** Implement configuration integration + end-to-end tests.
6. **Phase 6:** Add visualisation hooks.

### 17.3 Dependency Requirements

| Package | Status | Purpose |
|---|---|---|
| `gemmi` | Already in requirements | mmCIF parsing |
| `numpy` | Already used | Coordinate arrays, linear algebra |
| `scipy` | Already used | Spatial distance computations |
| `pandas` | Already used | Output tables |

No new external dependencies are required. `gemmi` is already used by
`af3inputbuilder/af3_builder/utils/cif_slicer.py` and is listed in the
repository's dependencies.

---

## 17. Generic Examples

### Example A: Single Protein, 3 Conditions, 10 Seeds

```
Project: Kinase activation study
Conditions: wild_type, mutant_A, mutant_B
Entity: protein (chain A, 320 residues)
Seeds: 10 per condition, 1 sample per seed

Expected outputs:
- structural_predictions.csv: 30 rows (3 conditions × 10 seeds)
- structural_metrics.csv: 30 × n_metrics rows
- structural_comparisons.csv: 20 × n_metrics rows (mutant_A vs wt, mutant_B vs wt, 10 seeds each)
- structural_effects.csv: 2 × n_metrics rows (one per mutant-vs-wt comparison)
```

### Example B: Protein + DNA, 2 Conditions

```
Project: Transcription factor binding study
Conditions: apo, dna_bound
Entities: protein (chain A), DNA strand 1 (chain B), DNA strand 2 (chain C)
Seeds: 5 per condition, 3 samples per seed

Expected outputs:
- structural_predictions.csv: 30 rows (2 conditions × 5 seeds × 3 samples)
- structural_metrics.csv: includes interface metrics (protein-DNA contacts)
- structural_comparisons.csv: 15 × n_metrics rows (matched-seed, apo vs dna_bound)
```

### Example C: Multi-Chain Complex

```
Project: Protein-protein interaction study
Conditions: alone, partner_X, partner_Y, partner_X_Y
Entities vary:
  - alone: protein A (chain A)
  - partner_X: protein A (chain A) + protein X (chain B)
  - partner_Y: protein A (chain A) + protein Y (chain C)
  - partner_X_Y: protein A (chain A) + protein X (chain B) + protein Y (chain C)
Seeds: 10 per condition

Expected outputs:
- structural_predictions.csv: 40 rows
- structural_metrics.csv: interface metrics applicable only where ≥ 2 chains
- structural_comparisons.csv: status = "not_comparable" for chains that differ
  between conditions (e.g. comparing chain B in partner_X vs chain C in partner_Y)
```

### Example D: Modified Protein

```
Project: PTM structural effects
Conditions: unmodified, phosphorylated, methylated
Entities: protein (chain A, 150 residues)
Seeds: 10 per condition

Expected outputs:
- structural_predictions.csv: 30 rows
- structural_metrics.csv: per-residue metrics may differ at modification sites
- structural_comparisons.csv: full comparisons
- regions defined: user may specify modification-site regions for targeted RMSD
```

---

## Appendix: mmCIF Field Mapping

How AF3 mmCIF fields map to the normalised data model:

| mmCIF field | NormalisedStructure field | Notes |
|---|---|---|
| `_atom_site.group_PDB` | — | Filtered: only ATOM/HETATM records |
| `_atom_site.id` | `AtomRecord.atom_id` | Preserved as-is |
| `_atom_site.type_symbol` | `AtomRecord.atom_type` | Element symbol |
| `_atom_site.label_atom_id` | `AtomRecord.atom_name` | Atom name within residue |
| `_atom_site.label_comp_id` | `AtomRecord.comp_id` | Residue/monomer name |
| `_atom_site.label_asym_id` | `AtomRecord.chain_id` | Chain identifier |
| `_atom_site.label_entity_id` | `AtomRecord.entity_id` | Entity identifier |
| `_atom_site.label_seq_id` | `AtomRecord.seq_id` | Sequence position (nullable) |
| `_atom_site.Cartn_x/y/z` | `AtomRecord.coords` | 3-tuple of floats |
| `_atom_site.B_iso_or_equiv` | `AtomRecord.b_factor` | B-factor / pLDDT proxy |
| `_atom_site.occupancy` | `AtomRecord.occupancy` | Occupancy |
| `_atom_site.auth_seq_id` | `AtomRecord.auth_seq_id` | Author sequence numbering |
| `_atom_site.auth_asym_id` | `AtomRecord.auth_asym_id` | Author chain ID |
| `_atom_site.pdbx_PDB_model_num` | `AtomRecord.model_num` | Model number |
| `_entity.id` | `EntityInfo.entity_id` | |
| `_entity.type` | `EntityInfo.entity_type` | |
| `_entity_poly.entity_id` | — | Used to set `polymer_type` |
| `_entity_poly.pdbx_strand_id` | `EntityInfo.chain_id` | |
| `_entity_poly.type` | `EntityInfo.polymer_type` | |
| `_struct_asym.entity_id` | — | Maps chain → entity |
| `_struct_asym.id` | `EntityInfo.chain_id` | |

---

*End of architecture document.*
