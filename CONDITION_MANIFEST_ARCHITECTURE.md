# Condition Manifest Architecture

## Overview

The **condition manifest** is a protein-agnostic registry system that serves as the authoritative description of all AF3 modeling conditions. It separates biological condition specification from AF3 JSON generation and from downstream analysis.

## Core Principle

```
Biological condition
        ↓
Condition metadata (manifest + junction tables)
        ↓
Component registries (proteins, modifications, nucleic acids, etc.)
        ↓
AF3 representation resolution
        ↓
JSON generation
```

**Adding a new protein is a data-entry operation, not a programming operation.**

## Directory Structure

```
condition_manifest/
├── __init__.py              # Package exports
├── registries.py            # Data classes and CSV loading for all registries
├── manifest.py              # Master condition manifest and junction tables
├── validation.py            # Generic validation system
├── inspection.py            # Design inspection and descriptive summaries
└── tests/
    └── test_manifest.py     # 46 tests covering all requirements

testdata/
├── pou2/registries/         # POU/OCT4 example dataset
│   ├── protein_registry.csv
│   ├── construct_registry.csv
│   ├── modification_registry.csv
│   ├── nucleic_acid_registry.csv
│   ├── ligand_registry.csv
│   ├── ion_registry.csv
│   ├── partner_registry.csv
│   ├── af3_compatibility_registry.csv
│   ├── residue_mapping_registry.csv
│   ├── master_condition_manifest.csv
│   ├── condition_modifications.csv
│   ├── condition_entities.csv
│   ├── condition_factors.csv
│   └── exclusivity_rules.csv
└── protein_x/registries/    # Second protein (proves generality)
    ├── ... (same CSV structure)
```

## Registry Files

### Protein Registry (`protein_registry.csv`)

| Field | Description |
|-------|-------------|
| `protein_id` | Primary key |
| `protein_name` | Human-readable name |
| `gene_name` | Gene symbol |
| `uniprot_id` | UniProt accession |
| `species` | Organism |
| `sequence_id` | Versioned sequence identifier |
| `sequence_version` | Sequence version |
| `notes` | Free text |

### Construct Registry (`construct_registry.csv`)

| Field | Description |
|-------|-------------|
| `construct_id` | Primary key |
| `protein_id` | Foreign key → protein_registry |
| `construct_name` | Human-readable name |
| `domain_id` | Domain identifier |
| `domain_name` | Domain name |
| `domain_start` | 1-based start position |
| `domain_end` | 1-based end position |
| `construct_sequence` | Amino acid sequence |
| `sequence_length` | Length of sequence |
| `source` | Reference source |
| `notes` | Free text |

### Modification Registry (`modification_registry.csv`)

| Field | Description |
|-------|-------------|
| `modification_id` | Primary key |
| `modification_name` | Human-readable name |
| `modification_class` | Category (phosphorylation, acetylation, etc.) |
| `base_residue` | Standard residue (S, T, K, etc.) |
| `modified_residue` | AF3 CCD code (SEP, TPO, ALY, etc.) |
| `chemical_description` | Chemical description |
| `molecular_formula` | Chemical formula |
| `formal_charge` | Charge |
| `evidence_level` | A_CONFIRMED_FUNCTIONAL, B_CONFIRMED_DETECTED, C_INDIRECT_OR_CONTEXTUAL, D_EXPLORATORY |
| `evidence_source` | Reference |
| `functional_evidence` | Functional characterization |
| `notes` | Free text |

### Nucleic Acid Registry (`nucleic_acid_registry.csv`)

Shared schema for DNA and RNA.

| Field | Description |
|-------|-------------|
| `entity_id` | Primary key |
| `entity_type` | "dna" or "rna" |
| `name` | Human-readable name |
| `sequence` | Nucleotide sequence |
| `sequence_version` | Version |
| `role` | Biological role |
| `source` | Reference source |
| `evidence_level` | Evidence classification |
| `notes` | Free text |

### Ligand Registry (`ligand_registry.csv`)

| Field | Description |
|-------|-------------|
| `ligand_id` | Primary key |
| `ligand_name` | Human-readable name |
| `ccd_code` | PDB CCD code |
| `smiles` | SMILES string |
| `inchi` | InChI string |
| `ligand_role` | Biological role |
| `source` | Reference source |
| `evidence_level` | Evidence classification |
| `custom_ccd_required` | "true"/"false" |
| `notes` | Free text |

### Ion Registry (`ion_registry.csv`)

| Field | Description |
|-------|-------------|
| `ion_id` | Primary key |
| `ion_name` | Human-readable name |
| `charge` | Ionic charge |
| `ccd_code` | PDB CCD code |
| `role` | Biological role |
| `concentration` | Experimental concentration |
| `notes` | Free text |

### Partner Registry (`partner_registry.csv`)

| Field | Description |
|-------|-------------|
| `partner_id` | Primary key |
| `partner_name` | Human-readable name |
| `protein_id` | Foreign key → protein_registry |
| `uniprot_id` | UniProt accession |
| `species` | Organism |
| `sequence_id` | Sequence identifier |
| `role` | Interaction role |
| `interaction_evidence` | Evidence type |
| `source` | Reference source |
| `notes` | Free text |

### AF3 Compatibility Registry (`af3_compatibility_registry.csv`)

Maps biological modifications to their AF3 representation.

| Field | Description |
|-------|-------------|
| `representation_id` | Primary key |
| `modification_id` | Foreign key → modification_registry |
| `entity_type` | "protein", "dna", "rna", "ligand" |
| `component_type` | AF3 representation type |
| `ccd_code` | PDB CCD code for AF3 |
| `custom_ccd_required` | "true"/"false" |
| `custom_ccd_id` | Custom CCD identifier |
| `covalent_bond_required` | "true"/"false" |
| `af3_status` | "supported", "unsupported", "experimental" |
| `verification_source` | How this was verified |
| `notes` | Free text |

### Residue Mapping Registry (`residue_mapping_registry.csv`)

Maps residue numbering between systems.

| Field | Description |
|-------|-------------|
| `mapping_id` | Primary key |
| `construct_id` | Foreign key → construct_registry |
| `sequence_position` | Full-length position (1-based) |
| `construct_position` | Construct/domain position (1-based) |
| `domain_position` | Domain position (1-based) |
| `reference_position` | Literature/database position |
| `reference_species` | Species for reference numbering |
| `reference_database` | Database name |
| `reference_accession` | Accession number |
| `numbering_mapping_status` | "exact", "offset", "gap", "uncertain" |
| `offset` | Integer offset if status is "offset" |
| `notes` | Free text |

## Master Manifest

### `master_condition_manifest.csv`

| Field | Description |
|-------|-------------|
| `condition_id` | Primary key |
| `condition_name` | Human-readable name |
| `condition_group` | Grouping (e.g., "baseline", "modified") |
| `parent_condition_id` | Reference to parent condition |
| `status` | "planned", "complete", "incomplete" |
| `description` | What this condition is |
| `biological_rationale` | Why this condition exists |
| `experimental_tier` | "primary", "secondary", "exploratory" |
| `experimental_priority` | Numeric priority |
| `notes` | Free text |

### `condition_modifications.csv` (junction table)

| Field | Description |
|-------|-------------|
| `condition_id` | FK → manifest |
| `modification_id` | FK → modification_registry |
| `sequence_position` | Position on the sequence |
| `construct_id` | FK → construct_registry |
| `stoichiometry` | Number of modifications |
| `evidence_level` | Evidence for this specific modification |
| `notes` | Free text |

### `condition_entities.csv` (junction table)

| Field | Description |
|-------|-------------|
| `condition_id` | FK → manifest |
| `entity_type` | "dna", "rna", "ligand", "ion", "partner" |
| `entity_id` | FK → appropriate registry |
| `stoichiometry` | Number of copies |
| `role` | Biological role in this condition |
| `notes` | Free text |

### `condition_factors.csv` (junction table)

| Field | Description |
|-------|-------------|
| `condition_id` | FK → manifest |
| `factor_name` | Experimental variable name |
| `factor_level` | Value of the variable |
| `factor_role` | "treatment", "control", "baseline" |
| `notes` | Free text |

### `exclusivity_rules.csv` (optional)

| Field | Description |
|-------|-------------|
| `rule_id` | Primary key |
| `exclusive_group` | Group identifier |
| `exclusive_key` | Composite constraint key |
| `entity_type` | What this applies to |
| `description` | Human-readable rule |
| `notes` | Free text |

## Usage

```python
from condition_manifest import (
    load_master_manifest,
    load_protein_registry,
    load_modification_registry,
    validate_manifest,
    inspect_manifest,
)

# Load registries
proteins = load_protein_registry(Path("testdata/pou2/registries/protein_registry.csv"))
mods = load_modification_registry(Path("testdata/pou2/registries/modification_registry.csv"))

# Load manifest
manifest = load_master_manifest(
    Path("testdata/pou2/registries/master_condition_manifest.csv"),
    modifications_path=Path("testdata/pou2/registries/condition_modifications.csv"),
    entities_path=Path("testdata/pou2/registries/condition_entities.csv"),
    factors_path=Path("testdata/pou2/registries/condition_factors.csv"),
)

# Validate
result = validate_manifest(manifest, protein_registry=proteins, modification_registry=mods)
if result.is_valid:
    print("Manifest is valid")
else:
    print(result.summary())

# Inspect design
insp = inspect_manifest(manifest)
print(insp.summary())

# Query
for cid in manifest.condition_ids:
    mods = manifest.get_modifications_for_condition(cid)
    entities = manifest.get_entities_for_condition(cid)
    factors = manifest.get_factors_for_condition(cid)
    print(f"{cid}: {len(mods)} mods, {len(entities)} entities, {len(factors)} factors")
```

## How to Add a New Protein

1. **Create a new directory** under `testdata/` (e.g., `testdata/my_protein/registries/`)
2. **Copy the CSV templates** from any existing dataset
3. **Fill in the registries** with your protein's information:
   - `protein_registry.csv` — your protein(s)
   - `construct_registry.csv` — your constructs/domains
   - `modification_registry.csv` — your modifications (PTMs, mutations, etc.)
   - `nucleic_acid_registry.csv` — any DNA/RNA involved
   - `ligand_registry.csv` — any ligands
   - `ion_registry.csv` — any ions
   - `partner_registry.csv` — any interacting proteins
   - `af3_compatibility_registry.csv` — how each modification maps to AF3
   - `residue_mapping_registry.csv` — residue numbering mappings
4. **Fill in the manifest**:
   - `master_condition_manifest.csv` — your experimental conditions
   - `condition_modifications.csv` — which modifications in which conditions
   - `condition_entities.csv` — which entities in which conditions
   - `condition_factors.csv` — experimental design factors
   - `exclusivity_rules.csv` — any mutual exclusivity constraints
5. **Run validation**:
   ```python
   from condition_manifest import load_master_manifest, validate_manifest
   manifest = load_master_manifest(Path("testdata/my_protein/registries/master_condition_manifest.csv"), ...)
   result = validate_manifest(manifest, ...)
   assert result.is_valid, result.summary()
   ```
6. **No code changes required.**

## Evidence Levels

| Level | Name | Description |
|-------|------|-------------|
| A | CONFIRMED_FUNCTIONAL | Directly demonstrated experimentally and functionally characterized |
| B | CONFIRMED_DETECTED | Experimentally detected but functional consequence unclear |
| C | INDIRECT_OR_CONTEXTUAL | Strong biological/contextual evidence but not directly demonstrated |
| D | EXPLORATORY | Computationally proposed or intentionally exploratory |

## Migration from Existing System

The existing `experiment_metadata.json` format (used by `af3_analysis/experiment_metadata.py`) can coexist with the manifest system. The manifest provides richer information but the simple JSON format remains useful for quick experiments.

To migrate:
1. The `condition_factors.csv` replaces the `attributes` section of `experiment_metadata.json`
2. The `master_condition_manifest.csv` replaces the `conditions` section
3. The registries provide entity information that was previously implicit

## Files NOT Modified

The following existing components were NOT modified:
- `af3_analysis/experiment_metadata.py` — remains available for backward compatibility
- `af3_analysis/visualization/` — continues to work with the existing schema system
- `af3_analysis/pipeline.py` — extraction pipeline unchanged
- `af3inputbuilder/af3_builder/` — builder unchanged
- `af3_analysis/schemas/` — existing schemas unchanged

## Generic Components

Everything in `condition_manifest/` is generic:
- No hard-coded protein names
- No hard-coded PTM names
- No hard-coded factor names
- No hard-coded experimental design assumptions
- No biological interpretation
- No statistical methodology selection

The only protein-specific content is in the **data files** (CSV registries), not in the code.
