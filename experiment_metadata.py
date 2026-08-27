"""
Experiment design metadata.

Loads, validates, and inspects experiment-design metadata from a JSON file.
The metadata describes conditions and their attributes without assuming any
specific biological entities.

Design principles
-----------------
* **No assumptions** about attribute names, types, or counts.
* **No inference** from condition-name strings.
* **Explicit** representation of every condition's attributes.
* **Validated** on load — unknown conditions, duplicates, and missing
  definitions produce explicit errors.
* **Descriptive only** — describes the design without selecting
  statistical models.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AttributeDefinition:
    """Definition of a single experimental attribute."""

    name: str
    attr_type: str  # "binary" or "categorical"
    description: str = ""

    def __post_init__(self) -> None:
        if self.attr_type not in ("binary", "categorical"):
            raise ValueError(
                f"Attribute '{self.name}' has invalid type '{self.attr_type}'; "
                f"expected 'binary' or 'categorical'."
            )


@dataclass(frozen=True)
class ConditionDefinition:
    """Definition of a single experimental condition."""

    condition_name: str
    label: str
    attributes: Dict[str, Any]  # attr_name -> value (bool for binary, str for categorical)


@dataclass
class ExperimentDesign:
    """Complete experiment design metadata.

    Attributes
    ----------
    experiment_id : str
        Unique identifier for the experiment.
    description : str
        Human-readable description.
    attributes : dict
        Mapping attribute name -> AttributeDefinition.
    conditions : dict
        Mapping condition name -> ConditionDefinition.
    """

    experiment_id: str
    description: str
    attributes: Dict[str, AttributeDefinition]
    conditions: Dict[str, ConditionDefinition]

    # -- Query helpers -------------------------------------------------------

    @property
    def condition_names(self) -> List[str]:
        """Sorted list of condition names."""
        return sorted(self.conditions.keys())

    @property
    def attribute_names(self) -> List[str]:
        """Sorted list of attribute names."""
        return sorted(self.attributes.keys())

    @property
    def n_conditions(self) -> int:
        return len(self.conditions)

    @property
    def n_attributes(self) -> int:
        return len(self.attributes)

    def get_attribute_values(self, attr_name: str) -> Set[Any]:
        """Return the set of observed values for an attribute across all conditions."""
        if attr_name not in self.attributes:
            raise KeyError(f"Attribute '{attr_name}' not defined in metadata.")
        return {
            cond.attributes.get(attr_name)
            for cond in self.conditions.values()
            if attr_name in cond.attributes
        }

    def get_attribute_type(self, attr_name: str) -> str:
        """Return the type of an attribute ('binary' or 'categorical')."""
        if attr_name not in self.attributes:
            raise KeyError(f"Attribute '{attr_name}' not defined in metadata.")
        return self.attributes[attr_name].attr_type

    def get_conditions_with_attribute(
        self, attr_name: str, value: Any = True
    ) -> List[str]:
        """Return condition names where attr_name == value."""
        return [
            name
            for name, cond in self.conditions.items()
            if cond.attributes.get(attr_name) == value
        ]

    def get_condition_label(self, condition_name: str) -> str:
        """Return the display label for a condition."""
        if condition_name not in self.conditions:
            return condition_name.replace("_", " ").title()
        return self.conditions[condition_name].label

    def get_attribute_label(self, attr_name: str) -> str:
        """Return a human-readable label for an attribute."""
        if attr_name in self.attributes:
            return self.attributes[attr_name].name
        return attr_name.replace("_", " ").title()

    def condition_attributes_matrix(self) -> Dict[str, Dict[str, Any]]:
        """Return a dict-of-dicts: condition_name -> {attr_name: value}."""
        return {
            name: dict(cond.attributes)
            for name, cond in self.conditions.items()
        }


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_experiment_design(path: str | Path) -> ExperimentDesign:
    """Load experiment design metadata from a JSON file.

    Parameters
    ----------
    path : str or Path
        Path to the JSON metadata file.

    Returns
    -------
    ExperimentDesign
        Validated experiment design.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file contains invalid or inconsistent metadata.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Experiment metadata file not found: {path}")

    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    return _parse_design(raw, source_path=str(path))


def _parse_design(raw: Dict[str, Any], source_path: str = "") -> ExperimentDesign:
    """Parse raw JSON into an ExperimentDesign, validating as we go."""
    # -- experiment_id --
    experiment_id = raw.get("experiment_id", "")
    if not experiment_id:
        raise ValueError(
            f"Metadata in '{source_path}' is missing required field 'experiment_id'."
        )

    description = raw.get("description", "")

    # -- attributes --
    raw_attrs = raw.get("attributes", {})
    if not raw_attrs:
        raise ValueError(
            f"Metadata in '{source_path}' has no 'attributes' defined. "
            f"At least one attribute is required."
        )

    attributes: Dict[str, AttributeDefinition] = {}
    for name, attr_def in raw_attrs.items():
        if not isinstance(attr_def, dict):
            raise ValueError(
                f"Attribute '{name}' in '{source_path}' must be an object, "
                f"got {type(attr_def).__name__}."
            )
        attributes[name] = AttributeDefinition(
            name=name,
            attr_type=attr_def.get("type", "binary"),
            description=attr_def.get("description", ""),
        )

    # -- conditions --
    raw_conds = raw.get("conditions", {})
    if not raw_conds:
        raise ValueError(
            f"Metadata in '{source_path}' has no 'conditions' defined. "
            f"At least one condition is required."
        )

    conditions: Dict[str, ConditionDefinition] = {}
    seen_labels: Set[str] = set()

    for cond_name, cond_def in raw_conds.items():
        if not isinstance(cond_def, dict):
            raise ValueError(
                f"Condition '{cond_name}' in '{source_path}' must be an object, "
                f"got {type(cond_def).__name__}."
            )

        label = cond_def.get("label", cond_name)
        if label in seen_labels:
            raise ValueError(
                f"Duplicate label '{label}' in '{source_path}' "
                f"(condition '{cond_name}')."
            )
        seen_labels.add(label)

        cond_attrs = cond_def.get("attributes", {})

        # Validate that all defined attributes are present
        missing = set(attributes.keys()) - set(cond_attrs.keys())
        if missing:
            raise ValueError(
                f"Condition '{cond_name}' in '{source_path}' is missing "
                f"attributes: {sorted(missing)}."
            )

        # Validate attribute values
        for attr_name, attr_value in cond_attrs.items():
            if attr_name not in attributes:
                raise ValueError(
                    f"Condition '{cond_name}' references undefined attribute "
                    f"'{attr_name}' in '{source_path}'."
                )
            attr_def = attributes[attr_name]
            if attr_def.attr_type == "binary" and not isinstance(attr_value, bool):
                raise ValueError(
                    f"Condition '{cond_name}', attribute '{attr_name}': "
                    f"expected boolean for binary attribute, "
                    f"got {type(attr_value).__name__} ({attr_value!r})."
                )
            if attr_def.attr_type == "categorical" and not isinstance(attr_value, str):
                raise ValueError(
                    f"Condition '{cond_name}', attribute '{attr_name}': "
                    f"expected string for categorical attribute, "
                    f"got {type(attr_value).__name__} ({attr_value!r})."
                )

        conditions[cond_name] = ConditionDefinition(
            condition_name=cond_name,
            label=label,
            attributes=dict(cond_attrs),
        )

    return ExperimentDesign(
        experiment_id=experiment_id,
        description=description,
        attributes=attributes,
        conditions=conditions,
    )


# ---------------------------------------------------------------------------
# Validation against observed data
# ---------------------------------------------------------------------------

def validate_metadata_against_data(
    design: ExperimentDesign,
    observed_condition_names: List[str],
) -> List[str]:
    """Check that observed condition names match the metadata.

    Parameters
    ----------
    design : ExperimentDesign
        The loaded experiment design.
    observed_condition_names : list of str
        Condition names found in the data (e.g., from seed_aggregated.csv).

    Returns
    -------
    list of str
        List of warning/error messages.  Empty if all checks pass.
    """
    messages: List[str] = []

    defined_names = set(design.conditions.keys())
    observed_set = set(observed_condition_names)

    # Conditions in data but not in metadata
    undefined = observed_set - defined_names
    if undefined:
        messages.append(
            f"WARNING: {len(undefined)} condition(s) in data but not in metadata: "
            f"{sorted(undefined)}"
        )

    # Conditions in metadata but not in data
    unused = defined_names - observed_set
    if unused:
        messages.append(
            f"WARNING: {len(unused)} condition(s) in metadata but not in data: "
            f"{sorted(unused)}"
        )

    # Duplicate labels
    labels = [c.label for c in design.conditions.values()]
    if len(labels) != len(set(labels)):
        dupes = [l for l in labels if labels.count(l) > 1]
        messages.append(
            f"ERROR: Duplicate labels found: {sorted(set(dupes))}"
        )

    return messages


# ---------------------------------------------------------------------------
# Design inspection
# ---------------------------------------------------------------------------

@dataclass
class DesignInspection:
    """Descriptive summary of an experimental design.

    This is purely descriptive — it does not select statistical models.
    """

    experiment_id: str
    n_conditions: int
    n_attributes: int
    attribute_names: List[str]
    attribute_types: Dict[str, str]  # attr_name -> type
    attribute_levels: Dict[str, Set[Any]]  # attr_name -> set of observed values
    observed_combinations: List[Dict[str, Any]]
    missing_combinations: List[Dict[str, Any]]
    is_complete_factorial: bool
    expected_n_conditions: int
    condition_labels: Dict[str, str]  # condition_name -> label

    def summary(self) -> str:
        """Return a human-readable summary string."""
        lines = [
            f"Experiment: {self.experiment_id}",
            f"Conditions: {self.n_conditions}",
            f"Attributes: {self.n_attributes}",
            "",
            "Attribute details:",
        ]
        for name in self.attribute_names:
            levels = sorted(str(v) for v in self.attribute_levels.get(name, set()))
            lines.append(
                f"  {name} ({self.attribute_types[name]}): "
                f"{len(levels)} level(s) = {levels}"
            )

        lines.append("")
        if self.is_complete_factorial:
            lines.append(
                f"Complete factorial design: {self.n_conditions} conditions "
                f"({self.expected_n_conditions} expected)"
            )
        else:
            lines.append(
                f"Incomplete design: {self.n_conditions} observed conditions "
                f"out of {self.expected_n_conditions} possible combinations"
            )
            if self.missing_combinations:
                lines.append(f"Missing combinations: {len(self.missing_combinations)}")

        lines.append("")
        lines.append("Observed conditions:")
        for cond_name in sorted(self.observed_combinations[0].keys()) if self.observed_combinations else []:
            pass  # filled below

        # Reformat: list each condition with its attributes
        for combo in sorted(self.observed_combinations, key=lambda c: str(c)):
            cond_name = combo.get("_condition_name", "?")
            label = self.condition_labels.get(cond_name, cond_name)
            attr_strs = [
                f"{k}={v}" for k, v in sorted(combo.items())
                if not k.startswith("_")
            ]
            lines.append(f"  {label} ({cond_name}): {', '.join(attr_strs)}")

        return "\n".join(lines)


def inspect_design(design: ExperimentDesign) -> DesignInspection:
    """Produce a descriptive inspection of the experiment design.

    Parameters
    ----------
    design : ExperimentDesign
        The loaded experiment design.

    Returns
    -------
    DesignInspection
        Descriptive summary.
    """
    import itertools

    attribute_names = design.attribute_names
    attribute_types = {
        name: design.attributes[name].attr_type for name in attribute_names
    }

    # Observed levels per attribute
    attribute_levels: Dict[str, Set[Any]] = {}
    for name in attribute_names:
        attribute_levels[name] = design.get_attribute_values(name)

    # Observed combinations
    observed_combinations: List[Dict[str, Any]] = []
    for cond_name in sorted(design.conditions.keys()):
        cond = design.conditions[cond_name]
        combo = {k: cond.attributes[k] for k in attribute_names}
        combo["_condition_name"] = cond_name
        observed_combinations.append(combo)

    # Expected combinations (full factorial)
    level_lists = []
    for name in attribute_names:
        levels = sorted(attribute_levels[name], key=str)
        level_lists.append([(name, v) for v in levels])

    all_combos = list(itertools.product(*level_lists))
    expected_n = len(all_combos)

    # Missing combinations
    observed_set = set()
    for combo in observed_combinations:
        key = tuple(
            (k, combo[k]) for k in attribute_names
        )
        observed_set.add(key)

    missing_combinations: List[Dict[str, Any]] = []
    for combo_tuple in all_combos:
        combo_dict = dict(combo_tuple)
        key = tuple((k, combo_dict[k]) for k in attribute_names)
        if key not in observed_set:
            missing_combinations.append(combo_dict)

    is_complete = len(missing_combinations) == 0

    condition_labels = {
        name: cond.label for name, cond in design.conditions.items()
    }

    return DesignInspection(
        experiment_id=design.experiment_id,
        n_conditions=design.n_conditions,
        n_attributes=design.n_attributes,
        attribute_names=attribute_names,
        attribute_types=attribute_types,
        attribute_levels=attribute_levels,
        observed_combinations=observed_combinations,
        missing_combinations=missing_combinations,
        is_complete_factorial=is_complete,
        expected_n_conditions=expected_n,
        condition_labels=condition_labels,
    )


__all__ = [
    "AttributeDefinition",
    "ConditionDefinition",
    "ExperimentDesign",
    "DesignInspection",
    "load_experiment_design",
    "validate_metadata_against_data",
    "inspect_design",
]
