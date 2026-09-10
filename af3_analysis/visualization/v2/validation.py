"""
V2 Validation & Manifest
========================
Post-generation validation checks and machine-readable manifest.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


def validate_figure(
    fig_name: str,
    output_path: Path,
    *,
    n_observations: int = 0,
    n_legend_entries: int = 0,
    max_legend: int = 12,
    expected_metrics: List[str] = None,
    plotted_metrics: List[str] = None,
    warnings: List[str] = None,
) -> Dict[str, Any]:
    """Run standard post-generation validation on a V2 figure.

    Returns
    -------
    dict with ``passed`` (bool) and ``messages`` (list of str).
    """
    messages: List[str] = []
    passed = True

    # 1. Output file exists and is non-empty
    if not output_path.exists():
        messages.append(f"FAIL: {fig_name} output file not found: {output_path}")
        passed = False
    elif output_path.stat().st_size == 0:
        messages.append(f"FAIL: {fig_name} output file is empty")
        passed = False

    # 2. Zero observations
    if n_observations == 0:
        messages.append(f"WARNING: {fig_name} has zero observations plotted")
        passed = False

    # 3. Legend size
    if n_legend_entries > max_legend:
        messages.append(
            f"WARNING: {fig_name} has {n_legend_entries} legend entries "
            f"(limit {max_legend})"
        )

    # 4. Missing expected metrics
    if expected_metrics and plotted_metrics:
        missing = set(expected_metrics) - set(plotted_metrics)
        if missing:
            messages.append(
                f"WARNING: {fig_name} missing expected metrics: {sorted(missing)}"
            )

    # 5. Image dimensions (reasonable range)
    if output_path.exists() and output_path.suffix == ".png":
        try:
            from PIL import Image
            with Image.open(output_path) as img:
                w, h = img.size
                if h > 16000:
                    messages.append(
                        f"WARNING: {fig_name} height is {h}px — may be too tall"
                    )
                if w > 16000:
                    messages.append(
                        f"WARNING: {fig_name} width is {w}px — may be too wide"
                    )
        except ImportError:
            pass  # PIL not available — skip dimension check
        except Exception as e:
            messages.append(f"NOTE: Could not read image dimensions: {e}")

    # 6. Carry through pre-existing warnings
    if warnings:
        for w in warnings:
            messages.append(w)

    return {"passed": passed, "messages": messages}


def write_manifest(
    save_path: Path,
    *,
    version: str = "v2",
    figures: List[Dict[str, Any]],
    environment_strategy: str = "pooled",
    confounded_by_tokenisation: bool = True,
    unit_of_analysis: str = "seed_level",
) -> Path:
    """Write a machine-readable V2 visualization manifest.

    Parameters
    ----------
    save_path : Path
        Directory to write ``visualization_manifest.json``.
    figures : list of dict
        Per-figure metadata.
    environment_strategy : str
        How environments are handled: ``"pooled"``, ``"faceted"``,
        ``"filtered"``.
    confounded_by_tokenisation : bool
        Whether mean pLDDT is confounded by tokenisation.
    unit_of_analysis : str
        Description of the analysis unit.

    Returns
    -------
    Path to the written manifest file.
    """
    manifest = {
        "version": version,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "environment_strategy": environment_strategy,
        "confounded_by_tokenisation": confounded_by_tokenisation,
        "unit_of_analysis": unit_of_analysis,
        "figures": figures,
    }

    manifest_path = save_path / "visualization_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)

    return manifest_path
