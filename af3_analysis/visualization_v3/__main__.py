"""
Command-line interface for the V3 visualization pipeline.

Usage:
    python -m af3_analysis.visualization_v3 <run_dir> [options]

The run directory must be an existing pipeline run directory containing
tables/. V3 writes exclusively to <run_dir>/v3/ and never modifies
existing pipeline outputs.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from .config import V3Config, load_v3_config
from .runner import run_v3_pipeline


def create_parser() -> argparse.ArgumentParser:
    """Create the V3 argument parser."""
    parser = argparse.ArgumentParser(
        prog="af3_analysis.visualization_v3",
        description=(
            "V3 structural visualization pipeline. Consumes an existing "
            "pipeline run directory (and optionally raw AF3 outputs and "
            "experiment metadata) and writes figures/tables/report to "
            "<run_dir>/v3/."
        ),
    )

    parser.add_argument(
        "run_dir",
        type=str,
        help="Existing pipeline run directory (contains tables/)",
    )

    parser.add_argument(
        "--raw-af3-root",
        type=str,
        default=None,
        help="Raw AF3 output root containing CIF files (required for "
             "structural figures)",
    )

    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        help="Path to experiment_metadata.json (enables labels and "
             "factorial contrasts)",
    )

    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to V3 JSON configuration file",
    )

    parser.add_argument(
        "--figures",
        type=str,
        default=None,
        help="Comma-separated figure IDs to generate (e.g. F01,F02,F12). "
             "Default: all enabled figures from config",
    )

    parser.add_argument(
        "--reference",
        type=str,
        default=None,
        help="Reference condition id (overrides config)",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if <run_dir>/v3 already contains output",
    )

    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress per-stage console logging",
    )

    return parser


def main(args: list = None) -> int:
    """CLI entry point."""
    parser = create_parser()
    parsed = parser.parse_args(args)

    run_dir = Path(parsed.run_dir).resolve()
    if not run_dir.is_dir():
        print(f"Run directory not found: {run_dir}", file=sys.stderr)
        return 1

    tables_dir = run_dir / "tables"
    if not tables_dir.is_dir():
        print(f"Run directory has no tables/ directory: {run_dir}",
              file=sys.stderr)
        print("V3 consumes an existing pipeline run directory.",
              file=sys.stderr)
        return 1

    # Load configuration
    if parsed.config:
        v3_config = load_v3_config(parsed.config)
    else:
        v3_config = V3Config()

    # CLI overrides
    if parsed.reference:
        # reference is a frozen dataclass field holding a dict; rebuild it
        v3_config = _with_reference(v3_config, parsed.reference)
    if parsed.figures:
        ids = [f.strip().upper() for f in parsed.figures.split(",") if f.strip()]
        v3_config = _with_figures(v3_config, ids)

    # Overwrite protection
    v3_dir = run_dir / "v3"
    if v3_dir.exists() and any(v3_dir.iterdir()) and not parsed.force:
        print(f"V3 output already exists: {v3_dir}", file=sys.stderr)
        print("Use --force to re-run.", file=sys.stderr)
        return 1

    # Logging to console (pipeline logs also go to the run's log file setup)
    logging.basicConfig(
        level=logging.WARNING if parsed.quiet else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    raw_af3_root = Path(parsed.raw_af3_root).resolve() if parsed.raw_af3_root else None
    if parsed.raw_af3_root and not raw_af3_root.is_dir():
        print(f"Raw AF3 root not found: {raw_af3_root}", file=sys.stderr)
        return 1

    metadata_path = Path(parsed.metadata).resolve() if parsed.metadata else None
    if parsed.metadata and not metadata_path.is_file():
        print(f"Metadata file not found: {metadata_path}", file=sys.stderr)
        return 1

    results = run_v3_pipeline(
        run_dir,
        raw_af3_root=raw_af3_root,
        experiment_metadata_path=metadata_path,
        v3_config=v3_config,
    )

    # Print summary
    summary = results.get("summary", {})
    print()
    print("V3 pipeline finished:")
    print(f"  Status:  {results.get('status', 'unknown')}")
    print(f"  Figures: {summary.get('n_figures_success', 0)} ok, "
          f"{summary.get('n_figures_skipped', 0)} skipped, "
          f"{summary.get('n_figures_failed', 0)} failed")
    print(f"  Tables:  {len(results.get('tables', {}))}")
    print(f"  Output:  {v3_dir}")

    if results.get("errors"):
        print("  Errors:")
        for err in results["errors"][:10]:
            print(f"    - {err}")

    return 0 if results.get("status") in ("complete", "completed_with_errors") else 1


def _with_reference(v3_config: V3Config, reference: str) -> V3Config:
    """Return a copy of the config with the reference condition set."""
    from dataclasses import replace
    return replace(v3_config, reference={"condition": reference})


def _with_figures(v3_config: V3Config, figure_ids: list) -> V3Config:
    """Return a copy of the config with only the given figures enabled."""
    from dataclasses import replace
    figures = {fig_id: (fig_id in figure_ids) for fig_id in list(v3_config.figures) or []}
    # Ensure unlisted figures are disabled when an explicit list is given
    all_ids = [
        "F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F09", "F10",
        "F11", "F12", "F13", "F14", "F15", "F16", "F17", "F18", "F19", "F20",
    ]
    for fig_id in all_ids:
        figures[fig_id] = fig_id in figure_ids
    return replace(v3_config, figures=figures)


if __name__ == "__main__":
    sys.exit(main())
