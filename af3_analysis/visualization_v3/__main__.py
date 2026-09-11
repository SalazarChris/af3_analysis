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

from .config import V3Config
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
        "--figures",
        type=str,
        default=None,
        help="Comma-separated figure IDs to generate (e.g. F01,F02,F12). "
             "Default: all figures enabled",
    )

    parser.add_argument(
        "--reference",
        type=str,
        default=None,
        help="Reference condition id (all structural differences are "
             "computed relative to this condition)",
    )

    parser.add_argument(
        "--site",
        action="append",
        default=None,
        metavar="LABEL:CHAIN:RESIDUE[:RADIUS]",
        help="Local-geometry site for figure F10, repeatable. "
             "Example: --site my_site:A:101:8.0",
    )

    parser.add_argument(
        "--region",
        action="append",
        default=None,
        metavar="LABEL:CHAIN:START:END",
        help="Domain/region definition for figure F11, repeatable. "
             "Example: --region my_domain:A:1:65",
    )

    parser.add_argument(
        "--auto-sites",
        action="store_true",
        help="Detect F10 sites automatically as concentrations of "
             "between-condition displacement (used only when --site is "
             "not given; output is labeled as predicted structural "
             "sites)",
    )

    parser.add_argument(
        "--auto-regions",
        action="store_true",
        help="Detect F11 regions automatically as connected components "
             "of the reference contact graph (used only when --region is "
             "not given; output is labeled as predicted structural "
             "regions)",
    )

    parser.add_argument(
        "--contact-distance",
        type=float,
        default=None,
        metavar="ANGSTROM",
        help="Contact-map distance cutoff in Angstrom (default: 8.0)",
    )

    parser.add_argument(
        "--min-coverage",
        type=float,
        default=None,
        metavar="FRACTION",
        help="Minimum coordinate coverage for valid comparisons "
             "(default: 0.80)",
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

    # Base configuration: defaults only, refined by CLI arguments below.
    # No external config file is used.
    v3_config = V3Config()

    # CLI arguments
    if parsed.reference:
        # reference is a frozen dataclass field holding a dict; rebuild it
        v3_config = _with_reference(v3_config, parsed.reference)
    if parsed.figures:
        ids = [f.strip().upper() for f in parsed.figures.split(",") if f.strip()]
        v3_config = _with_figures(v3_config, ids)

    if parsed.site:
        v3_config = _with_sites(v3_config, parsed.site)

    if parsed.region:
        v3_config = _with_regions(v3_config, parsed.region)

    if parsed.contact_distance is not None or parsed.min_coverage is not None:
        from dataclasses import replace
        structure = v3_config.structure
        if parsed.contact_distance is not None:
            structure = replace(structure, contact_distance=parsed.contact_distance)
        if parsed.min_coverage is not None:
            structure = replace(structure, minimum_coverage=parsed.min_coverage)
        v3_config = replace(v3_config, structure=structure)

    if parsed.auto_sites or parsed.auto_regions:
        v3_config = replace(
            v3_config,
            auto_sites=parsed.auto_sites,
            auto_regions=parsed.auto_regions,
        )

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


def _with_sites(v3_config: V3Config, site_specs: list) -> V3Config:
    """Return a copy of the config with parsed --site definitions.

    Each spec has the form LABEL:CHAIN:RESIDUE[:RADIUS]. RADIUS is optional
    and defaults to the pipeline default (8.0 Angstrom).
    """
    sites = []
    for spec in site_specs:
        parts = [p.strip() for p in spec.split(":")]
        if len(parts) not in (3, 4):
            raise SystemExit(
                f"Invalid --site spec '{spec}'. "
                "Expected LABEL:CHAIN:RESIDUE[:RADIUS] "
                "e.g. --site my_site:A:101:8.0"
            )
        label, chain, residue = parts[0], parts[1], parts[2]
        try:
            residue = int(residue)
        except ValueError:
            raise SystemExit(
                f"Invalid --site spec '{spec}': RESIDUE must be an integer"
            )
        site = {"label": label, "chain": chain, "residue": residue}
        if len(parts) == 4:
            try:
                site["radius"] = float(parts[3])
            except ValueError:
                raise SystemExit(
                    f"Invalid --site spec '{spec}': RADIUS must be a number"
                )
        sites.append(site)
    from dataclasses import replace
    return replace(v3_config, sites=sites)


def _with_regions(v3_config: V3Config, region_specs: list) -> V3Config:
    """Return a copy of the config with parsed --region definitions.

    Each spec has the form LABEL:CHAIN:START:END.
    """
    regions = []
    for spec in region_specs:
        parts = [p.strip() for p in spec.split(":")]
        if len(parts) != 4:
            raise SystemExit(
                f"Invalid --region spec '{spec}'. "
                "Expected LABEL:CHAIN:START:END "
                "e.g. --region my_domain:A:1:65"
            )
        label, chain, start, end = parts
        try:
            start, end = int(start), int(end)
        except ValueError:
            raise SystemExit(
                f"Invalid --region spec '{spec}': START and END must be integers"
            )
        if start > end:
            raise SystemExit(
                f"Invalid --region spec '{spec}': START must be <= END"
            )
        regions.append({"label": label, "chain": chain, "start": start, "end": end})
    from dataclasses import replace
    return replace(v3_config, regions=regions)


if __name__ == "__main__":
    sys.exit(main())
