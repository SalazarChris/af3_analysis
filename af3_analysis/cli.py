"""
Command-line interface for AF3 Confidence Analysis Pipeline.
"""

import argparse
import sys
from pathlib import Path

from .config import AnalysisConfig, load_config, create_run_directory, ConfigValidationError, ConfigPathError
from . import __version__


def create_parser() -> argparse.ArgumentParser:
    """Create the argument parser."""
    parser = argparse.ArgumentParser(
        prog="af3_analysis",
        description="AF3 Confidence Analysis Pipeline",
        epilog=f"Version {__version__}",
    )
    
    parser.add_argument(
        "config",
        type=str,
        help="Path to JSON configuration file",
    )
    
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print version and exit",
    )
    
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate configuration without creating run directory",
    )
    
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing run directory if present",
    )
    
    return parser


def cmd_validate(args: argparse.Namespace) -> int:
    """Validate configuration file."""
    try:
        config = load_config(args.config)
        print(f"Configuration valid: {args.config}")
        print(f"  Run ID: {config.run_id}")
        print(f"  Analysis mode: {config.analysis_mode}")
        print(f"  Output root: {config.output_root}")
        if config.raw_af3_root:
            print(f"  Raw AF3 root: {config.raw_af3_root}")
        if config.legacy_csv_root:
            print(f"  Legacy CSV root: {config.legacy_csv_root}")
        return 0
    except (ConfigValidationError, ConfigPathError) as e:
        print(f"Configuration error: {e}", file=sys.stderr)
        return 1


def cmd_run(args: argparse.Namespace) -> int:
    """Run the analysis pipeline."""
    try:
        config = load_config(args.config)
        
        run_dir = config.output_paths.run_dir
        if run_dir.exists() and not args.force:
            print(f"Run directory already exists: {run_dir}", file=sys.stderr)
            print("Use --force to overwrite.", file=sys.stderr)
            return 1
        
        create_run_directory(config)
        
        # Create run manifest
        manifest_path = config.output_paths.manifest_dir / "run_manifest.json"
        with open(manifest_path, "w") as f:
            import json
            json.dump(config.to_manifest(), f, indent=2)
        
        print(f"Run directory created: {run_dir}")
        print(f"  Manifest: {manifest_path}")
        print(f"  Tables: {config.output_paths.tables_dir}")
        print(f"  Figures: {config.output_paths.figures_dir}")
        print(f"  Reports: {config.output_paths.reports_dir}")
        print(f"  Logs: {config.output_paths.logs_dir}")
        return 0
        
    except (ConfigValidationError, ConfigPathError) as e:
        print(f"Configuration error: {e}", file=sys.stderr)
        return 1


def main(args: list = None) -> int:
    """Main entry point."""
    parser = create_parser()
    parsed = parser.parse_args(args)
    
    if parsed.version:
        print(__version__)
        return 0
    
    if parsed.validate:
        return cmd_validate(parsed)
    
    return cmd_run(parsed)


if __name__ == "__main__":
    sys.exit(main())
