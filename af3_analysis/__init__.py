"""
AF3 Confidence Analysis Pipeline
================================

A reproducible, registry-driven scientific analysis package for
AlphaFold 3 confidence metric evaluation.

Version: 0.1.0
"""

from importlib.metadata import version, PackageNotFoundError

__version__ = "0.1.0"

try:
    __version__ = version("af3_analysis")
except PackageNotFoundError:
    pass

from af3_analysis.config import create_config_interactive

__all__ = [
    "config",
    "errors",
    "logging_utils",
    "workflow",
    "schemas",
    "registry",
    "io",
    "preprocessing",
    "statistics",
    "structural",
    "visualization",
    "reporting",
    "pipeline",
    "create_config_interactive",
]
