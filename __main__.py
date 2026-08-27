"""
Entry point for running as a module: python -m af3_analysis
"""

from af3_analysis.cli import main
import sys

if __name__ == "__main__":
    sys.exit(main())
