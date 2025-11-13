#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas", "crispr-cas9-designer"]
# ///

"""Quick dependency check for CRISPR scoring script."""

try:
    import pandas as pd
    print(f" pandas {pd.__version__}")
    
    import crispr_cas9_designer
    print(f" crispr-cas9-designer loaded")
    
    print("\nAll dependencies available and importable")
    
except ImportError as e:
    print(f"Dependency issue: {e}")
    exit(1)
