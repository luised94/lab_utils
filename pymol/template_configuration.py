"""
PyMOL Shared Configuration
Settings that apply across all visualization sessions.
Modify these values to change behavior for all scripts.
"""

# ============================================================================
# IMAGE EXPORT SETTINGS
# ============================================================================
IMAGE_WIDTH = 1200
IMAGE_HEIGHT = 900
IMAGE_DPI = 300
RAY_TRACE = True
FLAG_OVERWRITE = True  # Set False to prevent overwriting existing images

# ============================================================================
# RENDERING QUALITY
# ============================================================================
BACKGROUND_COLOR = "white"
HASH_MAX = 200          # Higher = better ray-traced quality (50-250)
ANTIALIAS = 2           # Anti-aliasing level (0-4)
RAY_SHADOWS = 0         # 0=off, 1=on

# ============================================================================
# STANDARD OUTPUT
# ============================================================================
OUTPUT_DIR = "pymol_outputs"
