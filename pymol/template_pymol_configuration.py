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
FLAG_OVERWRITE = True

# ============================================================================
# RENDERING QUALITY
# ============================================================================
BACKGROUND_COLOR = "white"
HASH_MAX = 200
ANTIALIAS = 2
RAY_SHADOWS = 1

# Ray shadow fine-tuning
RAY_SHADOW = 0.2                    # Shadow darkness (0.5=default, lower=lighter)
RAY_SHADOW_DECAY_FACTOR = 0.1       # How fast shadow fades
RAY_SHADOW_DECAY_RANGE = 2          # Shadow fade range

# Fog and depth settings
RAY_TRACE_FOG = 0                   # Fog effect (0=off)
DEPTH_CUE = 1                       # Depth cueing (0=off)
FOG = 0                             # General fog (0=off)

# ============================================================================
# STANDARD OUTPUT
# ============================================================================
OUTPUT_DIR = "pymol_outputs"
