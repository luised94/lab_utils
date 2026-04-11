"""
PyMOL Visualization Session: ORC-Cdc6-DNA Complex (7MCA)
Purpose: Generate publication-quality images of AAA+ fold domains,
         suppressor mutations, and CDC6-ORC1 AAA+ interface.
Author: LEMR
Date: 2025-12-08
"""
from pymol import cmd
import pymol_configuration
import importlib
import sys
import helper_functions

# Force reload all custom modules
for module in [pymol_configuration, helper_functions]:
    importlib.reload(module)

# ============================================================================
# START
# ============================================================================
cmd.reinitialize()

# ============================================================================
# SECTION 1: SESSION-SPECIFIC CONFIGURATION
# ============================================================================

# Structure information
PDB_CODE = "7mca"

# --- Display settings (not covered by pymol_configuration) ---
CARTOON_SIDE_CHAIN_HELPER = 1
STICK_RADIUS = 0.3
SUPPRESSOR_SPHERE_SCALE = 2.0
DISTANCE_DASH_GAP = 0.5
DISTANCE_DASH_RADIUS = 0.1
OVERALL_ROTATION_TRANSPARENCY = 0.4
CDC6_ORC1_CARTOON_TRANSPARENCY = 0.3

# --- Chain colors ---
# NOTE: Full cbc color assignments for chains A-H can be added here once
# confirmed from a live session. Only the explicit non-cbc override is stored.
CHAIN_COLORS: dict[str, str] = {
    "I": "orange",    # Cdc6
}

# --- AAA+ fold domain residue ranges ---
# Keys are chain letters; values are PyMOL residue range strings.
# Iterated in Section 5 to build the 'aaafolds' selection.
AAAFOLDS_CHAIN_RANGES: dict[str, str] = {
    "A": "471-628",
    "B": "300-489",
    "C": "92-273",
    "D": "66-278",
    "E": "31-168",
    "I": "102-263",
}

# --- DNA-in-channel residue ranges ---
DNAINCHANNEL_CHAIN_RANGES: dict[str, str] = {
    "G": "10-31",
    "H": "56-76",
}

# --- Suppressor mutation views ---
# Defined ahead of SUPPRESSOR_RESIDUES so the list can reference them by name.
VIEW_SUPPRESSOR_ORC1_E495 = (
     0.426406741,   -0.417760849,   -0.802281439,
     0.807456851,    0.575538695,    0.129467249,
     0.407653809,   -0.703014195,    0.582741916,
    -0.002135545,    0.000229927,  -87.124885559,
   129.596054077,  138.158264160,  161.141448975,
    66.774276733,  107.287414551,  -20.000000000,
)
VIEW_SUPPRESSOR_ORC3_P481 = (
     0.276172251,   -0.958077013,    0.076264359,
     0.926068425,    0.286498696,    0.245602220,
    -0.257153094,    0.002797240,    0.966364503,
     0.000000000,    0.000000000,  -49.708568573,
   115.162040710,  118.990631104,   76.779479980,
    39.190612793,   60.226524353,  -20.000000000,
)
VIEW_SUPPRESSOR_ORC4_P225 = (
    -0.611087084,    0.467109442,    0.639039159,
     0.651989937,   -0.160774395,    0.740987599,
     0.448863804,    0.869451940,   -0.206302688,
    -0.000000000,    0.000000000,  -84.202651978,
   121.166198730,  117.056610107,  141.704757690,
    66.386016846,  102.019287109,  -20.000000000,
)
VIEW_SUPPRESSOR_ORC5_E104 = (
    -0.752292931,   -0.426453650,    0.502184510,
     0.615914285,   -0.725833297,    0.306290090,
     0.233883113,    0.539721191,    0.808698237,
     0.000000000,    0.000000000,  -72.639938354,
   101.370376587,  118.207054138,  130.884613037,
    57.269882202,   88.009994507,  -20.000000000,
)
VIEW_SUPPRESSOR_ORC6_E305 = (
     0.765645385,    0.437111735,   -0.471919209,
    -0.340094328,    0.897789836,    0.279806972,
     0.545994639,   -0.053733055,    0.836059928,
    -0.001210315,   -0.003005750,  -85.177246094,
   152.192550659,  111.454238892,   63.569774628,
    67.185455322,  103.247833252,  -20.000000000,
)

# --- Suppressor mutation constants ---
# One entry per suppressor. 'label' drives the output filename.
# 'view' is the close-up view for that suppressor's image spec.
SUPPRESSOR_RESIDUES: list[dict] = [
    {
        "chain":          "A",
        "residue_number": "495",
        "label":          "orc1_e495",
        "view":           VIEW_SUPPRESSOR_ORC1_E495,
    },
    {
        "chain":          "C",
        "residue_number": "481",
        "label":          "orc3_p481",
        "view":           VIEW_SUPPRESSOR_ORC3_P481,
    },
    {
        "chain":          "D",
        "residue_number": "225",
        "label":          "orc4_p225",
        "view":           VIEW_SUPPRESSOR_ORC4_P225,
    },
    {
        "chain":          "E",
        "residue_number": "104",
        "label":          "orc5_e104",
        "view":           VIEW_SUPPRESSOR_ORC5_E104,
    },
    {
        "chain":          "F",
        "residue_number": "305",
        "label":          "orc6_e305",
        "view":           VIEW_SUPPRESSOR_ORC6_E305,
    },
]

# Adjacent Orc1 residue shown alongside suppressor sticks in close-up views.
ORC1_ADJACENT_RESIDUE = "485"

# --- "In-the-way" residue ranges ---
# Residues hidden in the suppressor close-up views to clear sightlines.
# Each entry is (chain, residue_range_string).
INTHEWAY_RANGES: list[tuple[str, str]] = [
    ("D", "134-147"),
    ("D", "171-183"),
    ("B", "232-246"),
    ("B", "456-466"),
]

# --- Distance measurement constants ---
DISTANCE_SPECS: list[dict] = [
    {
        "name":  "4PSD",
        "atom1": "/7MCA/D/D/PRO`225/CB",
        "atom2": "/7MCA/H/H/DC`61/P",
    },
    {
        "name":  "5EKK",
        "atom1": "/7MCA/E/E/GLU`104/CD",
        "atom2": "/7MCA/A/A/LYS`362/CA",
    },
    {
        "name":  "5EKD",
        "atom1": "/7MCA/E/E/GLU`104/CD",
        "atom2": "/7MCA/G/G/DT`21/P",
    },
    {
        "name":  "1EKA",
        "atom1": "/7MCA/A/A/GLU`495/CD",
        "atom2": "/lig/K/D/AGS`2001/C5'",
    },
    {
        "name":  "6EKD",
        "atom1": "/7MCA/F/F/GLU`305/CD",
        "atom2": "/7MCA/G/G/DA`46/P",
    },
]

# --- CDC6-ORC1 AAA+ interface constants ---
C6O1: dict = {
    "selection_name":  "c6o1",
    "orc1_chain":      "A",
    "orc1_residue":    "616",
    "cdc6_chain":      "I",
    "cdc6_residues":   ["114", "224"],
    "transparency":    CDC6_ORC1_CARTOON_TRANSPARENCY,
}

# --- View matrices ---
VIEW_AAAFOLDS_OVERVIEW = (
    -0.527550519,   -0.201993212,    0.825157940,
    -0.841823280,    0.254780442,   -0.475835592,
    -0.114118114,   -0.945666373,   -0.304451793,
    -0.000000000,    0.000000000, -323.199737549,
   107.580001831,  120.430999756,  127.585998535,
   254.813293457,  391.586181641,  -20.000000000,
)
VIEW_SUPPRESSORS_OVERALL = (
    -0.494458616,   -0.225710392,    0.839386284,
    -0.869120657,    0.141693920,   -0.473868966,
    -0.011976970,   -0.963827670,   -0.266235441,
     0.000000000,    0.000000000, -509.880279541,
   117.735023499,  114.526748657,  114.030647278,
   401.993499756,  617.767089844,  -20.000000000,
)
VIEW_CDC6_ORC1_INTERFACE = (
     0.329257607,   -0.103709482,   -0.938521206,
     0.927203953,    0.223447070,    0.300590754,
     0.178538367,   -0.969175398,    0.169733956,
     0.000000000,    0.000000000,  -53.428050995,
   127.544395447,  152.950195312,  131.362655640,
    42.123081207,   64.733016968,  -20.000000000,
)

# ============================================================================
# SECTION 2: INITIALIZE SESSION
# ============================================================================
print("=" * 60)
print(f"PyMOL Visualization Session: {PDB_CODE}")
print("=" * 60)

# Apply shared settings
cmd.bg_color(pymol_configuration.BACKGROUND_COLOR)
cmd.set("hash_max", pymol_configuration.HASH_MAX)
cmd.set("antialias", pymol_configuration.ANTIALIAS)
cmd.set("ray_shadows", pymol_configuration.RAY_SHADOWS)
cmd.set("ray_shadow", pymol_configuration.RAY_SHADOW)
cmd.set("ray_shadow_decay_factor", pymol_configuration.RAY_SHADOW_DECAY_FACTOR)
cmd.set("ray_shadow_decay_range", pymol_configuration.RAY_SHADOW_DECAY_RANGE)
cmd.set("ray_trace_fog", pymol_configuration.RAY_TRACE_FOG)
cmd.set("depth_cue", pymol_configuration.DEPTH_CUE)
cmd.set("fog", pymol_configuration.FOG)

# Apply script-local display settings
cmd.set("cartoon_side_chain_helper", CARTOON_SIDE_CHAIN_HELPER)

# Load structure
print(f"\nLoading structure {PDB_CODE}...")
cmd.fetch(PDB_CODE, type="pdb")

# ============================================================================
# SECTION 3: STRUCTURE PREPARATION
# ============================================================================
helper_functions.prepare_structure()

# ============================================================================
# SECTION 4: CUSTOM CHAIN COLORING
# ============================================================================
print("\nColoring chains...")
cmd.util.cbc()
print("  Applied cbc base colors to all chains")
for chain, color in CHAIN_COLORS.items():
    selection = f"chain {chain}"
    if cmd.count_atoms(selection) > 0:
        cmd.color(color, selection)
        print(f"  Chain {chain}: {color}")
