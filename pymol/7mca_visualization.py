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

# --- Selection names ---
# Defined here so Section 6 lambdas reference names, not string literals.
AAAFOLDS_SELECTION     = "aaafolds"
DNAINCHANNEL_SELECTION = "dnainchannel"
SUPPRESSORS_SELECTION  = "suppressors"
NEAR_SELECTION         = "near"
INTHEWAY_SELECTION     = "intheway"

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

# ORC1 chain and adjacent residue shown alongside suppressor sticks in close-up views.
ORC1_CHAIN            = "A"
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

# ============================================================================
# SECTION 5: CREATE SELECTIONS AND DISTANCE OBJECTS
# ============================================================================
print("\n" + "=" * 60)
print("Creating selections")
print("=" * 60)

# --- AAA+ fold domains ---
print("\nBuilding 'aaafolds' selection...")
aaafolds_parts: list[str] = []
for chain, residue_range in AAAFOLDS_CHAIN_RANGES.items():
    aaafolds_parts.append(f"(chain {chain} and resi {residue_range})")
aaafolds_selection_string = " or ".join(aaafolds_parts)
cmd.select(AAAFOLDS_SELECTION, aaafolds_selection_string)
aaafolds_atom_count = cmd.count_atoms(AAAFOLDS_SELECTION)
if aaafolds_atom_count == 0:
    raise ValueError(f"ERROR: Selection '{AAAFOLDS_SELECTION}' is empty! Check AAAFOLDS_CHAIN_RANGES.")
print(f"  {AAAFOLDS_SELECTION}: {aaafolds_atom_count} atoms ({len(AAAFOLDS_CHAIN_RANGES)} chains)")

# --- DNA in channel ---
print("\nBuilding 'dnainchannel' selection...")
dnainchannel_parts: list[str] = []
for chain, residue_range in DNAINCHANNEL_CHAIN_RANGES.items():
    dnainchannel_parts.append(f"(chain {chain} and resi {residue_range})")
dnainchannel_selection_string = " or ".join(dnainchannel_parts)
cmd.select(DNAINCHANNEL_SELECTION, dnainchannel_selection_string)
dnainchannel_atom_count = cmd.count_atoms(DNAINCHANNEL_SELECTION)
if dnainchannel_atom_count == 0:
    raise ValueError(f"ERROR: Selection '{DNAINCHANNEL_SELECTION}' is empty! Check DNAINCHANNEL_CHAIN_RANGES.")
print(f"  {DNAINCHANNEL_SELECTION}: {dnainchannel_atom_count} atoms")

# --- Suppressor mutations ---
print("\nBuilding 'suppressors' selection...")
suppressors_parts: list[str] = []
for suppressor in SUPPRESSOR_RESIDUES:
    suppressors_parts.append(f"(chain {suppressor['chain']} and resi {suppressor['residue_number']})")
suppressors_selection_string = " or ".join(suppressors_parts)
cmd.select(SUPPRESSORS_SELECTION, suppressors_selection_string)
suppressors_atom_count = cmd.count_atoms(SUPPRESSORS_SELECTION)
if suppressors_atom_count == 0:
    raise ValueError(f"ERROR: Selection '{SUPPRESSORS_SELECTION}' is empty! Check SUPPRESSOR_RESIDUES.")
print(f"  {SUPPRESSORS_SELECTION}: {suppressors_atom_count} atoms ({len(SUPPRESSOR_RESIDUES)} residues)")

# --- Residues near suppressors ---
print("\nBuilding 'near' selection...")
cmd.select(NEAR_SELECTION, f"{SUPPRESSORS_SELECTION} expand 8")
near_atom_count = cmd.count_atoms(NEAR_SELECTION)
print(f"  {NEAR_SELECTION}: {near_atom_count} atoms")

# --- In-the-way residues (hidden in suppressor close-up views) ---
print("\nBuilding 'intheway' selection...")
intheway_parts: list[str] = []
for chain, residue_range in INTHEWAY_RANGES:
    intheway_parts.append(f"(chain {chain} and resi {residue_range})")
intheway_selection_string = " or ".join(intheway_parts)
cmd.select(INTHEWAY_SELECTION, intheway_selection_string)
intheway_atom_count = cmd.count_atoms(INTHEWAY_SELECTION)
print(f"  {INTHEWAY_SELECTION}: {intheway_atom_count} atoms")

# --- CDC6-ORC1 AAA+ interface ---
print("\nBuilding 'c6o1' selection...")
cdc6_residue_parts: list[str] = [f"resi {residue}" for residue in C6O1["cdc6_residues"]]
cdc6_residue_expression = " or ".join(cdc6_residue_parts)
c6o1_selection_string = (
    f"(chain {C6O1['orc1_chain']} and resi {C6O1['orc1_residue']}) or "
    f"(chain {C6O1['cdc6_chain']} and ({cdc6_residue_expression})) or "
    f"(lig and chain {C6O1['cdc6_chain']})"
)
cmd.select(C6O1["selection_name"], c6o1_selection_string)
c6o1_atom_count = cmd.count_atoms(C6O1["selection_name"])
if c6o1_atom_count == 0:
    raise ValueError(f"ERROR: Selection '{C6O1['selection_name']}' is empty! Check C6O1 constants.")
print(f"  {C6O1['selection_name']}: {c6o1_atom_count} atoms")

print(f"\nû Successfully created all selections")

# --- Distance objects ---
print("\n" + "=" * 60)
print("Creating distance objects")
print("=" * 60)
for distance_spec in DISTANCE_SPECS:
    cmd.distance(distance_spec["name"], distance_spec["atom1"], distance_spec["atom2"])
    print(f"  {distance_spec['name']}: {distance_spec['atom1']}  {distance_spec['atom2']}")

# Apply distance display settings and suppress labels
cmd.set("dash_gap", DISTANCE_DASH_GAP)
cmd.set("dash_radius", DISTANCE_DASH_RADIUS)
cmd.hide("labels")
print(f"\n  dash_gap={DISTANCE_DASH_GAP}, dash_radius={DISTANCE_DASH_RADIUS}, labels hidden")
print(f"\nû Successfully created {len(DISTANCE_SPECS)} distance objects")

# ============================================================================
# SECTION 6: DEFINE IMAGE SPECIFICATIONS
# ============================================================================
print("\n" + "=" * 60)
print("Defining image specifications")
print("=" * 60)

image_specs: list[dict] = []

# --- AAA+ fold domain overview ---
# reset_scene_to_base_cartoon shows cartoon for all polymer; the subsequent
# hide("everything") narrows display to only aaafolds and dnainchannel.
image_specs.append({
    "filename":    "7mca_aaafolds_overview.png",
    "description": "AAA+ fold domains with DNA in channel",
    "view":        VIEW_AAAFOLDS_OVERVIEW,
    "actions": [
        lambda: cmd.util.cbc(),
        lambda: helper_functions.reset_scene_to_base_cartoon(CHAIN_COLORS),
        lambda: cmd.hide("everything"),
        lambda: cmd.show("cartoon", AAAFOLDS_SELECTION),
        lambda: cmd.show("cartoon", DNAINCHANNEL_SELECTION),
        lambda: cmd.show("sticks", "lig"),
    ],
})

# --- Per-suppressor close-up views (generated from SUPPRESSOR_RESIDUES) ---
# Action list is identical for all five suppressors; view and filename vary per entry.
for suppressor in SUPPRESSOR_RESIDUES:
    image_specs.append({
        "filename":    f"7mca_suppressors_{suppressor['label']}.png",
        "description": f"Suppressor close-up: {suppressor['label']}",
        "view":        suppressor["view"],
        "actions": [
            lambda: cmd.util.cbc(),
            lambda: helper_functions.reset_scene_to_base_cartoon(CHAIN_COLORS),
            lambda: cmd.show("dashes"),
            lambda: cmd.set("stick_radius", STICK_RADIUS),
            lambda: cmd.show("sticks", f"{SUPPRESSORS_SELECTION} or (chain {ORC1_CHAIN} and resi {ORC1_ADJACENT_RESIDUE})"),
            lambda: cmd.hide("everything", INTHEWAY_SELECTION),
        ],
    })

# --- Suppressor overall rotation - plain ---
image_specs.append({
    "filename":        "7mca_suppressors_overall.png",
    "description":     "Full ORC-Cdc6 structure with suppressor sites as red spheres (rotation series)",
    "view":            VIEW_SUPPRESSORS_OVERALL,
    "rotation_angles": [0, 90, 180, 270],
    "actions": [
        lambda: cmd.util.cbc(),
        lambda: helper_functions.reset_scene_to_base_cartoon(CHAIN_COLORS),
        lambda: cmd.set("sphere_scale", SUPPRESSOR_SPHERE_SCALE),
        lambda: cmd.show("spheres", SUPPRESSORS_SELECTION),
        lambda: cmd.color("red", SUPPRESSORS_SELECTION),
    ],
})

# --- Suppressor overall rotation - 40% cartoon transparency ---
# Filename format: 7mca_suppressors_overall_trans40_{angle:03d}.png
# The render loop appends the angle to the base filename (everything before .png),
# yielding e.g. 7mca_suppressors_overall_trans40_000.png. This differs from the
# decisions document spec (angle before _trans40); Commit 4 can add an optional
# filename suffix key to the spec if the original ordering is required.
image_specs.append({
    "filename":        "7mca_suppressors_overall_trans40.png",
    "description":     "Full ORC-Cdc6 structure with suppressor sites as red spheres, transparent cartoon (rotation series)",
    "view":            VIEW_SUPPRESSORS_OVERALL,
    "rotation_angles": [0, 90, 180, 270],
    "actions": [
        lambda: cmd.util.cbc(),
        lambda: helper_functions.reset_scene_to_base_cartoon(CHAIN_COLORS),
        lambda: cmd.set("sphere_scale", SUPPRESSOR_SPHERE_SCALE),
        lambda: cmd.show("spheres", SUPPRESSORS_SELECTION),
        lambda: cmd.color("red", SUPPRESSORS_SELECTION),
        lambda: cmd.set("cartoon_transparency", OVERALL_ROTATION_TRANSPARENCY),
    ],
})

# --- CDC6-ORC1 AAA+ interface ---
# cartoon_transparency is reset to 0 globally before the selective 0.3 is applied,
# to clear any residual transparency from the preceding overall_trans40 spec.
image_specs.append({
    "filename":    "7mca_cdc6orc1_interface.png",
    "description": "CDC6-ORC1 AAA+ interface residues as sticks with transparent cartoon context",
    "view":        VIEW_CDC6_ORC1_INTERFACE,
    "actions": [
        lambda: cmd.util.cbc(),
        lambda: helper_functions.reset_scene_to_base_cartoon(CHAIN_COLORS),
        lambda: cmd.set("cartoon_transparency", 0),
        lambda: cmd.show("sticks", C6O1["selection_name"]),
        lambda: cmd.set("cartoon_transparency", C6O1["transparency"], PDB_CODE),
    ],
})

print(f"  Defined {len(image_specs)} image specifications")

# ============================================================================
# SECTION 7: RENDER ALL IMAGES
# ============================================================================
print("\n" + "=" * 60)
print("Rendering images")
print("=" * 60)

for i, spec in enumerate(image_specs, 1):
    print(f"\nImage spec {i}/{len(image_specs)}: {spec['filename']}")
    print(f"  Description: {spec['description']}")

    # Execute all scene-setup actions for this spec
    for action in spec['actions']:
        action()

    # Set the view for this spec
    cmd.set_view(spec['view'])

    # Rotation series or single image
    if 'rotation_angles' in spec:
        rotation_angles = spec['rotation_angles']
        base_filename = spec['filename'].replace('.png', '')
        print(f"  Generating {len(rotation_angles)} rotation views...")
        for angle_index, angle in enumerate(rotation_angles):
            angle_filename = f"{base_filename}_{angle:03d}.png"
            print(f"    Angle {angle}ø...")
            helper_functions.save_png(
                filename=angle_filename,
                width=pymol_configuration.IMAGE_WIDTH,
                height=pymol_configuration.IMAGE_HEIGHT,
                dpi=pymol_configuration.IMAGE_DPI,
                ray=pymol_configuration.RAY_TRACE,
                overwrite=pymol_configuration.FLAG_OVERWRITE,
                output_dir=pymol_configuration.OUTPUT_DIR,
            )
            is_last_angle = angle_index == len(rotation_angles) - 1
            if not is_last_angle:
                next_angle = rotation_angles[angle_index + 1]
                rotation_delta = next_angle - angle
                cmd.turn("y", rotation_delta)
    else:
        helper_functions.save_png(
            filename=spec['filename'],
            width=pymol_configuration.IMAGE_WIDTH,
            height=pymol_configuration.IMAGE_HEIGHT,
            dpi=pymol_configuration.IMAGE_DPI,
            ray=pymol_configuration.RAY_TRACE,
            overwrite=pymol_configuration.FLAG_OVERWRITE,
            output_dir=pymol_configuration.OUTPUT_DIR,
        )

# ============================================================================
# SECTION 8: SAVE SESSION AND FINAL SUMMARY
# ============================================================================
print("\n" + "=" * 60)
print("Saving session files")
print("=" * 60)

import os
session_file = os.path.join(pymol_configuration.OUTPUT_DIR, f"{PDB_CODE}_session.pse")
try:
    cmd.save(session_file)
    print(f"  Saved PyMOL session: {session_file}")
except Exception as e:
    print(f"  ERROR saving session: {e}")

print("\n" + "=" * 60)
print("SESSION COMPLETE")
print("=" * 60)
print(f"Structure: {PDB_CODE}")
print(f"Generated {len(image_specs)} image specs in: {pymol_configuration.OUTPUT_DIR}")
print("=" * 60)
