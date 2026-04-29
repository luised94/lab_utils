"""
PyMOL Visualization Session: TFIIA-TBP-DNA Complex (1NH2)
Purpose: Generate publication-quality images showing hydrophobic core
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
PDB_CODE = "1nh2"

# Chain colors (customize per structure)
CHAIN_COLORS = {
    "A": "yellow",    # TBP
    "B": "marine",    # TOA1
    "C": "marine",    # TOA1 second chain
    "D": "green",     # TOA2
}

# Selection parameters
DISTANCE_CUTOFF = 4.0
HYDROPHOBIC_RESIDUES = "ala+gly+val+ile+leu+phe+met+pro+trp"

# ROI constants - one dict per residue of interest
# All PyMOL selection strings for a given ROI are derived from these values
G16 = {
    "chain": "D",
    "residue_number": "16",
    "roi_selection": "roi_g16",
    "nearby_selection": "nearby_g16",
    "core_selection": "hydrophobic_core_g16",
    "contact_chain": "B",
    "nearby_contact_selection": "nearby_chain_b_g16",
}

V251 = {
    "chain": "C",
    "residue_number": "251",
    "roi_selection": "roi_v251",
    "nearby_selection": "nearby_v251",
    "core_selection": "hydrophobic_core_v251",
    "contact_chain": "C",
    "nearby_contact_selection": "nearby_chain_c_v251",
}

# ============================================================================
# ROI (RESIDUE OF INTEREST) SPECIFICATIONS
# ============================================================================
# Define multiple ROIs with explicit selection actions
# Each ROI creates: roi_NAME, nearby_NAME, hydrophobic_core_NAME selections
ROI_SPECS = [
    {
        "name": "g16",
        "description": "G16 hydrophobic core in TOA2",
        "roi_selection": G16["roi_selection"],
        "core_selection": G16["core_selection"],
        "actions": [
            lambda: cmd.select(
                G16["roi_selection"],
                f"chain {G16['chain']} and resi {G16['residue_number']}",
            ),
            lambda: cmd.select(
                G16["nearby_selection"],
                f"byres (all within {DISTANCE_CUTOFF} of {G16['roi_selection']})",
            ),
            lambda: cmd.select(
                G16["core_selection"],
                f"{G16['nearby_selection']} and resn {HYDROPHOBIC_RESIDUES} and not {G16['roi_selection']}",
            ),
            lambda: cmd.select(
                G16["nearby_contact_selection"],
                f"byres (chain {G16['contact_chain']} within {DISTANCE_CUTOFF} of {G16['core_selection']})",
            ),
            lambda: cmd.select(
                G16["core_selection"],
                f"{G16['core_selection']} or ({G16['nearby_contact_selection']} and resn {HYDROPHOBIC_RESIDUES})",
            ),
        ],
    },
    {
        "name": "v251",
        "description": "V251 in TOA1",
        "roi_selection": V251["roi_selection"],
        "core_selection": V251["core_selection"],
        "actions": [
            lambda: cmd.select(
                V251["roi_selection"],
                f"chain {V251['chain']} and resi {V251['residue_number']}",
            ),
            lambda: cmd.select(
                V251["nearby_selection"],
                f"byres (all within {DISTANCE_CUTOFF} of {V251['roi_selection']})",
            ),
            lambda: cmd.select(
                V251["core_selection"],
                f"{V251['nearby_selection']} and resn {HYDROPHOBIC_RESIDUES} and not {V251['roi_selection']}",
            ),
            lambda: cmd.select(
                V251["nearby_contact_selection"],
                f"byres (chain {V251['contact_chain']} within {DISTANCE_CUTOFF} of {V251['core_selection']})",
            ),
            lambda: cmd.select(
                V251["core_selection"],
                f"{V251['core_selection']} or ({V251['nearby_contact_selection']} and resn {HYDROPHOBIC_RESIDUES})",
            ),
        ],
    },
]

# --- 1nh2/TFIIA-TBP-DNA complex views ---
# View matrices (capture these with get_view() in PyMOL GUI)
VIEW_OVERALL = (
    0.074909329,
    0.990534842,
    -0.115008026,
    0.582228422,
    0.050188679,
    0.811473668,
    0.809565485,
    -0.127748311,
    -0.572958469,
    0.000000000,
    0.000000000,
    -246.975219727,
    5.988648415,
    2.612197876,
    0.841152191,
    194.940582275,
    299.009948730,
    -20.000000000,
)

# --- G16 views ---
# Previous view. Changed angle of view from overall.
VIEW_HYDROPHOBIC = (
    -0.275569350,
    0.849549353,
    -0.449802339,
    0.956983745,
    0.198277250,
    -0.211805791,
    -0.090754412,
    -0.488821685,
    -0.867647946,
    0.000000000,
    0.000000000,
    -148.670715332,
    -1.743014336,
    -25.174198151,
    -14.948152542,
    117.213119507,
    180.128311157,
    -20.000000000,
)

# Same angle as view overall.
VIEW_G16_HYDROPHOBIC = (
    0.074909329,
    0.990534842,
    -0.115008026,
    0.582228422,
    0.050188679,
    0.811473668,
    0.809565485,
    -0.127748311,
    -0.572958469,
    0.000082036,
    -0.000030905,
    -122.774742126,
    -3.036495209,
    -19.484025955,
    -28.640813828,
    73.339439392,
    173.608703613,
    -20.000000000,
)

# --- V251 views ---
# Intentional alias - both ROIs are presented from the same full-complex angle
VIEW_V251_OVERALL = VIEW_OVERALL

VIEW_V251_HYDROPHOBIC = (
    0.369137824,
    0.765394688,
    -0.527172744,
    0.796817005,
    0.031308878,
    0.603406489,
    0.478349477,
    -0.642804384,
    -0.598319829,
    -0.000008364,
    0.000012737,
    -81.610771179,
    4.645626068,
    0.182135344,
    -0.962601185,
    -24.776704788,
    188.000106812,
    -20.000000000,
)

# ============================================================================
# SECTION 2: INITIALIZE SESSION
# ============================================================================

print("=" * 60)
print(f"PyMOL Visualization Session: {PDB_CODE}")
print("=" * 60)

# Apply shared settings
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
for chain, color in CHAIN_COLORS.items():
    selection = f"chain {chain}"
    if cmd.count_atoms(selection) > 0:
        cmd.color(color, selection)
        print(f"  Chain {chain}: {color}")


# ============================================================================
# SECTION 5: CREATE CUSTOM SELECTIONS (FOR ALL ROIs)
# ============================================================================

print("\n" + "=" * 60)
print("Creating custom selections for all ROIs")
print("=" * 60)

for roi_spec in ROI_SPECS:
    print(f"\n{roi_spec['description']}:")

    # Execute all selection actions for this ROI
    for action in roi_spec["actions"]:
        action()

    # Verify the main selections were created successfully
    roi_selection_name = roi_spec["roi_selection"]
    core_selection_name = roi_spec["core_selection"]
    roi_atom_count = cmd.count_atoms(roi_selection_name)
    if roi_atom_count == 0:
        raise ValueError(
            f"ERROR: Selection {roi_selection_name} is empty! Check chain/residue in ROI_SPECS."
        )
    core_atom_count = cmd.count_atoms(core_selection_name)
    print(f"   {roi_selection_name}: {roi_atom_count} atoms")
    print(f"   {core_selection_name}: {core_atom_count} atoms")

print(f"\n Successfully created selections for {len(ROI_SPECS)} ROIs")
# ============================================================================
# SECTION 6: DEFINE IMAGE SPECIFICATIONS
# ============================================================================


image_specs = [
    # --- G16 images ---
    {
        "filename": "01_g16_overall_view.png",
        "description": "G16 full complex with ROI highlighted",
        "view": VIEW_OVERALL,
        "rotation_angles": [0, 90, 180, 270],
        "actions": [
            lambda: helper_functions.reset_scene_to_base_cartoon(CHAIN_COLORS),
            lambda: cmd.show("spheres", G16["roi_selection"]),
            lambda: cmd.color("red", G16["roi_selection"]),
            lambda: cmd.set("sphere_scale", 1.5, G16["roi_selection"]),
        ],
    },
    {
        "filename": "02_g16_hydrophobic_core_surface.png",
        "description": "G16 hydrophobic core as surface",
        "view": VIEW_G16_HYDROPHOBIC,
        "actions": [
            lambda: helper_functions.reset_scene_to_base_cartoon(CHAIN_COLORS),
            lambda: cmd.set("depth_cue", 0),
            lambda: cmd.set("fog", 0),
            lambda: cmd.show("spheres", G16["core_selection"]),
            lambda: cmd.set("sphere_scale", 1.0, G16["core_selection"]),
            lambda: cmd.color("grey50", G16["core_selection"]),
            lambda: cmd.set("transparency", 0.8, G16["core_selection"]),
        ],
    },
    # --- V251 images (NEW - add these) ---
    {
        "filename": "03_v251_overall_view.png",
        "description": "V251 supplemental - full complex with ROI highlighted",
        "view": VIEW_V251_OVERALL,  # Reuses same overall view
        "rotation_angles": [0, 90, 180, 270],
        "actions": [
            lambda: helper_functions.reset_scene_to_base_cartoon(CHAIN_COLORS),
            lambda: cmd.color(
                "marine", "obj 1nh2 and chain B"
            ),  # Reapply color to remove orange.
            lambda: cmd.color("green", "obj 1nh2 and chain D"),
            lambda: cmd.show("spheres", V251["roi_selection"]),
            lambda: cmd.color("magenta", V251["roi_selection"]),
            lambda: cmd.show("spheres", G16["roi_selection"]),
            lambda: cmd.color("red", G16["roi_selection"]),
            lambda: cmd.set("sphere_scale", 1.5, G16["roi_selection"]),
            lambda: cmd.set("sphere_scale", 1.5, V251["roi_selection"]),
        ],
    },
    {
        "filename": "04_v251_hydrophobic_core_surface.png",
        "description": "V251 supplemental - hydrophobic core with TBP residues 94-106 hidden",
        "view": VIEW_V251_HYDROPHOBIC,
        "actions": [
            lambda: helper_functions.reset_scene_to_base_cartoon(CHAIN_COLORS),
            lambda: cmd.show("spheres", V251["core_selection"]),
            lambda: cmd.set("sphere_scale", 1.0, V251["core_selection"]),
            lambda: cmd.color("grey50", V251["core_selection"]),
            lambda: cmd.color("magenta", V251["roi_selection"]),
            lambda: cmd.set("transparency", 0.3, V251["core_selection"]),
            # Special: hide TBP residues 94-106 from chain A
            lambda: cmd.hide("cartoon", "chain A and resi 94-106"),
        ],
    },
]
# ============================================================================
# SECTION 7: RENDER ALL IMAGES
# ============================================================================

print("\n" + "=" * 60)
print("Rendering images")
print("=" * 60)

for i, spec in enumerate(image_specs, 1):
    print(f"\nImage spec {i}/{len(image_specs)}: {spec['filename']}")
    print(f"  Description: {spec['description']}")

    # Execute actions
    for action in spec["actions"]:
        action()

    # Set view
    cmd.set_view(spec["view"])

    # Check if this should be a rotation series
    if "rotation_angles" in spec:
        angles = spec["rotation_angles"]
        base_filename = spec["filename"].replace(".png", "")
        print(f"  Generating {len(angles)} rotation views...")
        for angle_index, angle in enumerate(angles):
            angle_filename = f"{base_filename}_{angle:03d}.png"
            print(f"    Angle {angle}...")
            helper_functions.save_png(
                filename=angle_filename,
                width=pymol_configuration.IMAGE_WIDTH,
                height=pymol_configuration.IMAGE_HEIGHT,
                dpi=pymol_configuration.IMAGE_DPI,
                ray=pymol_configuration.RAY_TRACE,
                overwrite=pymol_configuration.FLAG_OVERWRITE,
                output_dir=pymol_configuration.OUTPUT_DIR,
            )
            is_last_angle = angle_index == len(angles) - 1
            if not is_last_angle:
                next_angle = angles[angle_index + 1]
                rotation_delta = next_angle - angle
                cmd.turn("y", rotation_delta)
    else:
        # Single image (no rotation)
        helper_functions.save_png(
            filename=spec["filename"],
            width=pymol_configuration.IMAGE_WIDTH,
            height=pymol_configuration.IMAGE_HEIGHT,
            dpi=pymol_configuration.IMAGE_DPI,
            ray=pymol_configuration.RAY_TRACE,
            overwrite=pymol_configuration.FLAG_OVERWRITE,
            output_dir=pymol_configuration.OUTPUT_DIR,
        )

# ============================================================================
# SECTION 8: CREATE CHAIN OBJECTS FOR APBS CALCULATIONS
# ============================================================================
print("\n" + "=" * 60)
print("Creating individual chain objects for APBS")
print("=" * 60)

# NOTE: The downstream APBS script (toa2_hydrophobic_core_pictures_after_apbs.py)
# references IMAGE_WIDTH, IMAGE_HEIGHT, IMAGE_DPI, RAY_TRACE, OUTPUT_DIR, and
# FLAG_OVERWRITE directly from this session's namespace. Both scripts must be
# run in the same PyMOL session for those bindings to be available.

# Use PDB_CODE directly as the object name - cmd.fetch() names the object after the PDB code
structure_object_name = PDB_CODE

# Create separate object for each chain
chain_objects = []
for chain in CHAIN_COLORS.keys():
    chain_obj_name = f"chain_{chain}"
    selection = f"{structure_object_name} and chain {chain}"

    if cmd.count_atoms(selection) > 0:
        cmd.create(chain_obj_name, selection)
        chain_objects.append(chain_obj_name)
        print(f"  Created: {chain_obj_name}")

# Hide chain objects (keep main structure visible)
for chain_obj in chain_objects:
    cmd.disable(chain_obj)

print(f"\n Created {len(chain_objects)} chain objects (currently hidden)")

print("\n" + "=" * 60)
print("NEXT STEPS: Run APBS on each chain object")
print("=" * 60)
print("Instructions:")
for i, chain_obj in enumerate(chain_objects, 1):
    print(f"\n  Run {i}: {chain_obj}")
    print(f"    1. Enable: {chain_obj}")
    print(f"    2. Plugins > APBS Tools")
    print(f"    3. Select '{chain_obj}' as molecule")
    print(f"    4. Run PDB2PQR > Run APBS > Visualize")

print("\n" + "=" * 60)
print("After APBS: Run the next script to visualize electrostatics")
print("=" * 60)

# ============================================================================
# SECTION 9: SAVE SESSION AND METADATA
# ============================================================================

print("\n" + "=" * 60)
print("Saving session files")
print("=" * 60)

# Save PyMOL session
import os

session_file = os.path.join(pymol_configuration.OUTPUT_DIR, f"{PDB_CODE}_session.pse")
try:
    cmd.save(session_file)
    print(f"  Saved PyMOL session: {session_file}")
except Exception as e:
    print(f"  ERROR saving session: {e}")

# ============================================================================
# SECTION 10: FINAL SUMMARY
# ============================================================================

print("\n" + "=" * 60)
print("SESSION COMPLETE")
print("=" * 60)
print(f"Structure: {PDB_CODE}")
print(f"Generated {len(image_specs)} images in: {pymol_configuration.OUTPUT_DIR}")
print("=" * 60)
