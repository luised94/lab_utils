"""
PyMOL Visualization Session: TFIIA-TBP-DNA Complex (1NH2)
Purpose: Generate publication-quality images showing hydrophobic core
Author: LEMR
Date: 2025-12-08
"""

from pymol import cmd
import pymol_configuration
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
TARGET_CHAIN = "D"
TARGET_RESI = 16

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

# View matrices (capture these with get_view() in PyMOL GUI)
VIEW_OVERALL = (
     0.074909329,    0.990534842,   -0.115008026,
     0.582228422,    0.050188679,    0.811473668,
     0.809565485,   -0.127748311,   -0.572958469,
     0.000000000,    0.000000000, -246.975219727,
     5.988648415,    2.612197876,    0.841152191,
   194.940582275,  299.009948730,  -20.000000000
)

VIEW_HYDROPHOBIC = (
    -0.275569350,    0.849549353,   -0.449802339,
     0.956983745,    0.198277250,   -0.211805791,
    -0.090754412,   -0.488821685,   -0.867647946,
     0.000000000,    0.000000000, -148.670715332,
    -1.743014336,  -25.174198151,  -14.948152542,
   117.213119507,  180.128311157,  -20.000000000
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

# Load structure
print(f"\nLoading structure {PDB_CODE}...")
cmd.fetch(PDB_CODE, type = 'pdb')

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
# SECTION 5: CREATE CUSTOM SELECTIONS
# ============================================================================

print("\nCreating custom selections...")

cmd.select("roi", f"chain {TARGET_CHAIN} and resi {TARGET_RESI}")
roi_count = cmd.count_atoms("roi")
print(f"  ROI (chain {TARGET_CHAIN}, resi {TARGET_RESI}): {roi_count} atoms")

if roi_count == 0:
    raise ValueError(f"ERROR: No atoms found for ROI!")

cmd.select("nearby", f"byres (all within {DISTANCE_CUTOFF} of roi)")
print(f"  Nearby residues: {cmd.count_atoms('nearby')} atoms")

cmd.select("hydrophobic_core", f"nearby and resn {HYDROPHOBIC_RESIDUES} and not roi")
print(f"  Hydrophobic core: {cmd.count_atoms('hydrophobic_core')} atoms")

# Expand again but select around chain b (Toa1) and get union.
cmd.select("nearby_chain_b", f"byres (chain B within {DISTANCE_CUTOFF} of hydrophobic_core)")
cmd.select("hydrophobic_core", f"hydrophobic_core or (nearby_chain_b and resn {HYDROPHOBIC_SET})")

# Color ROI red
cmd.color("red", "roi")
cmd.show("spheres", "roi")
cmd.set("sphere_scale", 1.5, "roi")
print("Colored ROI red")

# ============================================================================
# SECTION 6: DEFINE IMAGE SPECIFICATIONS
# ============================================================================

image_specs = [
    {
        'filename': '01_overall_view.png',
        'description': 'Full complex with ROI highlighted',
        'view': VIEW_OVERALL,
        'actions': [
            lambda: cmd.hide("everything"),
            lambda: cmd.show("cartoon", "polymer"),
            lambda: cmd.show("spheres", "roi"),
            lambda: cmd.color("red", "roi"),
            lambda: cmd.set("sphere_scale", 1.5, "roi"),
        ]
    },

    {
        'filename': '02_hydrophobic_core_surface.png',
        'description': 'Hydrophobic core shown as surface',
        'view': VIEW_HYDROPHOBIC,
        'actions': [
            lambda: cmd.hide("everything"),
            lambda: cmd.show("cartoon", "polymer"),
            lambda: cmd.show("spheres", "hydrophobic_core"),
            lambda: cmd.color("orange", "hydrophobic_core"),
            lambda: cmd.set("transparency", 0.8, "hydrophobic_core"),
        ]
    },
]

# ============================================================================
# SECTION 7: RENDER ALL IMAGES
# ============================================================================

print("\n" + "=" * 60)
print("Rendering images")
print("=" * 60)

for i, spec in enumerate(image_specs, 1):
    print(f"\nImage {i}/{len(image_specs)}: {spec['filename']}")
    print(f"  Description: {spec['description']}")

    # Execute actions
    for action in spec['actions']:
        action()

    # Set view
    cmd.set_view(spec['view'])

    # Save image
    helper_functions.save_png(
        filename=spec['filename'],
        width=pymol_configuration.IMAGE_WIDTH,
        height=pymol_configuration.IMAGE_HEIGHT,
        dpi=pymol_configuration.IMAGE_DPI,
        ray=pymol_configuration.RAY_TRACE,
        overwrite=pymol_configuration.FLAG_OVERWRITE,
        output_dir=pymol_configuration.OUTPUT_DIR
    )

# ============================================================================
# SECTION 8: SAVE SESSION AND METADATA
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
# SECTION 9: FINAL SUMMARY
# ============================================================================

print("\n" + "=" * 60)
print("SESSION COMPLETE")
print("=" * 60)
print(f"Structure: {PDB_CODE}")
print(f"Generated {len(image_specs)} images in: {pymol_configuration.OUTPUT_DIR}")
print("=" * 60)
