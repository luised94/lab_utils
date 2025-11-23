# Prepare tfiia structure for publication
# Show hydrophobic core near G16.
# Usage:
# File > Run script...
# Or place the script in the working directory (File > Working Directory > Change)
# run prepare_tfiia-tbp-dna_1nh2_structure.py
from pymol import cmd
import os
import sys

# Import user functions
import helper_functions
#from helper_functions import save_png

# ============================================================================
# START
# ============================================================================
cmd.reinitialize()

# ============================================================================
# CONFIGURATION - Modify these values for your specific case
# ============================================================================
FLAG_OVERWRITE = True
BACKGROUND_COLOR = "white"
HASH_MAX = 200  # For better ray-traced image quality

# Background color
cmd.bg_color(BACKGROUND_COLOR)

# Structure and mutation details
#STRUCTURE_FILE = "your_structure.pdb"  # or use fetch code
PDB_CODE = "1nh2"
TARGET_CHAIN = "D"
TARGET_RESI = 16
WILD_TYPE = "GLY"
MUTANT = "GLU"

# Distance for selecting hydrophobic core (Angstroms)
DISTANCE_CUTOFF = 4.0

# Hydrophobic amino acids definition
HYDROPHOBIC_SET = "ala+gly+val+ile+leu+phe+met+pro+trp"

# Representation for hydrophobic core: "surface", "spheres", or "sticks"
CORE_REPRESENTATION = "surface"  # Change as needed

# Chain color scheme (use PyMOL color names)
CHAIN_COLORS = {
    "A": "yellow", # tbp
    "B": "marine", # toa1
    "C": "marine", # toa1, second chain
    "D": "green" # toa2
    #"E": "orange", # DNA 1
    #"F": "orange", # DNA 2
}

# Mutation site colors
WT_COLOR = "yellow"
MUTANT_COLOR = "red"
CORE_COLOR = "sand"  # Hydrophobic core color

# Export settings
OUTPUT_DIR = "pymol_outputs"
IMAGE_WIDTH = 1200
IMAGE_HEIGHT = 900
RAY_TRACE = True  # Set to False for quick preview

# ============================================================================
# SETUP AND PREPARATION
# ============================================================================

print("=" * 60)
print("Starting PyMOL Visualization Script")
print("=" * 60)

# Create output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"Created output directory: {OUTPUT_DIR}")

# Load structure (uncomment appropriate line)
# cmd.load(STRUCTURE_FILE)
cmd.fetch(PDB_CODE)  # Or use fetch for PDB entries

# Get the object name (assumes first loaded object)
obj_name = cmd.get_object_list()[0]
print(f"Working with object: {obj_name}")

# Clean up structure
cmd.remove("solvent")
cmd.remove("hydrogen")
cmd.remove("resn HOH")
print("Removed solvent and waters")

# ============================================================================
# COLOR CHAINS
# ============================================================================

print("\nColoring chains...")
for chain, color in CHAIN_COLORS.items():
    selection = f"{obj_name} and chain {chain}"
    if cmd.count_atoms(selection) > 0:
        cmd.color(color, selection)
        print(f"  Chain {chain}: {color}")

# Set cartoon representation for proteins
cmd.show("cartoon", f"{obj_name} and polymer")
cmd.hide("lines", f"{obj_name} and polymer")

# ============================================================================
# IDENTIFY RESIDUE OF INTEREST AND HYDROPHOBIC CORE
# ============================================================================

print("\nIdentifying hydrophobic core...")

# Select residue of interest
cmd.select("roi", f"chain {TARGET_CHAIN} and resi {TARGET_RESI}")
roi_count = cmd.count_atoms("roi")
print(f"  Residue of interest (chain {TARGET_CHAIN}, resi {TARGET_RESI}): {roi_count} atoms")

if roi_count == 0:
    print("ERROR: Residue of interest not found!")
    raise ValueError(f"No atoms found for chain {TARGET_CHAIN} resi {TARGET_RESI}")

# Select nearby residues
cmd.select("nearby", f"byres ({obj_name} within {DISTANCE_CUTOFF} of roi)")
nearby_count = cmd.count_atoms("nearby")
print(f"  Nearby residues (within {DISTANCE_CUTOFF}): {nearby_count} atoms")

# Select hydrophobic residues in the core
cmd.select("hydrophobic_core", f"nearby and resn {HYDROPHOBIC_SET} and not roi")
core_count = cmd.count_atoms("hydrophobic_core")
print(f"  Hydrophobic core residues: {core_count} atoms")

# Expand again but select around chain b (Toa1) and get union.
cmd.select("nearby_chain_b", f"byres (chain B within {DISTANCE_CUTOFF} of hydrophobic_core)")
cmd.select("hydrophobic_core", f"hydrophobic_core or (nearby_chain_b and resn {HYDROPHOBIC_SET})")

# Color ROI red
cmd.color("red", "roi")
cmd.show("spheres", "roi")
cmd.set("sphere_scale", 1.5, "roi")
print("Colored ROI red")

# ============================================================================
# SCENE 1: OVERALL VIEW WITH 360 ROTATION
# ============================================================================

print("\n" + "=" * 60)
print("Creating Scene 1: Overall View with 360 Rotation")
print("=" * 60)

cmd.set_view((
     0.074909329,    0.990534842,   -0.115008026,
     0.582228422,    0.050188679,    0.811473668,
     0.809565485,   -0.127748311,   -0.572958469,
     0.000000000,    0.000000000, -246.975219727,
     5.988648415,    2.612197876,    0.841152191,
   194.940582275,  299.009948730,  -20.000000000 ))

# Better ray-tracing settings
cmd.set("hash_max", HASH_MAX)
cmd.set("antialias", 2)
cmd.set("ray_shadows", 1)
#cmd.set("ray_shadows", 1)
cmd.set("ray_shadow", 0.2)  # Default is 0.5, lower = lighter shadows
cmd.set("ray_shadow_decay_factor", 0.1)  # How fast shadow fades
cmd.set("ray_shadow_decay_range", 2)
cmd.set("ray_trace_fog", 0)
cmd.set("depth_cue", 0)
cmd.set("fog", 0)

# Generate 360 rotation images (every 90 degrees)
angles = [0, 90, 180, 270]
for angle in angles:
    print(f"Capturing view at {angle}...")

    helper_functions.save_png(
        filename = f"01_overall_view_{angle:03d}.png",
        width = IMAGE_WIDTH,
        height = IMAGE_HEIGHT,
        dpi = 300,
        ray = RAY_TRACE,
        overwrite = FLAG_OVERWRITE,
        output_dir = OUTPUT_DIR
    )

    if angle < 270:  # Don't rotate after last image
        cmd.turn("y", 90)

print(f"Generated {len(angles)} rotation views")

# See through the helix bundle.
cmd.set_view((
    -0.275569350,    0.849549353,   -0.449802339,
     0.956983745,    0.198277250,   -0.211805791,
    -0.090754412,   -0.488821685,   -0.867647946,
     0.000177816,   -0.000032090, -148.673400879,
     2.425498962,  -25.367309570,   -6.927139282,
   117.213119507,  180.128311157,  -20.000000000))
# View 2
#cmd.set_view((
#    -0.275569350,    0.849549353,   -0.449802339,
#     0.956983745,    0.198277250,   -0.211805791,
#    -0.090754412,   -0.488821685,   -0.867647946,
#     0.000000000,    0.000000000, -148.670715332,
#    -1.743014336,  -25.174198151,  -14.948152542,
#   117.213119507,  180.128311157,  -20.000000000))

# ============================================================================
# Take picture of four helix bundle, looking down the center
# ============================================================================
print("\n" + "=" * 60)
print("Taking four helix bundle hydrophobic_core sphere picture")
print("=" * 60)

cmd.hide("spheres", "roi")
cmd.set("sphere_scale", 0.8, "hydrophobic_core")
cmd.show("spheres", "hydrophobic_core")
cmd.color("orange", "hydrophobic_core")

helper_functions.save_png(
    filename = "tfiia_1nh2_hydrophobic_core_sphere.png",
    width = IMAGE_WIDTH,
    height = IMAGE_HEIGHT,
    dpi = 300,
    ray = RAY_TRACE,
    overwrite = FLAG_OVERWRITE,
    output_dir = OUTPUT_DIR
)

# ============================================================================
# CREATE INDIVIDUAL CHAIN OBJECTS FOR APBS CALCULATIONS
# ============================================================================

print("\n" + "=" * 60)
print("Creating individual chain objects for APBS")
print("=" * 60)

# Create separate object for each chain
chain_objects = []
for chain in CHAIN_COLORS.keys():
    chain_obj_name = f"chain_{chain}"
    selection = f"{obj_name} and chain {chain}"
    if cmd.count_atoms(selection) > 0:
        cmd.create(chain_obj_name, selection)
        chain_objects.append(chain_obj_name)
        print(f"Created: {chain_obj_name}")

# Disable chain objects (keep them hidden for now)
for chain_obj in chain_objects:
    cmd.disable(chain_obj)
print(f"\nCreated {len(chain_objects)} chain objects (currently hidden)")

print("\n" + "=" * 60)
print("NEXT STEPS: Run APBS on each chain object")
print("=" * 60)
print("For each chain object:")
for chain_obj in chain_objects:
    print(f"  1. Enable: {chain_obj}")
    print(f"  2. Plugins > APBS Tools")
    print(f"  3. Select '{chain_obj}' as molecule")
    print(f"  4. Run PDB2PQR > Run APBS > Visualize")
print("=" * 60)

print("\nPreparation script complete!")
print(f"Main object: {obj_name}")
print(f"Chain objects: {', '.join(chain_objects)}")
