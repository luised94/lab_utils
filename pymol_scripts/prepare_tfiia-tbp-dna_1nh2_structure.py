# Prepare tfiia structure for publication
# Show hydrophobic core near G16.
# Usage:
# File > Run script...
# Or place the script in the working directory (File > Working Directory > Change)
# run prepare_tfiia-tbp-dna_1nh2_structure.py
from pymol import cmd
import os
import sys
cmd.reinitialize()

# ============================================================================
# CONFIGURATION - Modify these values for your specific case
# ============================================================================
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

print("STOPPING HERE FOR REVIEW")
sys.exit()
