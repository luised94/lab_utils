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

# Structure and mutation details
#STRUCTURE_FILE = "your_structure.pdb"  # or use fetch code
PDB_CODE = "1nh2"
TARGET_CHAIN = "D"
TARGET_RESI = 16
WILD_TYPE = "GLY"
MUTANT = "GLU"

# Distance for selecting hydrophobic core (Angstroms)
DISTANCE_CUTOFF = 5.0

# Hydrophobic amino acids definition
HYDROPHOBIC_SET = "ala+gly+val+ile+leu+phe+met+pro+trp"

# Representation for hydrophobic core: "surface", "spheres", or "sticks"
CORE_REPRESENTATION = "surface"  # Change as needed

# Chain color scheme (use PyMOL color names)
CHAIN_COLORS = {
    "A": "marine", # tbp
    "B": "green", # toa1
    "C": "green", # toa1, second chain
    "D": "yellow" # toa2
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

print("STOPPING HERE FOR REVIEW")
sys.exit()
