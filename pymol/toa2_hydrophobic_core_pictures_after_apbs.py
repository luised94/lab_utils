# Take hydrophobic core pictures of TFIIA after preparation.
# Show surface of hydrophobic core.
# Usage:
# File > Run script...
# Or place the script in the working directory (File > Working Directory > Change)
# run toa2_hydrophobic_core_pictures_after_apbs.py

from pymol import cmd
import os
import sys

# Import user functions
import helper_functions

# Configuration
FLAG_OVERWRITE = True

# ============================================================================
# VALIDATION: Check Prerequisites Before Running
# ============================================================================

print("=" * 60)
print("Validating Prerequisites for Electrostatic Visualization")
print("=" * 60)

# Define expected chains (must match APBS run order)
# Order matters! First APBS run = chain D, second = chain B
EXPECTED_CHAINS = ["B", "D"] # B is toa1, fragment 1 and D is toa2.

# ============================================================================
# Debug: See what PyMOL has
# ============================================================================

print("\nAll objects:")
print(cmd.get_object_list())

print("\nAll names (including groups):")
print(cmd.get_names("all"))

print("\nJust objects:")
print(cmd.get_names("objects"))

print("\nPublic objects:")
print(cmd.get_names("public_objects"))

# ============================================================================
# Check that chain objects exist
# ============================================================================

all_objects = cmd.get_object_list()
missing_chains = []

for chain in EXPECTED_CHAINS:
    chain_obj = f"chain_{chain}"
    if chain_obj not in all_objects:
        missing_chains.append(chain_obj)

if missing_chains:
    print("ERROR: Missing chain objects:")
    for obj in missing_chains:
        print(f"  - {obj}")
    print("\nFAILURE: Run prepare_tfiia.py first to create chain objects")
    sys.exit(1)

print(f" All chain objects found: {EXPECTED_CHAINS}")

# ============================================================================
# Check for prepared objects (APBS creates preparedXX)
# ============================================================================

all_names = cmd.get_names("all")

# Find prepared objects
prepared_objects = [obj for obj in all_names if obj.startswith("prepared")]
prepared_objects.sort()  # prepared01, prepared02, etc.

print(f"\nFound {len(prepared_objects)} prepared objects: {prepared_objects}")

if len(prepared_objects) < len(EXPECTED_CHAINS):
    print(f"ERROR: Expected {len(EXPECTED_CHAINS)} APBS runs, found {len(prepared_objects)}")
    print("\nFAILURE: Complete APBS calculations for all chains")
    print("\nExpected order:")
    for i, chain in enumerate(EXPECTED_CHAINS, start=1):
        print(f"  Run {i}: chain_{chain}")
    sys.exit(1)

# ============================================================================
# Map prepared objects to chains and find associated maps
# ============================================================================

apbs_mapping = {}

print("\nAPBS Objects Found:")
for i, chain in enumerate(EXPECTED_CHAINS, start=1):
    run_num = f"{i:02d}"  # 01, 02, etc.

    prepared_obj = f"prepared{run_num}"
    map_obj = f"apbs_map{run_num}"
    ramp_obj = f"apbs_ramp{run_num}"

    # Check if objects exist
    if prepared_obj not in all_names:
        print(f"  WARNING: Missing {prepared_obj} for chain {chain}")
        continue

    if map_obj not in all_names:
        print(f"  WARNING: Missing {map_obj} for chain {chain}")
        continue

    apbs_mapping[chain] = {
        'prepared': prepared_obj,
        'map': map_obj,
        'ramp': ramp_obj if ramp_obj in all_names else None
    }

    print(f"  Chain {chain}: {prepared_obj}, {map_obj}, {ramp_obj}")

if len(apbs_mapping) < len(EXPECTED_CHAINS):
    print(f"\nERROR: Incomplete APBS data")
    print("\nFAILURE: Complete APBS calculations for all chains")
    sys.exit(1)

print(f" All APBS maps found for chains: {list(apbs_mapping.keys())}")

print("\n All prerequisites satisfied")
print("=" * 60)

cmd.set_view((
    -0.275569350,    0.849549353,   -0.449802339,
     0.956983745,    0.198277250,   -0.211805791,
    -0.090754412,   -0.488821685,   -0.867647946,
     0.000177816,   -0.000032090, -148.673400879,
     2.425498962,  -25.367309570,   -6.927139282,
   117.213119507,  180.128311157,  -20.000000000))


# ============================================================================
# ELECTROSTATIC SURFACE VISUALIZATION
# ============================================================================

print("\n" + "=" * 60)
print("Creating Electrostatic Surface Views")
print("=" * 60)

# Enable all necessary objects first
# Enable all objects programmatically
#print("\nEnabling all objects...")
#all_objects = cmd.get_object_list()
#for obj in all_objects:
#    cmd.enable(obj)
#print(f"Enabled {len(all_objects)} objects")
print("\nEnabling all objects...")
cmd.enable("1nh2")  # Original structure with hydrophobic_core selection
cmd.enable("chain_B")
cmd.enable("chain_D")
cmd.enable("prepared01")
cmd.enable("prepared02")

# Image counter for sequential naming
img_counter = 1

# ============================================================================
# View 1: Chain D interface (90 degrees)
# ============================================================================

print("\nView 1: Chain D interface (90 degrees)")

# Turn 90 degrees
#cmd.turn("y", 90)
cmd.set_view((
    -0.450920254,    0.852479875,    0.264468431,
    -0.212777704,    0.185093418,   -0.959405780,
    -0.866829515,   -0.488889456,    0.097926863,
    -0.000032587,   -0.000011193, -154.020599365,
    -2.366145849,  -16.602687836,   -7.165100098,
   105.576972961,  197.627319336,  -20.000000000 ))

# Hide everything
cmd.hide("everything")

# Show chain_D cartoon, save
cmd.show("cartoon", "chain_D")
helper_functions.save_png(
    filename=f"{img_counter:02d}_cartoon_D.png",
    width=IMAGE_WIDTH,
    height=IMAGE_HEIGHT,
    dpi=300,
    ray=RAY_TRACE,
    overwrite=FLAG_OVERWRITE,
    output_dir=OUTPUT_DIR
)
img_counter += 1

# Show prepared02 (chain_D electrostatics), save
cmd.show("surface", "prepared02")
helper_functions.save_png(
    filename=f"{img_counter:02d}_surface_D.png",
    width=IMAGE_WIDTH,
    height=IMAGE_HEIGHT,
    dpi=300,
    ray=RAY_TRACE,
    overwrite=FLAG_OVERWRITE,
    output_dir=OUTPUT_DIR
)
img_counter += 1

# Show chain_B cartoon, save
cmd.show("cartoon", "chain_B")
helper_functions.save_png(
    filename=f"{img_counter:02d}_surface_D_cartoon_B.png",
    width=IMAGE_WIDTH,
    height=IMAGE_HEIGHT,
    dpi=300,
    ray=RAY_TRACE,
    overwrite=FLAG_OVERWRITE,
    output_dir=OUTPUT_DIR
)
img_counter += 1

# Show hydrophobic core (selection exists in original 1nh2 object)
cmd.show("spheres", "hydrophobic_core")
cmd.set("sphere_scale", 0.6, "hydrophobic_core")
cmd.color("orange", "hydrophobic_core")
cmd.set("sphere_transparency", 0.2, "hydrophobic_core")
helper_functions.save_png(
    filename=f"{img_counter:02d}_surface_D_cartoon_B_hydrophobic.png",
    width=IMAGE_WIDTH,
    height=IMAGE_HEIGHT,
    dpi=300,
    ray=RAY_TRACE,
    overwrite=FLAG_OVERWRITE,
    output_dir=OUTPUT_DIR
)
img_counter += 1

# ============================================================================
# View 2: Chain B interface (270 degrees = turn -180 from current)
# ============================================================================

print("\nView 2: Chain B interface (270 degrees)")

# Turn -180 degrees
#cmd.turn("y", -180)
cmd.set_view((
     0.527584910,    0.820460081,    0.220202491,
     0.216749787,   -0.380644381,    0.898955524,
     0.821382821,   -0.426547915,   -0.378658950,
     0.000019502,   -0.000110969, -119.689163208,
    -2.949608326,  -20.330551147,  -17.304233551,
    75.965362549,  163.015731812,  -20.000000000 ))

# Hide everything
cmd.hide("everything")

# Show chain_B cartoon, save
cmd.show("cartoon", "chain_B")
helper_functions.save_png(
    filename=f"{img_counter:02d}_cartoon_B.png",
    width=IMAGE_WIDTH,
    height=IMAGE_HEIGHT,
    dpi=300,
    ray=RAY_TRACE,
    overwrite=FLAG_OVERWRITE,
    output_dir=OUTPUT_DIR
)
img_counter += 1

# Show prepared01 (chain_B electrostatics), save
cmd.show("surface", "prepared01")
helper_functions.save_png(
    filename=f"{img_counter:02d}_surface_B.png",
    width=IMAGE_WIDTH,
    height=IMAGE_HEIGHT,
    dpi=300,
    ray=RAY_TRACE,
    overwrite=FLAG_OVERWRITE,
    output_dir=OUTPUT_DIR
)
img_counter += 1

# Show chain_D cartoon, save
cmd.show("cartoon", "chain_D")
helper_functions.save_png(
    filename=f"{img_counter:02d}_surface_B_cartoon_D.png",
    width=IMAGE_WIDTH,
    height=IMAGE_HEIGHT,
    dpi=300,
    ray=RAY_TRACE,
    overwrite=FLAG_OVERWRITE,
    output_dir=OUTPUT_DIR
)
img_counter += 1

# Show hydrophobic core (selection exists in original 1nh2 object)
cmd.show("spheres", "hydrophobic_core")
cmd.set("sphere_scale", 0.6, "hydrophobic_core")
cmd.color("orange", "hydrophobic_core")
cmd.set("sphere_transparency", 0.2, "hydrophobic_core")
helper_functions.save_png(
    filename=f"{img_counter:02d}_surface_B_cartoon_D_hydrophobic.png",
    width=IMAGE_WIDTH,
    height=IMAGE_HEIGHT,
    dpi=300,
    ray=RAY_TRACE,
    overwrite=FLAG_OVERWRITE,
    output_dir=OUTPUT_DIR
)
img_counter += 1

print("\n" + "=" * 60)
print(f"Electrostatic visualization complete! Generated {img_counter - 1} images.")
print("=" * 60)
