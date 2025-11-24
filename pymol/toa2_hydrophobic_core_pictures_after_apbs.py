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

# Turn 90 degrees
# Hide everything.
# Show chain_D, save.
# Show run02, save.
# Show chain_B, save.
# Show chain_B hydrophobic, save.
# Turn y -180 degreers.
# Hide everything.
# Show chain_B, save.
# show run02, save.
# Show chain_D, save.
# Show chain_B, hydrophobic, save.

# ============================================================================
# Use the mapping in visualization
# ============================================================================

# Example usage:
# chain_D_map = apbs_mapping['D']['map']  # 'apbs_map01'
# chain_B_map = apbs_mapping['B']['map']  # 'apbs_map02'
