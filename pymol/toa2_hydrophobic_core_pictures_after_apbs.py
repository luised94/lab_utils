# ============================================================================
# VALIDATION: Check Prerequisites Before Running
# ============================================================================

print("=" * 60)
print("Validating Prerequisites for Electrostatic Visualization")
print("=" * 60)

# Define expected chains (must match APBS run order)
EXPECTED_CHAINS = ["D", "B"]  # Order matters! First APBS run = chain D, second = chain B

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
# Check that APBS maps exist and map them to chains
# ============================================================================

# Find all run objects (run01, run02, etc.)
run_objects = [obj for obj in all_objects if obj.startswith("run")]

if len(run_objects) < len(EXPECTED_CHAINS):
    print(f"ERROR: Expected {len(EXPECTED_CHAINS)} APBS runs, found {len(run_objects)}")
    print(f"  Found: {run_objects}")
    print("\nFAILURE: Complete APBS calculations for all chains")
    print("\nExpected order:")
    for i, chain in enumerate(EXPECTED_CHAINS, start=1):
        print(f"  Run {i}: chain_{chain}")
    sys.exit(1)

# Map runs to chains
apbs_mapping = {}

print("\nAPBS Objects Found:")
for i, chain in enumerate(EXPECTED_CHAINS, start=1):
    run_num = f"{i:02d}"  # 01, 02, etc.

    run_obj = f"run{run_num}"
    map_obj = f"apbs_map{run_num}"
    ramp_obj = f"apbs_ramp{run_num}"

    # Check if objects exist
    if run_obj not in all_objects:
        print(f"  WARNING: Missing {run_obj} for chain {chain}")
        continue

    if map_obj not in all_objects:
        print(f"  WARNING: Missing {map_obj} for chain {chain}")
        continue

    apbs_mapping[chain] = {
        'run': run_obj,
        'map': map_obj,
        'ramp': ramp_obj if ramp_obj in all_objects else None
    }

    print(f"  Chain {chain}: {run_obj}, {map_obj}, {ramp_obj}")

if len(apbs_mapping) < len(EXPECTED_CHAINS):
    print(f"\nERROR: Incomplete APBS data")
    print("\nFAILURE: Complete APBS calculations for all chains")
    print("\nExpected order:")
    for i, chain in enumerate(EXPECTED_CHAINS, start=1):
        print(f"  Run {i}: chain_{chain}")
    sys.exit(1)

print(f" All APBS maps found for chains: {list(apbs_mapping.keys())}")

print("\n All prerequisites satisfied")
print("=" * 60)

# ============================================================================
# Use the mapping in visualization
# ============================================================================

# Example usage:
# chain_D_map = apbs_mapping['D']['map']  # 'apbs_map01'
# chain_B_map = apbs_mapping['B']['map']  # 'apbs_map02'
