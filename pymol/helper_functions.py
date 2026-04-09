"""
PyMOL Helper Functions
User-defined utilities for molecular visualization and image export.
Place the function into the working directory of the pymol.
Add help_ prefix to new functions to mark as part of script.
To import in scripts:
# Import helper functions
import helper_functions # Without the file extesion.

"""

import os
from pymol import cmd
from typing import Optional


def save_png(
    filename: str,
    width: int,
    height: int,
    dpi: int = 300,
    ray: bool = True,
    overwrite: bool = False,
    output_dir: str = "pymol_outputs",
) -> bool:
    """
    Save PyMOL image with overwrite protection and automatic directory creation.

    Args:
        filename (str): Output filename (e.g., "my_image.png")
        width (int): Image width in pixels
        height (int): Image height in pixels
        dpi (int): Resolution in dots per inch (default: 300)
        ray (bool): Whether to use ray tracing (default: True)
        overwrite (bool): Whether to overwrite existing files (default: False)
        output_dir (str): Output directory (default: "pymol_outputs")

    Returns:
        bool: True if image was saved, False if skipped (file exists, overwrite=False)

    Example:
        save_png("view1.png", width=1200, height=900, ray=True)
    """
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"  Created directory: {output_dir}")

    filepath = os.path.join(output_dir, filename)

    # Check if file exists and overwrite setting
    if os.path.exists(filepath) and not overwrite:
        print(f"  SKIP (already exists): {filename}")
        return False

    # Save the image
    try:
        cmd.png(filepath, width=width, height=height, dpi=dpi, ray=ray)
        print(f"  SAVED: {filename}")
        return True
    except Exception as e:
        print(f"  ERROR saving {filename}: {e}")
        return False


def prepare_structure() -> bool:
    """
    Standard structure preparation applied to all visualization sessions.

    **IMPORTANT: Call this AFTER loading your structure** with cmd.fetch() or cmd.load()

    This function:
    - Removes hydrogen atoms and solvent molecules
    - Creates standard selections: 'lig' (ligands), 'bb' (backbone)
    - Sets cartoon representation for proteins
    - Shows ligands as sticks (if present)
    - Applies carbon coloring scheme
    - Orients view around backbone

    Returns:
        bool: True if preparation succeeded, False if no structure loaded

    Raises:
        RuntimeError: If no structure is loaded in PyMOL session

    Example:
        cmd.fetch("1nh2")
        prepare_structure()  # Now ready for custom coloring and selections
    """
    print("\nPreparing structure...")

    # Check if any structure is loaded
    objects = cmd.get_object_list()
    if len(objects) == 0:
        error_msg = (
            "ERROR: No structure loaded in PyMOL session!\n"
            "  -> Call cmd.fetch(pdb_code) or cmd.load(filename) BEFORE prepare_structure()\n"
            "  -> Example: cmd.fetch('1nh2')"
        )
        print(error_msg)
        raise RuntimeError("No structure loaded. Cannot prepare.")

    print(f"  Found {len(objects)} object(s): {', '.join(objects)}")

    # Remove clutter
    cmd.remove("hydrogen")
    cmd.remove("solvent")
    print("  Removed hydrogen and solvent")

    # Create standard selections
    try:
        cmd.select("lig", "organic")
        if cmd.count_atoms("lig") > 0:
            print(f"  Created 'lig' selection: {cmd.count_atoms('lig')} atoms")
        else:
            print("  No ligands found (selection 'lig' is empty)")
    except Exception as e:
        print(f"  Could not create ligand selection: {e}")

    try:
        cmd.select("bb", "backbone")
        bb_count = cmd.count_atoms("bb")
        if bb_count > 0:
            print(f"  Created 'bb' selection: {bb_count} atoms")
        else:
            print("  WARNING: No backbone found - is this a protein structure?")
    except Exception as e:
        print(f"  Could not create backbone selection: {e}")

    # Set default display
    cmd.hide("everything")
    cmd.show("cartoon", "polymer")
    print("  Set cartoon representation")

    # Show ligand if present
    if cmd.count_atoms("lig") > 0:
        cmd.show("sticks", "lig")
        cmd.color("black", "lig")
        print("  Showing ligand as sticks (black)")

    # Color carbons by heteroatom colors
    cmd.util.cnc()
    print("  Applied carbon coloring scheme (C=element-color, heteroatoms colored)")

    # Orient around backbone
    if cmd.count_atoms("bb") > 0:
        cmd.orient("bb")
        print("  Oriented view around backbone")
    else:
        # Fallback to orienting around everything
        cmd.orient()
        print("  Oriented view around entire structure (no backbone found)")

    # Clear selections
    cmd.deselect()

    print("   Preparation complete")
    return True


def reset_scene_to_base_cartoon(
    chain_colors: dict[str, str], cartoon_selection: str = "polymer"
) -> None:
    """
    Reset PyMOL scene to a known base state before spec-specific rendering.

    Hides all representations, restores chain colors from the provided mapping,
    then shows cartoon for the given selection. Called at the start of each image
    spec action list to prevent color or representation state from leaking between
    render passes.

    Args:
        chain_colors (dict[str, str]): Mapping of chain identifier to PyMOL color
            name (e.g., {"A": "yellow", "B": "marine"}). Applied to all chains.
        cartoon_selection (str): PyMOL selection string for cartoon representation
            (default: "polymer").

    Returns:
        None

    Example:
        reset_scene_to_base_cartoon(CHAIN_COLORS)
        reset_scene_to_base_cartoon(CHAIN_COLORS, cartoon_selection="chain A")
    """
    cmd.hide("everything")
    for chain, color in chain_colors.items():
        cmd.color(color, f"chain {chain}")
    cmd.show("cartoon", cartoon_selection)
