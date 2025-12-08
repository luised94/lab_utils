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
    output_dir: str = "pymol_outputs"
) -> bool:
    """
    Save PyMOL image with overwrite protection.

    Args:
        filename (str): Output filename
        width (int): Image width in pixels
        height (int): Image height in pixels
        dpi (int): Resolution in dots per inch (default: 300)
        ray (bool): Whether to use ray tracing (default: True)
        overwrite (bool): Whether to overwrite existing files (default: False)
        output_dir (str): Output directory (default: "pymol_outputs")

    Returns:
        bool: True if image was saved, False if skipped
    """
    filepath = os.path.join(output_dir, filename)

    if os.path.exists(filepath) and not overwrite:
        print(f"  SKIP (already exists): {filename}")
        return False
    else:
        cmd.png(filepath, width=width, height=height, dpi=dpi, ray=ray)
        print(f"  SAVED: {filename}")
        return True

def prepare_structure():
    """
    Standard structure preparation applied to all sessions.
    Cleans structure, sets default display, handles ligands.
    Call this after loading your structure with cmd.fetch() or cmd.load()
    """
    from pymol import cmd

    print("\nPreparing structure...")

    # Remove clutter
    cmd.remove("hydrogen")
    cmd.remove("solvent")
    print("  Removed hydrogen and solvent")

    # Create standard selections
    try:
        cmd.select("lig", "organic")
        if cmd.count_atoms("lig") > 0:
            print(f"  Found ligand: {cmd.count_atoms('lig')} atoms")
    except:
        pass

    try:
        cmd.select("bb", "backbone")
    except:
        pass

    # Set default display
    cmd.hide("everything")
    cmd.show("cartoon", "polymer")
    print("  Set cartoon representation")

    # Show ligand if present
    if cmd.count_atoms("lig") > 0:
        cmd.show("sticks", "lig")
        cmd.color("black", "lig")
        print("  Showing ligand as sticks")

    # Color carbons by heteroatom colors
    cmd.util.cnc()
    print("  Applied carbon coloring scheme")

    # Orient around backbone
    if cmd.count_atoms("bb") > 0:
        cmd.orient("bb")
        print("  Oriented around backbone")

    cmd.deselect()
    print("  Preparation complete")
