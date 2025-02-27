from pymol import cmd


def coloram(selection="all"):

    """
    AUTHOR
    Tamas Hegedus @ hegelab.org

    DESCRIPTION
    Colors structures according to b-factors with AlphaMissense colors

    USAGE
    coloram sele

    PARAMETERS

    sele (string)
    The name of the selection/object to color with AlphaMissense colors. Default: all
    """

    cmd.set_color("benign", [161, 176, 222])
    cmd.set_color("ambiguous", [230,230,230])
    cmd.set_color("patho0", [245, 125, 129])
    cmd.set_color("patho1", [165, 13, 18])

    cmd.color("benign", f"({selection}) and b < 0.34")
    cmd.color("ambiguous", f"({selection}) and (b = 0.34 or b > 0.34)")
    cmd.color("patho0", f"({selection}) and (b = 0.564 or (b > 0.564 and b < 0.78))")
    cmd.color("patho1", f"({selection}) and (b = 0.78 or b > 0.78)")

cmd.extend("coloram", coloram)
cmd.auto_arg[0]["coloram"] = [cmd.object_sc, "object", ""]
