import os

# Define the path based on the current user's home directory
user_path = os.path.expanduser("~")
path = ""
cmd.reinitialize()
if user_path == "C:\\Users\\Luis":\
    path = "C:/Users/Luis/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff/"\
elif user_path == "C:\\Users\\liusm":\
    path = "C:/Users/liusm/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff/"

# Load the PDB file
cmd.load(os.path.join(path, "reference_data/7mca.pdb"))

# Run the PyMOL script
cmd.run(os.path.join(path, "script/niceify2.pml"))

# Set PyMOL properties
cmd.set("cartoon_side_chain_helper", "on")\
cmd.color("orange", "c. I and 7MCA")\
cmd.set("stick_radius", 0.3)

# Define selections
cmd.select("suppressors", "(c. A and resi 495) or (c. C and resi 481) or (c. D and resi 225) or (c. E and resi 104) or (c. F and resi 305)")\
cmd.select("near", "suppressors expand 8")\
cmd.show("sticks", "suppressors or (c. A and resi 485)")\
cmd.select("intheway", "(c. D and resi 134-147) or (c. D and resi 171-183) or (c. B and resi 232-246) or (c. B and resi 456-466)")\
cmd.hide("everything", "intheway")

# Define distances
cmd.distance("4PSD", "/7MCA/D/D/PRO`225/CB", "/7MCA/H/H/DC`61/P")\
cmd.distance("5EKK", "/7MCA/E/E/GLU`104/CD", "/7MCA/A/A/LYS`362/CA")\
cmd.distance("5EKD", "/7MCA/E/E/GLU`104/CD", "/7MCA/G/G/DT`21/P")\
cmd.distance("1EKA", "/7MCA/A/A/GLU`495/CD", "/lig/K/D/AGS`2001/C5'")\
cmd.distance("6EKD", "/7MCA/F/F/GLU`305/CD", "/7MCA/G/G/DA`46/P")

# Set rendering properties
cmd.hide("labels")\
cmd.set("dash_gap", 0.5)\
cmd.set("dash_radius", 0.1)\
cmd.set("ray_trace_fog", 0)\
cmd.set("ray_shadows", 0)\
cmd.unset("depth_cue")\
cmd.bg_color("white")\
cmd.set("antialias", 2)\
cmd.set("hash_max", 300)

# Define views and render images
views = [\
    {"name": "O1_supp.png",\ 
    "view": (0.426406741,   -0.417760849, -0.802281439, 0.807456851,\
    0.575538695, 0.129467249, 0.407653809, -0.703014195,\
    0.582741916, -0.002135545, 0.000229927, -87.124885559,\
    129.596054077, 138.158264160, 161.141448975, 66.774276733,\
    107.287414551, -20.000000000)},\
    {"name": "O3_supp.png",\ 
    "view": (0.276172251,   -0.958077013,    0.076264359,\
     0.926068425,    0.286498696,    0.245602220,\
    -0.257153094,    0.002797240,    0.966364503,\
     0.000000000,    0.000000000,  -49.708568573,\
   115.162040710,  118.990631104,   76.779479980,\
    39.190612793,   60.226524353,  -20.000000000 )},\
    {"name": "O4_supp.png",\
    "view": (\
    -0.611087084,    0.467109442,    0.639039159,\
     0.651989937,   -0.160774395,    0.740987599,\
     0.448863804,    0.869451940,   -0.206302688,\
    -0.000000000,    0.000000000,  -84.202651978,\
   121.166198730,  117.056610107,  141.704757690,\
    66.386016846,  102.019287109,  -20.000000000 )},\
    {"name": "O5_supp.png",\ 
    "view": (\
    -0.752292931,   -0.426453650,    0.502184510,\
     0.615914285,   -0.725833297,    0.306290090,\
     0.233883113,    0.539721191,    0.808698237,\
     0.000000000,    0.000000000,  -72.639938354,\
   101.370376587,  118.207054138,  130.884613037,\
    57.269882202,   88.009994507,  -20.000000000 )},\
    {"name": "O6_supp.png",\ 
    "view": (\
     0.765645385,    0.437111735,   -0.471919209,\
    -0.340094328,    0.897789836,    0.279806972,\
     0.545994639,   -0.053733055,    0.836059928,\
    -0.001210315,   -0.003005750,  -85.177246094,\
   152.192550659,  111.454238892,   63.569774628,\
    67.185455322,  103.247833252,  -20.000000000)}]

#Take picture of suppressors
for view_info in views:\
    cmd.set_view(view_info["view"])\
    cmd.ray(2400, 2400)\
    cmd.png(os.path.join(path, "output", view_info["name"]), 0, 0, -1, ray=0)

#Overall ORC/Cdc6 structure
view = [\
    -0.494458616, -0.225710392, 0.839386284,\
    -0.869120657, 0.141693920, -0.473868966,\
    -0.011976970, -0.963827670, -0.266235441,\
    0.000000000, 0.000000000, -509.880279541,\
    117.735023499, 114.526748657, 114.030647278,\
    401.993499756, 617.767089844, -20.000000000\
]

cmd.set_view(view)

# Define the number of times to run the commands
num_iterations = 4
degrees = 0
# Define the PNG file name with a suffix based on degrees and iteration
# Generate the PNG image
# Increment degrees for the next iteration
# Turn the object by 90 degrees for the next iteration
for iteration in range(num_iterations):\
    cmd.ray(2400, 2400)\
    png_file = os.path.join(path, "output", f"ORC_Cdc6_complete_{degrees}.png")\
    cmd.png(png_file, 0, 0, -1, ray=0)\
    degrees += 90\
    cmd.turn("y", 90)

# Execute PyMOL commands
cmd.set("sphere_scale", 2)
cmd.show("spheres", "suppressors")
cmd.color("red", "suppressors")

# Define the number of times to run the commands
num_iterations = 4
degrees = 0
# Define the PNG file name with a suffix based on degrees and iteration
# Generate the PNG image
# Increment degrees for the next iteration
# Turn the object by 90 degrees for the next iteration
for iteration in range(num_iterations):\
    cmd.ray(2400, 2400)\
    png_file = os.path.join(path, "output", f"ORC_Cdc6_complete_{degrees}_red.png")\
    cmd.png(png_file, 0, 0, -1, ray=0)\
    degrees += 90\
    cmd.turn("y", 90)

# Additional rendering with transparency
cmd.set("cartoon_transparency", 0.4)\

# Define the number of times to run the commands
num_iterations = 4
degrees = 0
# Define the PNG file name with a suffix based on degrees and iteration
# Generate the PNG image
# Increment degrees for the next iteration
# Turn the object by 90 degrees for the next iteration
for iteration in range(num_iterations):\
    cmd.ray(2400, 2400)\
    png_file = os.path.join(path, "output", f"ORC_Cdc6_complete_{degrees}_trans40_red.png")\
    cmd.png(png_file, 0, 0, -1, ray=0)\
    degrees += 90\
    cmd.turn("y", 90)
