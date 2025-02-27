
cmd.reinitialize()
import os

user_path = os.path.expanduser("~")

if user_path == "C:\\Users\\Luis":\
    cmd.load("C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Projects\\automate-the-boring-stuff\\reference_data\\5zr1.pdb")\
    cmd.load("C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Projects\\automate-the-boring-stuff\\reference_data\\4xgc.pdb")\
    cmd.run("C:\\Users\\Luis\\Dropbox (MIT)\\Lab\Projects\\automate-the-boring-stuff\\script\\niceify2.pml")\
elif user_path == "C:\\Users\\liusm":\
    cmd.load("C:\\Users\\liusm\\Dropbox (MIT)\\Lab\\Projects\\automate-the-boring-stuff\\reference_data\\5zr1.pdb")\
    cmd.load("C:\\Users\\liusm\\Dropbox (MIT)\\Lab\\Projects\\automate-the-boring-stuff\\reference_data\\4xgc.pdb")\
    cmd.run("C:\\Users\\liusm\\Dropbox (MIT)\\Lab\Projects\\automate-the-boring-stuff\\script\\niceify2.pml")


initial_object = '5zr1'
# Get a list of unique chain identifiers in the initial object
chain_list = cmd.get_chains(initial_object)
cmd.select('ars', f'(c. G or c. H) and 5zr1')
for chain in chain_list:\
  if chain != "G" and chain != "H":\
    cmd.select(f'{chain}nearDNA', f'ars expand 8 and 5zr1 and c. {chain}')

initial_object = '4xgc'
for chain in chain_list:\
  print(chain)\
  new_object_name = f"{initial_object}_{chain}"\
  cmd.create(new_object_name, f"{initial_object} and chain {chain}")
# Loop through the chain identifiers
# Create a new object name by combining the initial object name and the chain identifier
# Extract the chain into a new object
initial_object = '5zr1'
for chain in chain_list:\
  print(chain)\
  new_object_name = f"{initial_object}_{chain}"\
  cmd.create(new_object_name, f"{initial_object} and chain {chain}")

# Replace 'reference_object' with the name of the object you want to align to
reference_object = '4xgc'

# Define a list of dictionaries, each containing residue range, chain, and object information
alignment_info = [\
    {'res_range1': '8-223', 'chain1': 'D', 'object1': '4xgc', 'res_range2': '49-292', 'chain2': 'D', 'object2': '5zr1'},\
    {'res_range1': '550-798', 'chain1': 'A', 'object1': '4xgc', 'res_range2': '421-718', 'chain2': 'A', 'object2': '5zr1'},\
    {'res_range1': '124-599', 'chain1': 'E', 'object1': '4xgc', 'res_range2': '14-393', 'chain2': 'E', 'object2': '5zr1'},\
    {'res_range1': '18-379', 'chain1': 'C', 'object1': '4xgc', 'res_range2': '7-381', 'chain2': 'C', 'object2': '5zr1'},\
    {'res_range1': '283-612', 'chain1': 'B', 'object1': '4xgc', 'res_range2': '246-613', 'chain2': 'B', 'object2': '5zr1'}\
]
#{'res_range1': '514-591', 'chain1': 'B', 'object1': '4xgc', 'res_range2': '471-561', 'chain2': 'B', 'object2': '5zr1'}\

# Loop through the alignment information
# Create selections for the specific regions you want to superimpose
# Superimpose the current chain to the reference object '4xgc'
for info in alignment_info:\
  print(info['res_range1'])\
  res_range1 = info['res_range1']\
  chain1 = info['chain1']\
  object1 = info['object1']\
  res_range2 = info['res_range2']\
  chain2 = info['chain2']\
  object2 = info['object2']\
  cmd.select('region1', f'{object1} and c. {chain1} and resi {res_range1}')\
  cmd.select('region2', f'{object2}_{chain2} and resi {res_range2}')\
  cmd.super('region2', 'region1')

residue_list = [("D", 225), ("C", 481), ("A", 495), ("E", 104), ("F", 305)]  # Add more residues as needed

# Loop through the list
# Create the selection for a specific residue and chain
# Expand the selection
# Orient the view around the expanded selection
# Wait for user input to proceed to the next residue
# import time
for chain, residue in residue_list:\
    residue_str = str(residue)\
    cmd.select(f"sofr_{chain}_{residue_str}", f"5zr1_{chain} and resi {residue_str}")\
    cmd.select(f"near_{chain}_{residue_str}", f"sofr_{chain}_{residue_str} expand 8")\
    
cmd.select("electro", "5zr1_D or 5zr1_B")
cmd.select("electro2", '4xgc_B or 4xgc_D')

alignment_info = [\
    {'res_range1': '514-616', 'chain1': 'B', 'object1': '4xgc', 'res_range2': '505-520', 'chain2': 'B', 'object2': '5zr1'}\\
]

for info in alignment_info:\
  res_range1 = info['res_range1']\
  chain1 = info['chain1']\
  object1 = info['object1']\
  res_range2 = info['res_range2']\
  chain2 = info['chain2']\
  object2 = info['object2']\
  cmd.select('region1', f'{object1} and c. {chain1} and resi {res_range1}')\
  cmd.select('region2', f'{object2}_{chain2} and resi {res_range2}')\
  cmd.super('region2', 'region1')

# from pmg_tk.startup.apbs_gui.creating import pdb2pqr_cli
# from pmg_tk.startup.apbs_gui.electrostatics import map_new_apbs

# cmd.fetch("1ubq")
# 
# pdb2pqr_cli("prepared01", "1ubq", options=["--ff", "amber"])
# map_new_apbs("apbs_map01", "prepared01")
# 
# cmd.ramp_new("apbs_ramp01", "apbs_map01", [-5, 0, 5])
# cmd.set("surface_ramp_above_mode", 1, "prepared01")
# cmd.set("surface_color", "apbs_ramp01", "prepared01")
# cmd.show("surface", "prepared01")
#{'res_range1': '514-591', 'chain1': 'B', 'object1': '4xgc', 'res_range2': '471-561', 'chain2': 'B', 'object2': '5zr1'}\

# Loop through the alignment information
# Create selections for the specific regions you want to superimpose
# Superimpose the current chain to the reference object '4xgc'
for info in alignment_info:\
  print(info['res_range1'])\
  res_range1 = info['res_range1']\
  chain1 = info['chain1']\
  object1 = info['object1']\
  res_range2 = info['res_range2']\
  chain2 = info['chain2']\
  object2 = info['object2']\
  cmd.select('region1', f'{object2} and c. {chain2} and resi {res_range2}')\
  cmd.select('region2', f'{object1}_{chain1} and resi {res_range1}')\
  cmd.super('region2', 'region1')

chain_list = cmd.get_chains('4xgc')
for chain in chain_list:\
  if chain != "F" and chain != "G":\
    cmd.select(f'{chain}nearDNA', f'ars expand 6 and 4xgc_{chain}')

cmd.select('D223to227', 'resi 223-227 and 5zr1_D')
cmd.select('near225', 'D223to227 expand 8 and 5zr1_B')
# Get a list of unique chain identifiers in the initial object
# chain_list = cmd.get_chains(initial_object)
# residue_range 
# Loop through the chain identifiers
# for chain in chain_list:\
#     # Create a new object name by combining the initial object name and the chain identifier
#     chain_object_name = f"{initial_object}_{chain}"
#     
#     # Select the specific regions you want to superimpose for the current chain
#     cmd.select('region1', f'{reference_object} and c. {chain} and resi 8-223')
#     cmd.select('region2', f'{chain_object_name} and c. {chain} and resi 49-292')
#     
#     # Superimpose the current chain to the reference object
#     cmd.super('region2', 'region1')
#     
#     # Delete the temporary selections
#     cmd.delete('region1')
#     cmd.delete('region2')
