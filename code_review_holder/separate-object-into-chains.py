

#run C:/Users/liusm/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff/script/separate-object-into-chains.py '5zr1'
# import sys
# # Replace 'your_initial_object' with the name of your initial object
cmd.set("stick_radius" , .3)

initial_object = '5zr1'

# Get a list of unique chain identifiers in the initial object
chain_list = cmd.get_chains(initial_object)

# Loop through the chain identifiers
# Create a new object name by combining the initial object name and the chain identifier
# Extract the chain into a new object
for chain in chain_list: \
  print(chain) \
  new_object_name = f"{initial_object}_{chain}" \
  cmd.create(new_object_name, f"{initial_object} and chain {chain}")
