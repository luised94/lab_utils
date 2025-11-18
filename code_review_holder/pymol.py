#
cmd.fetch('5zr1')
cmd.seq_view(1)
# Select chains G and H
cmd.select('ARS_DNA', '5zr1 and (chain G or chain H)')
# Hide the selected chains
cmd.hide('everything', 'ARS_DNA')
cmd.select('ORC', '5zr1 and (chain A or chain B or chain C or chain D or chain E or chain F)')

# Define the chain-residue range pairs
chain_residue_ranges = {
    'A': (10, 20),
    'B': (30, 40),
    'C': (50, 60)
}

# Select specified residue ranges for each chain
for chain, (start_res, end_res) in chain_residue_ranges.items():
    selection_name = f'chain_{chain}_res_{start_res}_{end_res}'
    selection_string = f'my_structure and chain {chain} and res {start_res}-{end_res}'
    cmd.select(selection_name, selection_string)
