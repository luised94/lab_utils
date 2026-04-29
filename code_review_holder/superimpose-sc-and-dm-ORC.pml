
#7JPO
#ORC-O1AAA: Human Origin Recognition Complex (ORC) with dynamic/unresolved ORC2 WH
#7JPP
#ORC-O2WH: Human Origin Recognition Complex (ORC) with dynamic/unresolved ORC1 AAA+ domain
#7JPQ
#ORC-O2-5: Human Origin Recognition Complex (ORC) with subunits 2,3,4,5
#7JPR
#ORC-OPEN: Human Origin Recognition Complex (ORC) in an open conformation
#7JPS
#ORC-DNA: Human Origin Recognition Complex (ORC) with DNA bound in the core
reinitialize
#Have to open structures manually
# Load the first structure
#fetch 4xgc, "4xgc", type=pdb

# Load the second structure
#fetch 5zr1, "4xgc", type=pdb

sele c. A and 5zr1 
create orc1, sele

# Select the specific regions you want to superimpose
select region1, 4xgc and c. D and resi 8-223
select region2, 5zr1 and c. D and resi 49-292

# Superimpose the selected regions
super region2, region1


hide everything, 5zr1

set_view (\
     0.946230829,    0.322090656,   -0.030079708,\
     0.196780011,   -0.646900117,   -0.736748159,\
    -0.256758302,    0.691214800,   -0.675497890,\
     0.000000000,    0.000000000, -136.858642578,\
   140.274871826,  125.757171631,  129.013458252,\
   107.900398254,  165.816894531,  -20.000000000 )

sele dnanear, dna expand 5 and ((c. B or c. A) and 4xgc)
orient dnanear
set_view (\
     0.851048648,   -0.183372393,    0.492031991,\
     0.382366151,   -0.425801396,   -0.820055783,\
     0.359882444,    0.886038125,   -0.292260438,\
    -0.000000000,    0.000000000, -186.364059448,\
   149.004013062,  129.379318237,  130.774337769,\
   146.930847168,  225.797271729,  -20.000000000 )
sele o4sofr, c. D and resi 225 and 5zr1
sele o4near, o4sofr expand 8 and 
orient near

set_view (\
    -0.331043869,    0.914366484,   -0.233120322,\
     0.203296944,   -0.172138184,   -0.963867724,\
    -0.921454132,   -0.366473675,   -0.128902197,\
     0.000000000,    0.000000000,  -60.750503540,\
   143.643005371,  127.049209595,  104.828681946,\
    47.896160126,   73.604843140,  -20.000000000 )
    
sele c. B and 4xgc    
sele o1nearo2, sele expand 5 and (c. A and 4xgc) 
  
# Superimpose the selected regions
super region1, region2

# Select the specific regions you want to superimpose
select region1, 4xgc and c. A and resi 550-798
select region2, 5zr1 and c. A and resi 421-718
super region1, region2
set_view (\
     0.519189894,    0.124325894,    0.845567644,\
     0.737576544,    0.434644938,   -0.516789436,\
    -0.431766838,    0.891974092,    0.133963212,\
    -0.000000000,    0.000000000,  -89.461975098,\
   141.611175537,  123.853889465,  105.002700806,\
    70.532501221,  108.391448975,  -20.000000000 )



sele o1sofr, c. A and resi 495 and 5zr1
sele o1near, o1sofr expand 8
orient o1near

# Select the specific regions you want to superimpose
select region1, 4xgc and c. E and resi 124-599
select region2, 5zr1 and c. E and resi 14-393

# Superimpose the selected regions
super region1, region2

sele o5sofr, c. E and resi 104 and 5zr1
sele o5near, o5sofr expand 8
orient o5near

# Select the specific regions you want to superimpose
select region1, 4xgc and c. C and resi 18-379
select region2, 5zr1 and c. C and resi 7-381

# Superimpose the selected regions
super region1, region2

sele o3sofr, c. C and resi 481 and 5zr1
sele o3near, o3sofr expand 8
orient o3near


# Select the specific regions you want to superimpose
select region1, 4xgc and c. B and resi 283-612
select region2, 5zr1 and c. B and resi 246-613

# Superimpose the selected regions
super region1, region2
super c. B and 5zr1, c. B and 4xgc

sele O2C1, 5zr1 and resi 471-561 and c. B
sele O2C2, 4xgc and resi 514-591 and c. B

super O2C1, O2C2

set_view (\
     0.765575647,   -0.600903273,    0.229800045,\
    -0.541925609,   -0.794843674,   -0.273018122,\
     0.346710503,    0.084480882,   -0.934160173,\
     0.000000000,    0.000000000,  -89.584808350,\
   141.684921265,  124.003890991,  105.036018372,\
    70.629348755,  108.540267944,  -20.000000000 )

super (c. B and 5zr1), (c. B and 4xgc)


reinitialize
#reload structures
sele dna, c. G or c. H


sele o1dnanear, dna expand 6 and 5zr1 and c. A 
sele o2dnanear, dna expand 6 and 5zr1 and c. B
sele o3dnanear, dna expand 6 and 5zr1 and c. C
sele o4dnanear, dna expand 6 and 5zr1 and c. D
sele o5dnanear, dna expand 6 and 5zr1 and c. E

#sele dnanear, dna expand 6 and 5zr1 and (c. A or c. B or c. C or c. D or c. E)
# Define a list of residues with their corresponding chains
residue_list = [("D", 225), ("C", 481), ("A", 495), ("E", 104), ("F", 305)]  # Add more residues as needed

# Loop through the list
for chain, residue in residue_list:
    residue_str = str(residue)
    
    # Create the selection for a specific residue and chain
    cmd.select(f"sofr_{chain}_{residue_str}", f"c. {chain} and resi {residue_str} and 5zr1")
    
    # Expand the selection
    cmd.select(f"near_{chain}_{residue_str}", f"sofr_{chain}_{residue_str} expand 8")
    
    # Orient the view around the expanded selection
    cmd.orient(f"near_{chain}_{residue_str}")
    
    # Wait for user input to proceed to the next residue
    input(f"Press Enter to select the next residue, or Ctrl+C to exit...")
    
# Clear the selections for the next iteration
    cmd.delete(f"o4sofr_{chain}_{residue_str}")
    cmd.delete(f"o4near_{chain}_{residue_str}")    
# Apply the superimposition to the entire second structure
transform "structure2", position

# Save the superimposed structure
save "superimposed.pdb", "structure2"

# Delete the temporary selections
delete region1
delete region2
