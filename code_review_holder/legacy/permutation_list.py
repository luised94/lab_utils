# initialize lists
from os import sep
from unicodedata import name

#Use for loops to iterate through arrays that contain labels for samples and get name combinations of samples. \
#The order of labeling is (rescue-copy allele)_(suppressors)_(cell cycle)_(auxin)_(antibody)_(number)
#Create array with conditions of samples
rescue_label = ["none", "WT","RA"] #Says which rescue copy is in sample if applicable
supp_label = ["none", "1EK", "4PS", "T2GE"] #Says which suppressor mutation is in sample if applicable
cell_cycle_label= ["aF", "Noc"] #Stands for alpha factor and nocodazole
auxin_label=["yes", "no"] #Yes means sample was treated with auxin
antibody_label=["Myc", "V5", "UM174"] #MYC is for degraded Orc4, V5 is for rescue copy, UM174 IPs MCM

# create empty list to store the combinations
sample_names = []
storage = []
storage_1=[]
#Generate labels for Input samples
for label in rescue_label[0:2]:
    name_=label+"_"+supp_label[0]
    print(name_)
    storage.append(name_)

for label in supp_label:
    name_=rescue_label[2]+"_"+label
    print(name_)
    storage.append(name_)
print(storage)

for sample in storage:
    name_=sample+"_"+cell_cycle_label[0]
    print(name_)
    storage_1.append(name_)

storage =[]
for sample in storage_1:
    name_=sample+"_"+auxin_label[1]
    print(name_)
    storage.append(name_)

storage_1=[]
for sample in storage:
    name_=sample+"_"+"input"
    print(name_)
    storage_1.append(name_)

sample_names= sample_names + storage_1
print(sample_names)


#Generates labels for none-rescue sample
#Generate first sample then succescively add labels
none_sample = rescue_label[0]+"_"+supp_label[0]
storage = []
storage_1 =[]
for label in cell_cycle_label:
    name_=none_sample+"_"+label
    print(name_)
    storage.append(name_)

print(storage)

for label in auxin_label:
    for sample in storage:
        name_=sample+"_"+label
        print(name_)
        storage_1.append(name_)

storage = []

for label in antibody_label:
    for sample in storage_1:
        name_=sample+"_"+label
        print(name_)
        storage.append(name_)
len(storage)

sample_names= sample_names + storage
print(sample_names)

#Generate labels for WT rescue samples
WT_sample = rescue_label[1]+"_"+supp_label[0]

storage = []
storage_1 =[]

for label in cell_cycle_label:
    name_=WT_sample+"_"+label
    print(name_)
    storage.append(name_)

print(storage)

for sample in storage:
        name_=sample+"_"+auxin_label[0]
        print(name_)
        storage_1.append(name_)

storage = []

for label in antibody_label:
    for sample in storage_1:
        name_=sample+"_"+label
        print(name_)
        storage.append(name_)
len(storage)

sample_names= sample_names + storage
print(sample_names)
len(sample_names)
#Generate labels for 4R samples

storage = []
storage_1 =[]
for label in supp_label:
    name_=rescue_label[2]+"_"+label
    print(name_)
    storage.append(name_)

for label in cell_cycle_label:
    for sample in storage:
        name_=sample+"_"+label
        print(name_)
        storage_1.append(name_)
print(storage_1)

storage=[]
for sample in storage_1:
        name_=sample+"_"+auxin_label[0]
        print(name_)
        storage.append(name_)
print(storage)


storage_1 = []
for label in antibody_label:
    for sample in storage:
        name_=sample+"_"+label
        print(name_)
        storage_1.append(name_)
len(storage_1)

sample_names= sample_names + storage_1
print(sample_names)

len(sample_names)

print(*sample_names, sep="\n")


#Create file. Add write sample name while adding number to it. Not generalizable way to add number
f = open('CHIP_sample_names.txt', 'w')
count = 1
for names in sample_names:
    if count >= 10:
        f.write(names+"_"+"0"+str(count)+"\n")
    else: 
        f.write(names+"_"+"00"+str(count)+"\n")
    count = count + 1

f.close()
