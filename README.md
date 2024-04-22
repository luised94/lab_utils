# lab_utils
Code used for laboratory analysis

# NAMING CONVETION

Directories - area_of_analysis/NUM_descriptiveName/script.ext

area_of_analysis: snake_case, the biological area of inquiry relative to the code inside, usually related to the technique or the type of data
NUM: Three Digit Integer, Number that serves as unique ID but is related to the order in which the scripts inside the directory are usually run (dependence between the directories)
descriptiveName: camelCase, describes as concise as possible the purpose of the script

Scripts - NUM_programminglanguage_placetorun_descriptiveName.ext
Scripts will be under most approapriate diretory according to its function, biological area, technique and type of data
NUM: Three Digit Integer, Number that serves as unique ID but is related to the order in which the scripts inside the directory are usually run (dependence between the directories)
programminglanguage: Name of programming language (e. g. R, bash, etc) This way they can be sorted by this as well. Should be consistent with extension.
placetorun: Whether to run on node or slurm. Node just means in interactive mode (srun --pty bash). Slurm means to run with sbatch. 
descriptiveName: camelCase, describes as concise as possible the purpose of the script

# STICKY_NOTES.md
Notes I take while developing the scripts. 

# TODO

1. MED- See Random R questions thread to create a list of packages that have been used and document use in scripts.
2. HIGH- Create genome track plotting code.
3. MED- Create rsync tool for plots of experiments
4. LOW- Output environment info and git info to NGS data directory. HIGH once done. 
5. HIGH- Analyze Eaton data to plot it with my data
6. MED- Verify indexing is approapriate for R and bash scripts. 
