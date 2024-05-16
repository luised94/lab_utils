# lab_utils
Code used for laboratory analysis

Each directory represents a specific type of analysis, usually related to a technique or type of data. It also containts documentation and a directory of deprecated code.

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

# LOGGING
Verify output with vim ~/data/<dir>/logs/*_9004526_*_1.out

# STICKY_NOTES.md
Notes I take while developing the scripts. 
# COMPLETED
1. Basic code to filter and align via sbatch using SLURM_ARRAY_TASK_ID to index arrays.
2. Basic code to plot genomeTracks.
3. Wrapper for sbatch scripts.
4. Setting up R 4.2.0 on other systems.
5. Analyze Eaton data to plot it with my data
6. LOW- Create a script that consolidates logs files based on JOB_ID.

# WARNINGS
Should probably create tests for verifying these stats.

1. Verify fastq files for character values - and _ . They could break the naming. If this happens, just use awk to grab the last field.
2. For old fastq files, filtering by length is very different. Todays (2024) sequencings  are higher quality and longer. Right now, the conditional is for Eaton data, but needs to be based on average length of dataset.

# TODO

8. HIGH- Quality control for fastq and BAM!
5. HIGH- Make plots prettier, normalize data, incorporate origin features
1. MED- See Random R questions thread to create a list of packages that have been used and document use in scripts.
2. LOW- Output environment info and git info to NGS data directory. HIGH once done. 
4. MED- Verify indexing is approapriate for R and bash scripts. 
6. HIGH- Get peaks from bigwig files, then compare between samples
7. HIGH- Correlation matrix, heatmaps, feature correlations, 
9. MED- Create through documentation for all repositories
10. MED- Analyze other code in paperpile directory to redesing and incorporate into code repository
11. LOW- Add better formatting for echo commands. 
12. MED- Use renv to manage R packages for each module.
