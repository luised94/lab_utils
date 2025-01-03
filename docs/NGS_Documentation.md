#STATUS: REMOVE.
# NGS Analysis Pipeline Documentation

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Folder Structure](#folder-structure)
5. [Workflow Steps](#workflow-steps)
6. [Configuration](#configuration)
7. [Troubleshooting](#troubleshooting)
8. [References](#references)

## Introduction
Brief overview of the pipeline's purpose and capabilities.

## Installation
Instructions for setting up the environment and dependencies.

## Quick Start
Guide to running the quickstart script and basic usage.

## Folder Structure
Overview of the repository's organization.

## Workflow Steps
### Before you start
Ensure you downloaded reference genome and indexed it. Ensure you have downloaded feature files from other scientific articles or databases.

### Experiment analysis
1. Design and carry out an experiment with next-generation sequencing as its read out. (genetic screen or chromatin immunoprecipitation analysis)
2. Run 000_setupExperimentDir.R with proper categories, filter function, order, control experiments. (Directory name will be the day the BMC accepts request.)
3. Run 001_downloadBMCData.sh using the bmc server and directory. (See BMC email for instructions.)
4. Run 000_consolidateFastq.sh to consolidate fastq files into a single file. (Makes the files very large.)
5. Run 003_updateSampleGrid.R to ensure the sample table is properly loaded and the sample_ID column is added. (Requires consolidated fastq files to obtain sample IDs.)
6. In preparation: Run qualityControl scripts to ensure the fastq files are consistent and high quality. Use 002_filterFastq.sh for any preprocessing such as trimming and removing low quality reads.
7. Align fastq files to reference genome with slurmWrapper.sh and 003_alignFastq.sh.
8. In preparation: Run quality control on bam mapping files.
9. Generate bigwig files for all downstream analysis.
10. Plot genomic tracks.
11. Inpreparation: Determine peaks.
12. Inpreparation: Determine motifs.
13. Inpreparation: Perform factor analysis and comparisons with other experimental datasets and between samples.

### 1. Data Download
- Purpose: Description of this step
- Main Scripts:
- `script1.sh`: Brief description
- `script2.R`: Brief description
- Input: Required input files/formats
- Output: Generated output files/formats
- Usage: Example commands

### 2. FASTQ Processing
(Follow the same structure as above)

### 3. Alignment
(Follow the same structure as above)

### 4. Peak Calling
(Follow the same structure as above)

(Continue for each numbered folder)

## Configuration
Explanation of any configuration files or important parameters.

## Troubleshooting

1. Verify fastq files for character values - and _ . They could break the naming. If this happens, just use awk to grab the last field.
2. For old fastq files, filtering by length is very different. Todays (2024) sequencings  are higher quality and longer. Right now, the conditional is for Eaton data, but needs to be based on average length of dataset.

# TODO
### HIGH
- 004_processFeaturesToRdsAndBed.R: Normalize feature files and incorporate into genomic tracks. Need to add the code that outputs the reorganized files to Bed. Not sure if I should do RDS and something else.
### MED
### LOW

Each documentation for each subdirectory has a # TODO section 
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
13. Decide whether to use a script template or plain text to generate the categorical data for analysis. This would then be read into a script that would plot, compare peaks etc. 
14. Add an indication of file name generated by scripts with command run such that I can use grep to quickly identify the parameters of each file. 

# COMPLETED
1. Basic code to filter and align via sbatch using SLURM_ARRAY_TASK_ID to index arrays.
2. Basic code to plot genomeTracks.
3. Wrapper for sbatch scripts.
4. Setting up R 4.2.0 on other systems.
5. Analyze Eaton data to plot it with my data
6. Create a script that consolidates logs files based on JOB_ID.
7. Updated scripts in 000_generalSetup with a config file setup and output scp command.
## References
Relevant papers, tools, or external resources used.
