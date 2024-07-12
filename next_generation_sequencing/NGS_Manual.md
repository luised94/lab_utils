# Next-Generation Sequencing

## Overview
This document describes the scripts in the next-generation sequencing module of my lab_utils code repository. 
See README.md for specifications about the lab_utils repository.
Many of the scripts in this module depend on the slurm wrapper file found in the linux cluster module (../linux_cluster/000_sh_node_slurmWrapper.sh). Definitely take a look at the documentation for that script.
## Scripts
List of scripts in this module with brief descriptions.

## Workflow
High-level description of the module's workflow.

## Usage
General usage instructions for the module.

## Script Details
### 000_bioMicroCenterData
Most of the sequencing is carried out by the MIT BMC Core Facility. Therefore, the first steps of the project involve setting up sample data for submission. 
#### 001_processBMCSampleGridDataCSV.R
- **Purpose**: Read in the sample data from the BMC submission form, process to rename the columns, create the short name and sample_ID columns.
- **Input**: 
1. directory_to_process as argument
2. sample grid information as csv. directory_to_process will be used to find this input.
- **Output**: csv of sample grid information with new colnames, sample_ID and short_name columns.
- **Usage**: ./001_processBMCSampleGridDataCSV.R <directory_to_process>
- **Parameters**: Not applicable
- **Dependencies**: No dependencies
- **Notes**: Any important caveats or considerations

### 001_referenceGenomes
Download the reference genomes that are relevant to my labwork. These are required for all of the sequencing analysis. They are typically designed to be run once before starting the analysis.
#### 001_downloadReferenceGenomes.sh
- **Purpose**: Download reference genomes for read alignment.
- **Input**: Should work without inputs since the files are hard coded. 
- **Output**: Directories with genome as fna file inside a subdirectory.
- **Usage**: ./001_downloadReferenceGenomes.sh
- **Parameters**: Modify the accessions array to download other genomes. base_url may change if the API specifications change.
- **Dependencies**: NCBI datasets API, curl
- **Notes**: Some of the parameters are hard coded but they can be adjusted. The accessions were determined by searching in the database manually. It creates the REFGEN genome in place. The directory structure is weird. Finally, each genome comes with a README.md and the user should go through with it.

#### 002_reorganizeReferenceGenomesDirs.sh
- **Purpose**: Process the directories downloaded by 001_downloadReferenceGenomes.sh to remove the subdirectories, place everything in the top directory and rename the directory to a more descriptive. 
- **Input**: Output generated from 001_downloadReferenceGenomes.sh. Must be run inside REFGENS directory.
- **Output**: Renamed directory with all files in the top directory.
- **Usage**: (Inside the REFGENS directory) ./002_sh_node_reorganizeReferenceGenomesDirs.sh
- **Parameters**: Requires the REFGENS directory from 001_downloadReferenceGenomes.sh
- **Dependencies**: 001_downloadReferenceGenomes.sh
- **Notes**: Needs to run inside the REFGENS directory. 

#### 003_modifyS288cFastaHeaders.sh
- **Purpose**: Create a backup of the S288C genome then use awk to rename the fasta headers to UCSC format (chr <roman-numeral>)
- **Input**: Output from 002_reorganizeReferenceGenomesDirs. Only the Saccharomyces genome is required.
- **Output**: New fna file for S288C genome with modified headers.
- **Usage**: Not currently a script, run each line from the command line.
- **Parameters**: Requires previous script output.
- **Dependencies**: 002_reorganizeReferenceGenomesDirs.sh
- **Notes**: I never turned this into a script and just ran the lines on the command line cause I wasnt sure if it would work. 


#### 004_bt2buildRefGenomes.sh
- **Purpose**: Create the bowtie2 build index for all the reference genomes using slurm.
- **Input**: Reorganized genome directories with fna file. S288C should have been fasta modified.
- **Output**: bowtie2 build index inside the directories. <nameofgenome>_index
- **Usage**: Use through slurm wrapper. 
- **Parameters**: Genome directories
- **Dependencies**: 002_reorganizeReferenceGenomesDirs.sh, 003_modifyS288cFastaHeaders.sh, slurm_wrapper.sh, bowtie2/2.3.5.1
- **Notes**: Follows the usual pattern for most of my slurm scripts. Slurm header, create directories, echo information to log file for inspection, create array, index array using TASK_ID variable, perform command on selected element of array, echo more info to log file. 

### 002_controlData
#### 001_downloadEatonData.sh
- **Purpose**: Brief description of the script's function
- **Input**: Required input files/formats
- **Output**: Generated output files/formats
- **Usage**: Example command(s)
- **Parameters**: Description of important parameters
- **Dependencies**: Any required libraries or tools
- **Notes**: Any important caveats or considerations

### 003_fastqProcessing
#### [Script Name 1]
- **Purpose**: Brief description of the script's function
- **Input**: Required input files/formats
- **Output**: Generated output files/formats
- **Usage**: Example command(s)
- **Parameters**: Description of important parameters
- **Dependencies**: Any required libraries or tools
- **Notes**: Any important caveats or considerations

## Troubleshooting
Common issues and their solutions.

## References
Relevant papers, tools, or external resources.
