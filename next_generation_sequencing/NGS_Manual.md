# Next-Generation Sequencing

## Overview
This document describes the scripts in the next-generation sequencing module of my lab_utils code repository. 
See README.md for specifications about the lab_utils repository.
Many of the scripts in this module depend on the slurm wrapper file found in the linux cluster module (../linux_cluster/000_sh_node_slurmWrapper.sh). Definitely take a look at the documentation for that script.

## Scripts
List of scripts in this module with brief descriptions.

## Workflow
This assumes you are using a similar workflow to the used in the lab_utils repository.
If this is the first time starting a project, and you do not have any of the reference files or controls, then start with generalSetup and 002_controlData.
If you have reference, then create the experiment directory, and download the data from the BMC. (See generalSetup)
To do a quick start, you can run 003_alignFastq along with the slurm wrapper script to generate BAM files. Afterwards, run 003_generateCoverageFiles.sh with the slurm wrapper to generate bigwig files that can be used to plot tracks. This can give you a clear idea of the result of a CHIP-seq experiment.

## Usage
General usage instructions for the module.

## Script Details

### generalSetup
Create the directories and download the data from the BMC.

#### 000_directoryCreation.sh
- **Purpose**: Read in the sample data from the BMC submission form, process to rename the columns, create the short name and sample_ID columns.
- **Input**: 
1. directory_to_process as argument
2. sample grid information as csv. directory_to_process will be used to find this input.
- **Output**: csv of sample grid information with new colnames, sample_ID and short_name columns.
- **Usage**: ./001_processBMCSampleGridDataCSV.R <directory_to_process>
- **Parameters**: Not applicable
- **Dependencies**: No dependencies
- **Notes**: Any important caveats or considerations

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
Scripts to download control data that will be used to compare to my samples.

#### 001_downloadEatonData.sh
- **Purpose**: Download data from Eaton 2010 paper. 
- **Input**: directory to download files to.
- **Output**: directory with fastq files for ORC CHIP samples arrested in G2.
- **Usage**: $./001_downloadEatonData.sh EatonBel
- **Parameters**: Directory, mapfile lines denote the files to concatenate and accession numbers to download.
- **Dependencies**: curl, sra ftp api 
- **Notes**: Hardcoded filenames and accessions to download. Could convert into a more generic file to download and name fastq files with replicates

#### 002_downloadFeatureData.sh
- **Purpose**: Download data from Eaton 2010 paper. 
- **Input**: directory to download files to.
- **Output**: directory with fastq files for ORC CHIP samples arrested in G2.
- **Usage**: $./001_downloadEatonData.sh EatonBel
- **Parameters**: Directory, mapfile lines denote the files to concatenate and accession numbers to download.
- **Dependencies**: curl, sra ftp api 
- **Notes**: Hardcoded filenames and accessions to download. Could convert into a more generic file to download and name fastq files with replicates

#### 003_downloadHawkinsTimingData.R
- **Purpose**: Download hawkins timing data for categorical analysis and track plotting.
- **Input**: No input required.
- **Output**: hawkins timing data as an xlsx file. 
- **Usage**: $Rscript 003_downloadHawkinsTimingData.R
- **Parameters**: Hawkins Timing URL
- **Dependencies**: curl package, hawkins timing url
- **Notes**: URL is hardcoded into the script. For some reason, downloading from the url did not work when using curl from the command line, which is why I had to use the R curl package. Still need to process the files into a similar format to use them in categorical and track plotting analysis.

### 003_fastqProcessing
Scripts to perform processing and align fastq files.

#### 000_consolidateFastq.sh
- **Purpose**: Concatenate files if BMC provides them as two separate fastq files.
- **Input**: Fastq files from a BMC sequencing run and directory to process.
- **Output**: Fastq files in fastq directory.
- **Usage**: $./003_fastqProcessing/000_sh_node_consolidateFastq.sh EatonBel
- **Parameters**: Fastq files, UNIQUE_ID parsing
- **Dependencies**: 
- **Notes**: Processing to get unique IDs dependes on cut which is very vulnerable. Need to test new method for concatenation.

#### 002_filterFastq.sh
- **Purpose**: Filter fastq files using fastp for length and quality.
- **Input**: Directory with fastq files.
- **Output**: Fastq files in fastq directory.
- **Usage**: Use via slurm wrapper. 
- **Parameters**: Fastq files, fastp conditions.
- **Dependencies**: fastp/0.20.0
- **Notes**: Filtering depends on whether file is from Eaton paper. Need to see if I can make it more dependent on the length distribution or average distribution.

#### 003_alignFastq.sh
- **Purpose**: Align the filtered fastq files to all genomes in the REFGENS directory 
- **Input**: Directory with processed fastq files.
- **Output**: Index and sorted bam files for all fastq and genome pairs.
- **Usage**: Use via slurm wrapper. 
- **Parameters**: Must get number of tasks correct.
- **Dependencies**: bowtie2/2.3.5.1, samtools/1.10
- **Notes**: Filtering depends on whether file is from Eaton paper. Need to see if I can make it more dependent on the length distribution or average distribution.

#### 004_qualityControlFastq.sh
- **Purpose**: Use fastqc to get quality control statistics.
- **Input**: Fastq files, pre and post processing
- **Output**: Fastqc directories with the html file that can be opened with browser.
- **Usage**: Use via slurm wrapper with number of tasks.
- **Parameters**: Must get number of tasks correct.
- **Dependencies**: fastqc/0.11.5, samtools/1.10
- **Notes**: Fastqc files are not very useful since they are text. You can grab sections using tools such as grep and awk but they must be visualized to be useful.

#### 005_unzipFastqc.sh
- **Purpose**: Unzip the produced files from fastqc command.
- **Input**: Directory to processs with fastqc directories.
- **Output**: Directory with fastqc files unzipped which can be parsed and used for other analysis.
- **Usage**: Use via slurm wrapper with number of tasks.
- **Parameters**: Must get number of tasks correct.
- **Dependencies**: unzip
- **Notes**: Have to run from the command line since the find and unzip command doesnt work as part of the script.

#### 006_parseFastqc.R
- **Purpose**: Parse the unzipped directory of the fastqc file into distinct sections that can be accessed.
- **Input**: Directory to processs with fastqc directories with unzipped directories.
- **Output**: A collection of tab-delimited files that represents all the modules of the fastqc command. 
- **Usage**: ./006_parseFastqc.R <dir>
- **Parameters**: Directories with unzipped files. 
- **Dependencies**: Output from previous scripts, particularly 005_unzipFastqc.sh
- **Notes**: I wrote this custom script to have more control over plotting specific plots. Considering substituting these for official R packages but will likely commit to this setup since I will understand how to use the data instead of having to figure out other packages that may or may not plot the way I want.

### 004_bamProcessing
Scripts to handle and process bam files.

#### 001_qualityControlBam.sh
- **Purpose**: Use samtools to determine reads mapped and other metrics for the bam quality control.
- **Input**: Directory, requires bam files.
- **Output**: Tab-delimited text files with the metrics.
- **Usage**: 001_qualityControlBam.sh <dir>
- **Parameters**: Directories with required bam files.
- **Dependencies**: samtools/1.10
- **Notes**: Very basic statistics for the BAM files. Need to find out some other methods such as with deeptools.

#### 002_checkQcBam.R
- **Purpose**: Use R to load the quality control files from 001_qualityControlBam.sh and generate the report plots.
- **Input**: Directory, requires files generated by qualityControlBam.sh
- **Output**: Plot images for reports and documentation.
- **Usage**: 002_checkQcBam.R <dir>
- **Parameters**: Directories with required bam files.
- **Dependencies**: samtools/1.10
- **Notes**: Still needs a lot of work.

#### 003_generateCoverageFiles.sh
- **Purpose**: Generate the bigwig files that be used for plotting tracks and calling peaks.
- **Input**: Directory, requires bam files.
- **Output**: Output bigwig (.bw) files.
- **Usage**: 003_generateCoverageFiles.sh <dir>, use via slurm wrapper.
- **Parameters**: Directories with required bam files.
- **Dependencies**: samtools/1.10
- **Notes**:

### 005_genomeTracks
Scripts to create track plots for genomic data.

#### 001_createGenomeTrack.R
- **Purpose**: Plot genome tracks as svg.
- **Input**: Directory, requires bigwig files.
- **Output**: SVG files
- **Usage**: 001_qualityControlBam.sh 
- **Parameters**: Parameters are hard coded for now.
- **Dependencies**: QuasR, GenomicAlignments, Gviz, rtracklayer, ShortRead, tidyverse
- **Notes**: Need to figure out how to polish the plots and how to make this more modular.

## Troubleshooting
Common issues and their solutions.

## References
Relevant papers, tools, or external resources.
