# {{Module}}

## Overview
This document describes the scripts in the next-generation sequencing module of my lab_utils code repository. 
See README.md for specifications about the lab_utils repository.

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

## Troubleshooting
Common issues and their solutions.

## References
Relevant papers, tools, or external resources.
