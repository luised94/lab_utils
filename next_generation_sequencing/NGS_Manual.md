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

#### [Script Name 2]
...

## Troubleshooting
Common issues and their solutions.

## References
Relevant papers, tools, or external resources.
