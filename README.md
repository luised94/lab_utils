# lab_utils
## Disclaimer
AI tools (Claude, Perplexity, etcetera) are used in the project via chat interfaces.
Usually, they help with brainstorming and alternative suggestions but some of the code may be directly taken from the suggestions.

## Overview
Code used for laboratory analysis of flow cytometry and next-generation sequencing data provided by my institution's core facility (Biomicro Center).

Most important scripts are found in core_scripts/ directory.

## Installation
```{bash}
git clone https://github.com/luised94/lab_utils.git
```

## Quick Start
See docs/bmc_chip_seq_analysis_instructions.qmd for sequencing analysis.

### 2025-08-27
README and documentation needs updating. Quick start is the easiest way to start for now.

### 2025-10-20
Through development of the collection of scripts, I have gone back and forth wiht different organization styles. Right now, I have settled on just sourcing files during an interactive session since that captures the environment after the script since I can then try code on the variables and do quick debugging sessions.

Not all of the scripts are updated to this style yet!

## Dependencies
Because I perform the next-generation sequencing analysis using my institution's cluster, I have to use the version of the tools that are installed there for the most part. 

For this reason, I use R 4.2.0 to perform the analyses.

The analysis are done locally (WSL2 on Windows) or in a linux computing cluster. The linux cluster is the condition that dictates what dependencies are used, especially for the next generation sequencing analysis.

1. R 4.2.0
2. bash 4.2
3. Command line utils
4. bowtie2/2.3.5.1
5. fastp/0.20.0
6. fastqc/0.11.5
7. deeptools/3.0.1
8. gatk
9. python/2.7
10. miniforge
11. macs2
12. picard
13. java

## Usage Examples
I assume scripts will be run while the current working directory is the repository root.

Most bash scripts can be used by running the script from the command line.
```{bash}
./script.sh <args>
Rscript script.R <args>
```

Most R scripts can be used by running the script from the command line.
```{bash}
R --no-save
# Start up message should show.
> source("core_scripts/<script_name>.R")
```

## LOGGING
Most scripts output some sort of log file (stdout and stderr) that can be inspected with a text editor. The log files can usually be verified with vim ~/data/<dir>/logs/*_9004526_*_1.out.

## TROUBLESHOOTING 
Each documentation section has a troubleshooting section that lets the user know about common errors that could be encountered, such as the scripts depending on the name of the files.

## STICKY_NOTES.md
Notes I take while developing the scripts. 

## TAGS 
IGNORE THIS SECTION.
I have a set of tags that I try to use to put marks on code for future reference. The form of the tags is <comment><TAG>. recursive (-r) grep can be used to find the tags.

TODO: Tasks that I have to complete for that particular code file. 
HOWTO: Designates different code snippets for reference when I want to see how to do a particular thing.
FIXME: Highlight areas that need fixing.
NOTE: Add important notes or explanations.
BUG: Mark known bugs or issues.
OPTIMIZE: Indicate areas that could be optimized for better performance.
REFACTOR: for code that needs refactoring
TEST: for testing purposes
