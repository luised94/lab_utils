# lab_utils

## Overview
Code used for laboratory analysis

Each directory represents a specific type of analysis, usually related to a technique or type of data. It also containts documentation and a directory of deprecated code.

## Configuration
Because I perform the next-generation sequencing analysis using my institution's cluster, I have to use the version of the tools that are installed there for the most part. 
For this reason, I use R 4.2.0 to perform the analyses.
See 001_setupR/000_installingR4.2.0 to see what to run to install the appropriate version of R. It isnt in the style of a script yet. Need to update this.
## Quick Start

## Documentation 
Each folder has a documentation file that goes through the workflow and provides a general description of each script.
There is are also sections for:
- TODO: Lets me know what I have to implement. 
- TROUBLESHOOTING: Lets the user what to watch out for such as how a files name could affect its processing. 

## NAMING CONVETION
The folders and scripts follow a naming convention. Each one is described in this section.

### Directories 
- area_of_analysis/NUM_descriptiveName/script.ext
area_of_analysis: snake_case, the biological area of inquiry relative to the code inside, usually related to the technique or the type of data
NUM: Three Digit Integer, Number that serves as unique ID but is related to the order in which the scripts inside the directory are usually run (dependence between the directories)
descriptiveName: camelCase, describes as concise as possible the purpose of the scripts inside the directory.

### Scripts 
- NUM_descriptiveName.ext
Scripts will be under most approapriate diretory according to its function, biological area, technique and type of data. 

Scripts - NUM_programminglanguage_placetorun_descriptiveName.ext

NUM: Three Digit Integer, Number that serves as unique ID but is related to the order in which the scripts inside the directory are usually run (dependence between the directories)
descriptiveName: camelCase, describes as concise as possible the purpose of the script

## Dependencies
The analysis are done locally or in a linux computing cluster. The linux cluster is the condition that dictates what dependencies are used, especially for the next generation sequencing analysis.

1. R 4.2.0
2. Command line utils
3. bowtie2
4. fastp 
5. fastqc 
6. deeptools ( cluster version requires python 2.7)
7. gatk
## Installation
```{bash}
git clone https://github.com/luised94/lab_utils.git
```

## Usage Examples
Most scripts can be used by running the script from the command line.
```{bash}
./script.sh <args>
Rscript script.R <args>
```
However I do have a script that serves as a wrapper for scripts that should be run in the linux cluster. The slurm wrapper script requires the directory to be analyzed, the number of tasks and the script to run. Mostly used for next generation sequencing analysis.

To use the slurm wrapper script:
```{bash}
~/lab_utils/linux_cluster/000_slurmWrapper.sh <array> <script_to_run> <directory_to_process>
```
The script displays a useful usage message that describes the arguments in more detail.

## LOGGING
Most scripts output some sort of log file (stdout and stderr) that can be inspected with a text editor. The log files can usually be verified with vim ~/data/<dir>/logs/*_9004526_*_1.out.

## TAGS 
I have a set of tags that I try to use to put marks on code for future reference. The form of the tags is <comment><TAG>. recursive (-r) grep can be used to find the tags.
## LOGGING
Most scripts output some sort of log file (stdout and stderr) that can be inspected with a text editor. The log files can usually be verified with vim ~/data/<dir>/logs/*_9004526_*_1.out

## TAGS 
I have a set of tags that I try to use to put marks on code for future reference. The form of the tags is <comment><TAG>. grep can be used to find the tags.

TODO: Tasks that I have to complete for that particular code file. 
HOWTO: Designates different code snippets for reference when I want to see how to do a particular thing.
FIXME: Highlight areas that need fixing.
NOTE: Add important notes or explanations.
BUG: Mark known bugs or issues.
OPTIMIZE: Indicate areas that could be optimized for better performance.
REFACTOR: for code that needs refactoring
TEST: for testing purposes

## TROUBLESHOOTING 
Each documentation section has a troubleshooting section that lets the user know about common errors that could be encountered, such as the scripts depending on the name of the files.

## TODO
Each documentation for each subdirectory has a # TODO section 

## STICKY_NOTES.md
Notes I take while developing the scripts. 
