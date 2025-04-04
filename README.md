# lab_utils
## TODO
Cross-correlation Analysis: Using phantompeakqualtools for strand cross-correlation (NSC/RSC metrics)
Peak Statistics: Using bedtools to analyze peak width distributions, distances between peaks, and genomic feature overlaps
Coverage Comparisons: Using deepTools multiBigwigSummary and plotCorrelation to compare signal profiles
Enrichment Analysis: Using GREAT or similar tools for genomic region enrichment
Motif Analysis: Using MEME-ChIP or HOMER for motif discovery in peaks
Peak Conservation: Analyzing conservation scores within peaks using phyloP/phastCons

## Overview
Code used for laboratory analysis
Most important scripts are found in core_scripts/ directory.

## Configuration
Because I perform the next-generation sequencing analysis using my institution's cluster, I have to use the version of the tools that are installed there for the most part. 
For this reason, I use R 4.2.0 to perform the analyses.

## Dependencies
The analysis are done locally or in a linux computing cluster. The linux cluster is the condition that dictates what dependencies are used, especially for the next generation sequencing analysis.

1. R 4.2.0
2. Command line utils
3. bowtie2/2.3.5.1
4. fastp/0.20.0
5. fastqc/0.11.5
6. deeptools/3.0.1
7. gatk
8. python/2.7
9. miniforge
10. macs2
11. picard
12. java
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

## LOGGING
Most scripts output some sort of log file (stdout and stderr) that can be inspected with a text editor. The log files can usually be verified with vim ~/data/<dir>/logs/*_9004526_*_1.out.

## TAGS 
I have a set of tags that I try to use to put marks on code for future reference. The form of the tags is <comment><TAG>. recursive (-r) grep can be used to find the tags.

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

## STICKY_NOTES.md
Notes I take while developing the scripts. 
