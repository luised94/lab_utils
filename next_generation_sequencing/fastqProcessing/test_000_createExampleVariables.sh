#!/bin/bash
#Description: Initialize variables and modules to test in interactive node. 
#USAGE:Run source ~/data/lab_utils/next_generation_sequencing/FTQPRC/test_000_createExampleVariables.sh <dir>'

#SETUP
DIR_TO_PROCESS="$1"

DIR_TO_PROCESS="$HOME/data/$DIR_TO_PROCESS"
REFGENOME_DIR="$HOME/data/REFGENS"

#MODULE_LOAD
module purge
module load gnu/5.4.0
module load bowtie2/2.3.5.1
module load samtools/1.10
module load fastqc

#INITIALIZE_ARRAY
mapfile -t FASTQ_PATHS < <(find "${DIR_TO_PROCESS}" -type f -name "processed_*.fastq" )
mapfile -t GENOME_PATHS < <(find "${REFGENOME_DIR}" -type f -name "*_refgenome.fna")
mapfile -t BAM_PATHS < <(find "${DIR_TO_PROCESS}" -type f -name "*.bam")
