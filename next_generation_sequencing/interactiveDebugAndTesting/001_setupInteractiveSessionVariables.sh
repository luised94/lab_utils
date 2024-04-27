#!/bin/bash
#Description: Initialize variables and modules to test in interactive node. 
#USAGE:Run source ~/data/lab_utils/next_generation_sequencing/FTQPRC/test_000_createExampleVariables.sh <dir>'

#SETUP
DIR_TO_PROCESS="$HOME/data/$1"
REFGENOME_DIR="$HOME/data/REFGENS"

#MODULE_LOAD
module purge
module load gnu/5.4.0
module load bowtie2/2.3.5.1
module load samtools/1.10
module load fastqc
module load python/2.7.13
module load deeptools/3.0.1
module load r/4.2.0

#INITIALIZE_ARRAY
mapfile -t FASTQ_PATHS < <(find "${DIR_TO_PROCESS}" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \) | sort )
mapfile -t PROCESSED_FASTQ < <(find "${DIR_TO_PROCESS}" -type f -name "processed_*.fastq" | sort )
mapfile -t GENOME_PATHS < <(find "${REFGENOME_DIR}" -type f -name "*_refgenome.fna" )
mapfile -t BAM_PATHS < <(find "${DIR_TO_PROCESS}" -type f -name "*.bam" | sort)
mapfile -t BIGWIG_PATHS < <(find "${DIR_TO_PROCESS}" -type f -name "*S288C.bw" | sort )
