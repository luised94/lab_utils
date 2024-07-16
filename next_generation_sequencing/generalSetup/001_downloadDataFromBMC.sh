#!/bin/bash
#DESCRIPTION: Download the fastq file from the BMC servers. Ignore bams, htmls, fastqc etc. Only the BMC processed fastq pipeline.
#USAGE:

bmc_server=$1
bell_lab_directory=$2
srun rsync -av -n --include '*/' --include '*.fastq' --exclude '*' /net/${bmc_server}/data/bmc/public/Bell/${bell_lab_directory}/ ~/data/${bell_lab_directory}/fastq/
cd ~/data/${bell_lab_directory}/fastq/
find . -type f -name "*.fastq" -exec mv {} . \;
find . -type d \( -name "*D24*" -o -name "infosite*" \) -exec rm -rf {} \;
rm ~/data/${bell_lab_directory}/fastq/*unmapped*.fastq
