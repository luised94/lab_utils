#!/bin/bash
#USAGE: From anywhere, run '~/data/lab_utils/next_generation_sequencing/FTQPRC/000_directoryCreation.sh'
#ALTERNATIVE: node_bt2build_refgenomes.sh was created as a workaround but is unnecessary now. Uses for loop for bowtie2-build. File moved to archive.

data_dir=$HOME/data/
mapfile -t bel_dirs < <(find "$data_dir" -maxdepth 1 -type d -name "*Bel*")

for dir in "${bel_dirs[@]}"; do
  echo "Creating subdirs for $dir"	
  mkdir -p "$dir"/{peak,fastq,alignment,qualityControl,processedFastq,bigwig,plots,logs,documentation}
done

