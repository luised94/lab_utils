#!/bin/bash

DIRECTORY="~/data/250930Bel/fastq"
mapfile -t fastq_files < <(find "$DIRECTORY" -type f -name "*.fastq")

cat ${fastq_files[6]} | head -n 8
cat ${fastq_files[7]} | head -n 8

cat ${fastq_files[6]} | wc -l
cat ${fastq_files[7]} | wc -l

# Example of filename: 250930Bel_D25-12496-2_2_sequence.fastq
for fastq_file in "${fastq_files[@]}"; do
  if [[ "$fastq_file" =~ "*unmapped*" ]]; then
    echo "Skipping unmapped fastq file"
    next
  fi
  echo "$fastq_file"

done
