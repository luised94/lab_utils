#!/bin/bash



#Define arrays
DIR_TO_PROCESS="$HOME/data/240304Bel"
REFGENOME_DIR="$HOME/data/REFGENS"

mapfile -t fastq_paths < <(find "${DIR_TO_PROCESS}" -type f -name "processed_*.fastq" )
mapfile -t genomes_paths < <(find "${REFGENOME_DIR}" -type f -name "*_refgenome.fna")

#Total number if jobs is the product of the number of genomes and the number of FASTQ files
total_jobs=$(( ${#genomes_paths[@]} * ${#fastq_paths[@]} ))

# Print header
echo "SLURM_ARRAY_TASK_ID | GENOME_INDEX | FASTQ_INDEX | GENOME_PATH | FASTQ_PATH"

# Loop through each SLURM task ID
for (( task_id=1; task_id<=total_jobs; task_id++ )); do
    # Calculate the genome and FASTQ indices
    #Slurm is one-based but bash is 0-based
    genome_index=$(( (task_id - 1) / ${#fastq_paths[@]} ))
    fastq_index=$(( (task_id - 1) % ${#fastq_paths[@]} ))
    fastq_ID=$(echo "${fastq_paths[$fastq_index]}" | cut -d_ -f3 )
    genome_name=$( echo "${genomes_paths[$genome_index]}" | cut -d_ -f1 | rev | cut -d/ -f1 | rev )

    # Print the indices and corresponding paths
    echo "$task_id | $genome_index | $fastq_index" #| $(basename ${genomes_paths[$genome_index]%_refgenome.fna}_index) | $(basename ${fastq_paths[$fastq_index]})" 
    echo "${fastq_ID}_${genome_name}.sam"
done
