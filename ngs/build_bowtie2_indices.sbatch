#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --nice=10000
#SBATCH --exclude=c[5-22]
#SBATCH --output=/home/%u/data/logs/slurm_%A_%a.out
#SBATCH --error=/home/%u/data/logs/slurm_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luised94@mit.edu
#===============================================================================
# TITLE: Build Bowtie2 Indices for Reference Genomes
# DESCRIPTION: Creates bowtie2 indices for genome files
# USAGE: 
# Find out how many genomes you have
#   ls ~/data/reference_genomes/*_genome.fna | wc -l
# Then run:
#   sbatch --array=0-N build_bowtie2_indices.sh
#     where N is (number of genomes - 1)
#===============================================================================
set -euo pipefail

# Configuration
reference_genomes_directory="$HOME/data/reference_genomes"
genome_file_pattern="*_genome.fna"
threads_count="$SLURM_CPUS_PER_TASK"

echo "========================================"
echo "Bowtie2 Index Building"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Reference directory: $reference_genomes_directory"
echo "Threads: $threads_count"
echo ""

# Find all genome files
mapfile -t genome_files < <(find "$reference_genomes_directory" -maxdepth 1 -name "$genome_file_pattern" | sort)

# Check if any genome files found
genome_count="${#genome_files[@]}"
if [[ $genome_count -eq 0 ]]; then
    echo "ERROR: No genome files found matching pattern: $genome_file_pattern"
    echo "Search directory: $reference_genomes_directory"
    exit 1
fi

echo "Found $genome_count genome files:"
for genome_file in "${genome_files[@]}"; do
    echo "  - $(basename "$genome_file")"
done
echo ""

# Select genome for this array task
if [[ $SLURM_ARRAY_TASK_ID -ge $genome_count ]]; then
    echo "ERROR: Array task ID ($SLURM_ARRAY_TASK_ID) exceeds genome count ($genome_count)"
    echo "Submit with: sbatch --array=0-$((genome_count - 1)) $0"
    exit 1
fi

selected_genome="${genome_files[$SLURM_ARRAY_TASK_ID]}"
echo "Processing genome $((SLURM_ARRAY_TASK_ID + 1)) of $genome_count"
echo "Genome file: $(basename "$selected_genome")"
echo ""

# Build index base name (remove _genome.fna suffix)
genome_basename=$(basename "$selected_genome")
index_base="${selected_genome%_genome.fna}_index"
index_base_name=$(basename "$index_base")

echo "Index base name: $index_base_name"
echo ""

# Check if index already exists
first_index_file="${index_base}.1.bt2"
if [[ -f "$first_index_file" ]]; then
    echo "Bowtie2 index already exists for this genome"
    echo "Found: $(basename "$first_index_file")"
    echo "Skipping index build"
    exit 0
fi

# Build bowtie2 index
echo "Starting bowtie2-build..."
echo "Command: bowtie2-build --threads $threads_count $selected_genome $index_base"
echo ""

start_time=$(date +%s)

if bowtie2-build --threads "$threads_count" "$selected_genome" "$index_base"; then
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))

    echo ""
    echo "========================================"
    echo "Index build completed successfully"
    echo "========================================"
    echo "Elapsed time: ${elapsed_time} seconds"
    echo "Index files created:"
    ls -lh "${index_base}".*.bt2
else
    echo ""
    echo "ERROR: bowtie2-build failed for: $genome_basename"
    exit 1
fi
