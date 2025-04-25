#!/usr/bin/env bash
################################################################################
# Remove duplicates and blacklist reads in bam files
################################################################################
# Purpose: Generate intermediate bam files for raw bam files after deduplication and blacklist files
# Usage: Run as script.
# $./003_deduplicate_and_deblacklist_bam_files.sh
# DEPENDENCIES: bash, data from 250207Bel BMC experiment copied to preprocessing_test directory
# OUTPUT: Files 
# NOTES: Updated the files being used in analysis, went from the data in the 241010Bel to the data in 250207Bel
# AUTHOR: LEMR
# DATE: 2025-04-17
# UPDATE: 2025-02-07
################################################################################
SCRIPT_NAME=$(basename "$0")
CURRENT_TIMESTAMP=$(date +'%Y-%m-%d %H:%M:%S')
echo "=========================================="
echo "Script: $SCRIPT_NAME"
echo "Start Time: $CURRENT_TIMESTAMP"
echo "=========================================="
#set -euo pipefail

# Reusable component start
# === Cluster Environment Setup ===
# Check if running on head node (should use interactive job instead)
if [[ "$(hostname)" == "luria" ]]; then
    echo "Error: This script should not be run on the head node" 1>&2
    echo "Please run: srun --pty bash" 1>&2
    exit 1
fi

# Check if running inside a Slurm allocation
if [[ -z "${SLURM_JOB_ID}" ]]; then
    echo "Error: This script must be run inside a Slurm allocation" 1>&2
    echo "Please run: srun --pty bash" 1>&2
    exit 1
fi
# Reusable component end

# Load required modules
module load picard/2.18.26 java samtools bedtools python/2.7.13 deeptools/3.0.1
export PICARD_JAR="/home/software/picard/picard-2.18.26/picard.jar"

if [[ ! -f "$PICARD_JAR" ]]; then
  echo "[ERROR] Picard JAR missing: $PICARD_JAR" 1>&2
  exit 10
fi

#========================================
# Cluster config
#THREADS=8

BLACKLIST_BED_FILE="$HOME/data/feature_files/20250423_merged_saccharomyces_cerevisiae_s288c_blacklist.bed"
# Check if blacklist file exists once before the loops
if [[ ! -f "$BLACKLIST_BED_FILE" ]];
then
    echo "Blacklist file not found: $BLACKLIST_BED_FILE. Blacklist runs will be skipped."
    exit 1
fi

# Output config
OUTDIR="$HOME/data/preprocessing_test"
SUB_DIRS=("align" "predictd" "peaks" "coverage" "metrics")
# Create directory structure
mkdir -p "${SUB_DIRS[@]/#/$OUTDIR/}"

# Initialization of array with file paths
# Store find results in a variable
#mapfile -d '' FILES < <(find ~/data/preprocessing_test/align -maxdepth 1 -type f -name "*_raw.bam" -print0)
mapfile -t FILES < <( find ~/data/preprocessing_test/align -maxdepth 1 -type f -name "*_raw.bam" | sort )

# Step 2: Check if any files were found
# See ./lab_utils/reference_code/pipeline_completion/duplicate_files_for_testing.sh to see what samples should be initialized.
# Alternatively verify the directory with the find command.
if [[ "${#FILES}"  -eq 0 ]];
then
    echo "No *_raw.bam files found in ~/data/preprocessing_test/align" >&2
    exit 1
fi

##############################################
# Preprocessing Workflow
##############################################
# === Step 1: Duplicate Removal ===
# === Step 1: Index Deduplicated BAMs ===
echo -e "\n=== Deduplicating and removing reads in blacklist ==="
for filepath in "${FILES[@]}";
do
  filename=$(basename "$filepath")
  sample_name=${filename%.bam}
  deduplicated_bam="$OUTDIR/align/${sample_name}_deduped.bam"
  index_file_of_deduplicated_bam="${deduplicated_bam}.bai"
  blacklist_filtered_bam="${deduplicated_bam%.bam}_blFiltered.bam"
  index_file_of_blacklist_filtered_bam="${blacklist_filtered_bam}.bai"

  echo "--- Sample Information ---"
  echo "Filepath: $filename"
  echo "Sample_name: $sample_name"
  echo "Deduplicated bam: $deduplicated_bam"
  echo "Deduplicated bam index file: $index_file_of_deduplicated_bam"
  echo "Blacklist bam: $blacklist_filtered_bam"
  echo "Blacklist bam index file: $index_file_of_blacklist_filtered_bam"

  # === Deduplicate Bam files ===
  if [[ ! -f "$deduplicated_bam" ]];
  then
    java -jar $PICARD_JAR MarkDuplicates \
        I="$filepath" \
        O="$deduplicated_bam" \
        M="$OUTDIR/metrics/${sample_name}_dup_metrics.txt" \
        REMOVE_DUPLICATES=true
  else
    echo "Skipping existing: $deduplicated_bam"
  fi

  # Validate output
  if [[ ! -s "$deduplicated_bam" ]];
  then
    echo "ERROR: Failed to create $deduplicated_bam" >&2
    exit 2
  fi

  echo "Indexing deduped BAM..."
  if [[ ! -f "$index_file_of_deduplicated_bam" ]];
  then
    echo "Indexing: $deduplicated_bam"
    if ! samtools index "$deduplicated_bam" 2>&1;
    then
      echo "[ERROR] Indexing failed for ${deduplicated_bam} " >&2
      exit 3
    fi
    echo "File $deduplicated_bam indexed succesfully"
  else
    echo "Index exists: $index_file_of_deduplicated_bam"
  fi

  # === Filter reads in blacklist regions ===
  if [[ ! -f "$blacklist_filtered_bam" ]];
  then
    alignmentSieve --bam "$deduplicated_bam" \
      --blackListFileName "$BLACKLIST_BED_FILE" \
      --outputBamFile "$blacklist_filtered_bam" \
      --outFileFormat BAM \
      --numberOfProcessors auto
  else
    echo "Skipping existing: $blacklist_filtered_bam"
  fi

  # Validate output
  if [[ ! -s "$blacklist_filtered_bam" ]];
  then
    echo "ERROR: Failed to create $blacklist_filtered_bam" >&2
    exit 2
  fi

  echo "Indexing blacklist filtered BAM..."
  if [[ ! -f "$index_file_of_blacklist_filtered_bam" ]];
  then
    echo "Indexing: $blacklist_filtered_bam"
    samtools index "$blacklist_filtered_bam"
    if ! samtools index "$blacklist_filtered_bam" 2>&1;
    then
      echo "[ERROR] Indexing failed for ${blacklist_filtered_bam} " >&2
      exit 3
    fi
    echo "File $blacklist_filtered_bam indexed succesfully"
  else
    echo "Index exists: $index_file_of_blacklist_filtered_bam"
  fi

done # end deduplicate and blacklist filter loop
