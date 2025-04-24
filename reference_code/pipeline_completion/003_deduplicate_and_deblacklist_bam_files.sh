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

# === Cluster Environment Setup ===
#if [[ "$(hostname)" != "luria" ]]; then
#    echo "Error: This script must be run on luria cluster" 1>&2
#    exit 1
#fi

# Load required modules
module load picard/2.18.26 java samtools bedtools
export PICARD_JAR="/home/software/picard/picard-2.18.26/picard.jar"

if [[ ! -f "$PICARD_JAR" ]]; then
  echo "[ERROR] Picard JAR missing: $PICARD_JAR" 1>&2
  exit 10
fi

# Initialize Conda and MACS2 environment
CONDA_ROOT=~/miniforge3

if [[ -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
  . "$CONDA_ROOT/etc/profile.d/conda.sh"
else
  echo "ERROR: Conda not found at $CONDA_ROOT" >&2
  exit 1
fi

# Activate MACS2 environment with validation
conda activate macs2_env 2> /dev/null

# Verify MACS2 installation
if ! command -v macs2 &> /dev/null; then
  echo "ERROR: MACS2 not installed. Install with:" >&2
  echo "conda install -c bioconda macs2=2.2.7.1" >&2
  exit 2
fi

echo "MACS2 environment ready in $(which python)"
#source ~/lab_utils/core_scripts/setup_conda_and_macs2.sh || exit 1
#======================================== 
# Cluster config
#THREADS=8
# Genome data
GENOME_FASTA="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C/SaccharomycescerevisiaeS288C_refgenome.fna"
if [[ ! -f $GENOME_FASTA ]]; then
  echo "$GENOME_FASTA does not exist."
  exit 1
fi

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

# Initialization of array with keys and file paths
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

#SELECTED_SAMPLE_KEYS=$(echo "${!SAMPLES[@]}" | sed 's/ /\n/g' | grep -v "input")
#mapfile -t SELECTED_SAMPLE_KEYS < <(echo "${!SAMPLES[@]}" | tr ' ' '\n' | grep -E '^test|^reference')

##############################################
# Preprocessing Workflow
##############################################
# === Step 1: Duplicate Removal ===
# === Step 1: Index Deduplicated BAMs ===
echo -e "\n=== Deduplicating and removing reads in blacklist ==="
for filepath in "${FILES[@]}";
do
  filename=$(basename "$filepath")
  key=${filename%.bam}
  deduplicated_bam="$OUTDIR/align/${key}_deduped.bam"
  index_file_of_deduplicated_bam="${deduplicated_bam}.bai"
  blacklist_filtered_bam="${deduplicated_bam%.bam}_blFiltered.bam"
  index_file_of_blacklist_filtered_bam="${blacklist_filtered_bam}.bai"

  echo "--- Sample Information ---"
  echo "Filepath: $filename"
  echo "Key: $key"
  echo "Deduplicated bam: $deduplicated_bam"
  echo "Deduplicated bam index file: $index_file_of_deduplicated_bam"
  echo "Blacklist bam: $blacklist_filtered_bam"
  echo "Blacklist bam index file: $index_file_of_blacklist_filtered_bam"

  # === Deduplicate Bam files ===
  #if [[ -f "$deduplicated_bam" ]];
  #then
  #  echo "Skipping existing: $deduplicated_bam"
  #else
  #  java -jar $PICARD_JAR MarkDuplicates \
  #      I="$filepath" \
  #      O="$deduplicated_bam" \
  #      M="$OUTDIR/metrics/${key}_dup_metrics.txt" \
  #      REMOVE_DUPLICATES=true
  #fi

  ## Validate output
  #if [[ ! -s "$deduplicated_bam" ]];
  #then
  #  echo "ERROR: Failed to create $deduplicated_bam" >&2
  #  exit 2
  #fi

  #echo "Indexing deduped BAM..."
  #if [[ -f "$index_file_of_deduplicated_bam" ]];
  #then
  #  echo "Index exists: $index_file_of_deduplicated_bam"
  #else
  #  echo "Indexing: $deduplicated_bam"
  #  samtools index "$deduplicated_bam"
  #  exit_code=$?
  #  if [[ $exit_code -eq 0 ]] && [[ -f "$index_file_of_deduplicated_bam" ]]; then
  #    echo "File $deduplicated_bam indexed succesfully"
  #  else
  #    echo "[ERROR] Indexing failed for ${deduplicated_bam} (exit $exit_code)" >&2
  #    exit 3
  #  fi
  #fi

  ## === Filter reads in blacklist regions ===
  #if [[ -f "$blacklist_filtered_bam" ]];
  #then
  #  echo "Skipping existing: $blacklist_filtered_bam"
  #else
  #  alignmentSieve --bam "$deduplicated_bam" \
  #    --blackListFileName "$BLACKLIST_BED_FILE" \
  #    --outputBamFile "$blacklist_filtered_bam" \
  #    --outFileFormat BAM \
  #    --numberOfProcessors auto
  #fi

  ## Validate output
  #if [[ ! -s "$blacklist_filtered_bam" ]];
  #then
  #  echo "ERROR: Failed to create $blacklist_filtered_bam" >&2
  #  exit 2
  #fi

  #echo "Indexing blacklist filtered BAM..."
  #if [[ -f "$index_file_of_blacklist_filtered_bam" ]];
  #then
  #  echo "Index exists: $index_file_of_blacklist_filtered_bam"
  #else
  #  echo "Indexing: $blacklist_filtered_bam"
  #  samtools index "$blacklist_filtered_bam"
  #  exit_code=$?
  #  if [[ $exit_code -eq 0 ]] && [[ -f "$blacklist_filtered_bam" ]];
  #  then
  #    echo "File $blacklist_filtered_bam indexed succesfully"
  #  else
  #    echo "[ERROR] Indexing failed for ${blacklist_filtered_bam} (exit $exit_code)" >&2
  #    exit 3
  #  fi
  #fi

done # end deduplicate and index bam for loop
