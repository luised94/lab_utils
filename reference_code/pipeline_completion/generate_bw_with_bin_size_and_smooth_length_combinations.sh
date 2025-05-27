#!/usr/bin/env bash
################################################################################
# Calculate bw files with different bin sizes and smoothLength combinations
################################################################################
# Purpose: Generate bw to plot genomic tracks
# Usage: Run as script.
# $./006_calculate_bw_coverage_for_preprocessed_bam_files.sh
# DEPENDENCIES: bash, data from 250207Bel BMC experiment and scripts 001-005
# OUTPUT: Duplicate files with renaming for easier identification and isolation from data directories
# NOTES: Updated the files being used in analysis, went from the data in the 241010Bel to the data in 250207Bel
# AUTHOR: LEMR
# DATE: 2025-02-07
# UPDATE: 2025-04-22
################################################################################
SCRIPT_NAME=$(basename "$0")
CURRENT_TIMESTAMP=$(date +'%Y-%m-%d %H:%M:%S')
echo "=========================================="
echo "Script: $SCRIPT_NAME"
echo "Start Time: $CURRENT_TIMESTAMP"
echo "=========================================="
module load python/2.7.13 deeptools/3.0.1
echo "Loaded deeptools module..."
#set -euo pipefail

# Reusable component start ===
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
# Reusable component end ===

# Output config
OUTDIR="$HOME/data/preprocessing_test"
OUTPUT_DIR="$OUTDIR/coverage/"
SUB_DIRS=("align" "predictd" "peaks" "coverage")
# Create directory structure
mkdir -p "${SUB_DIRS[@]/#/$OUTDIR/}"

# Define normalization methods
# Raw, RPKM, and CPM
# leave as raw
declare -a NORMALIZATION_METHOD=("raw")  
THREADS=8 # int, number of cores
#GENOME_SIZE=12000000  # 1.2e7 in integer form, bp
# Test multiple combinations of bin size and lengths
BIN_SIZES=(10 25 50 100 200 500)
SMOOTH_LENGTHS=(25 50 100 200 500 1000 2000)

# Initialization of array with keys and file paths
# Create an array of files
mapfile -t FILES < <(find ~/data/preprocessing_test/align -maxdepth 1 -type f -name "*_raw*.bam" | sort)

# Check if any files were found
if [ ${#FILES[@]} -eq 0 ];
then
  echo "No *.bam files found in ~/data/preprocessing_test/align" >&2
  exit 1
fi

# Process the files if found
for filepath in "${FILES[@]}"; do
  filename=$(basename "$filepath")
  key=${filename%.bam}
  sample_type=$(echo "$key" | cut -d'_' -f1)
  sample_name=$(echo "$key" | cut -d'_' -f1-2)
  # TODO: Need to adjust this. Maybe determine the amount of underscores or just take the last available.
  bam_type=$(echo "$key" | awk -F_ '{ print $NF }')
  #bam_type=$(echo "$key" | cut -d'_' -f3 )

  echo "---Sample information---"
  echo "  Input Path: $filepath"
  echo "  Filename: $filename"
  echo "  Sample_name: $sample_name"
  echo "  Sample type: $sample_type"
  echo "  Bam type: $bam_type"

  for norm_method in "${NORMALIZATION_METHOD[@]}";
  do
    norm_suffix=${norm_method,,}
    base_flags=(
        --bam "$filepath"
        --numberOfProcessors $(( THREADS / 2))
        #--effectiveGenomeSize "$GENOME_SIZE"
        #--ignoreDuplicates
    )
    # Add normalization flag only if not raw (case-insensitive comparison)
    [[ "$norm_suffix" != "raw" ]] && base_flags+=(--normalizeUsing "$norm_method")
    echo "  ---Normalization information---"
    echo "    Normalization method norm_suffix: $norm_suffix"

    # --- Case 1: Process WITHOUT blacklist ---
    #output_name_no_bl="${OUTPUT_DIR}${key}_${norm_suffix}_noBlacklist.bw"
    for bin_size in "${BIN_SIZES[@]}";
    do
      echo "    ---Bin size information ---"
      echo "    Bin size: $bin_size"
      #+=(--binSize "$bin_size")
      flags_with_bin=("${base_flags[@]}" --binSize "$bin_size")
      for smooth_length in "${SMOOTH_LENGTHS[@]}"
      do
        output_name="${OUTPUT_DIR}${key}_${norm_suffix}_${bin_size}_${smooth_length}.bw"
        echo "    ---Smooth Length information ---"
        echo "    Smooth Length: $smooth_length"
        echo "    Current output name: $output_name"
        flags_with_smooth_length=("${flags_with_bin[@]}" --smoothLength "$smooth_length")
        #base_flags+=(--smoothLength "$smooth_length")
        echo "    Output without blacklist: $output_name"
        if [[ -f "$output_name" ]];
        then
            echo "Output file already exists. Skipping: $output_name"
        else
          # Create flags specifically for the NO blacklist run, adding the output file
          flags_no_bl=("${flags_with_smooth_length[@]}" -o "$output_name")
          echo "Executing coverage calculation without blacklist..."
          # Print the command that would be executed
          echo "COMMAND: bamCoverage ${flags_no_bl[*]}"
          # Execute with all flags properly expanded
          #bamCoverage "${flags_no_bl[@]}"
        fi
      done
    done
    # --- Case 2: Process with blacklist ---
    # Blacklist processing was moved to after deduplication
    #if [[ "$BLACKLIST_EXISTS" == true ]]; then
    #  output_name_with_bl="${OUTPUT_DIR}${key}_${norm_suffix}_withBlacklist.bw"
    #  if [[ -f "$output_name_with_bl" ]] ; then
    #      echo "Output file already exists. Skipping: $output_name_with_bl"
    #  else
    #    echo "    Output with blacklist: $output_name_with_bl"
    #    echo "Executing coverage calculation with blacklist..."
    #    flags_with_bl=("${base_flags[@]}" -o "$output_name_with_bl" --blackListFileName "$BLACKLIST_BED_FILE")
    #    # Print the command that would be executed
    #    echo "COMMAND: bamCoverage ${flags_with_bl[*]}"
    #    # Execute with all flags properly expanded
    #    bamCoverage "${flags_with_bl[@]}"
    #  fi
    #else
    #  echo "Blacklist file not found. Skipping: $output_name_with_bl"
    #fi

  done # End of norm_method loop
done # End of filepath loop

# Create array of bigwig files
#declare -a BIGWIG_FILES
#mapfile -t BIGWIG_FILES < <(find "$OUTDIR/coverage" -name "*.bw" | sort)
#
## Validate array
#if [[ ${#BIGWIG_FILES[@]} -eq 0 ]];
#then
#  echo "ERROR: No bigwig files found in $OUTDIR/coverage" >&2
#  echo "Verify the directory or 003_calculate_bw_coverage_for_preprocessed_bam_files.sh" >&2
#  #exit 3
#fi
echo "Processed ${#FILES[@]} bam files"
echo "Finished generating bigwig files..."
