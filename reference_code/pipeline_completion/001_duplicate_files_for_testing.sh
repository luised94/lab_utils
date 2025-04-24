#!/usr/bin/env bash
################################################################################
# Duplicate files for testing pipeline completion
################################################################################
# Purpose: Define files that will be duplicated from their location in the ~/data/<experiment_id>/. Files will be renamed to set the metadata.
# Usage: Run as script.
# $./001_duplicate_files_for_testing.sh
# DEPENDENCIES: bash, data from 250207Bel BMC experiment
# OUTPUT: Duplicate files with renaming for easier identification and isolation from data directories
# NOTES: Updated the files being used in analysis, went from the data in the 241010Bel to the data in 250207Bel
# AUTHOR: LEMR
# DATE: 2025-04-17
# UPDATE: 2025-02-07
################################################################################
echo "Starting script for copying bam and fastq files for pipeline testing"
declare -A SAMPLES=(
  # WT_NONE_NOCO_HM1108_2_positive
  ['test_001']="$HOME/data/250207Bel/alignment/processed_126050_sequence_to_S288C_sorted.bam"
  # 4R_4PS_NOCO_HM1108_2_treatment
  ['test_002']="$HOME/data/250207Bel/alignment/processed_126058_sequence_to_S288C_sorted.bam"
  # WT_NONE_ALPHA_UM174_2_positive
  ['test_003']="$HOME/data/250207Bel/alignment/processed_126065_sequence_to_S288C_sorted.bam"
  # WT_NONE_NOCO_UM174_2_positive
  ['test_004']="$HOME/data/250207Bel/alignment/processed_126066_sequence_to_S288C_sorted.bam"
  # 4R_NONE_ALPHA_Input_2_input
  ['test_005']="$HOME/data/250207Bel/alignment/processed_126042_sequence_to_S288C_sorted.bam"

  # WT_NONE_ALPHA_Input_2_input
  ['input_001']="$HOME/data/250207Bel/alignment/processed_126041_sequence_to_S288C_sorted.bam"
  # Eaton2010
  ['reference_001']="$HOME/data/100303Bel/alignment/processed_034475_sequence_to_S288C_sorted.bam"
)

# Output config
OUTDIR="$HOME/data/preprocessing_test"
SUB_DIRS=("align" "predictd" "peaks" "coverage" "fastq" "quality_control" "plots")

# Create directory structure
mkdir -p "${SUB_DIRS[@]/#/$OUTDIR/}"

declare -i missing_files=0
for sample_type in "${!SAMPLES[@]}";
do
    input="${SAMPLES[$sample_type]}"
    copy_path="$OUTDIR/align/${sample_type}_raw.bam"
    echo "Processing $sample_type..."
    echo "Input: $input"
    echo "Copy: $copy_path"

    if [[ ! -e "$input" ]];
    then
        echo "[ERROR] Missing $sample_type file: $input" 1>&2
        missing_files=1
    else
        echo "Will copy. Bam file found: $input"
        # Copy if it does not exist
        [[ ! -e "$copy_path" ]] && { cp "$input" "$copy_path"; }
    fi
done

(( missing_files )) && {
  echo "Aborting: Missing bam files";
  echo "Use another experiment for testing or download Eaton 2010 data.";
  exit 1;
}

echo "All bam files copied."
mapfile -t BAM_FILES < <(find "$OUTDIR/align" -type f -name "*.bam")
if [[ ${#SAMPLES} -eq ${#BAM_FILES} ]];
then
  echo "[WARNING] Number of samples in array does not equal number of bam files in the directory."
  echo "[SUGGESTION] Double check the files in the directory. Could be leftover from previous analysis."
fi

echo "Indexing raw bam files..."
mapfile -t RAW_BAM < <(find "$OUTDIR/align" -type f -name "*_raw.bam")
module add samtools
for bam in "${RAW_BAM[@]}"; do samtools index "$bam" ; done
mapfile -t BAI_FILES < <(find "$OUTDIR/align" -type f -name "*.bam.bai")
if [[ ${#BAI_FILES} -eq ${#BAM_FILES} ]];
then
  echo "[WARNING] Number of index files does not equal number of bam files in the directory."
  echo "[SUGGESTION] Double check the files in the directory. Could be leftover from previous analysis."
fi

echo "Script complete."
echo "See file for the for loop"
