#!/usr/bin/env bash
################################################################################
# Duplicate files for testing pipeline completion
################################################################################
# Purpose: Define files that will be duplicated from their location in the ~/data/<experiment_id>/. Files will be renamed to set the metadata.
# Usage: Run as script.
# $./duplicate_files_for_testing.sh
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

    if [[ ! -e "$input" ]]; then
        echo "[ERROR] Missing $sample_type file: $input" 1>&2
        missing_files=1
    else
        echo "Will copy. Bam file found: $input"
        [[ ! -e "$copy_path" ]] && { cp "$input" "$copy_path"; }
    fi
done

(( missing_files )) && { echo "Aborting: Missing bam files"; echo "Use another experiment for testing or download Eaton 2010 data."; exit 1; }
echo "All bam files copied."

#declare -A FASTQ_FILES=(
#    ['test']="$HOME/data/250207Bel/fastq/consolidated_245018_sequence.fastq"
#    ['input']="$HOME/data/250207Bel/fastq/consolidated_245003_sequence.fastq"
#    ['reference']="$HOME/data/100303Bel/fastq/consolidated_034475_sequence.fastq"
#)
#
#declare -i missing_files=0
#for sample_type in "${!FASTQ_FILES[@]}";
#do
#    input="${FASTQ_FILES[$sample_type]}"
#    copy_path="$OUTDIR/fastq/${sample_type}_raw.fastq"
#    echo "Processing $sample_type..."
#    echo "Input: $input"
#    echo "Copy: $copy_path"
#
#    if [[ ! -e "$input" ]]; then
#        echo "[ERROR] Missing $sample_type file: $input" 1>&2
#        missing_files=1
#    else
#        echo "Will copy. Fastq file found: $input"
#        [[ ! -e $copy_path ]] && { cp $input $copy_path; }
#    fi
#done
#
#(( missing_files )) && { echo "Aborting: Missing fastq files"; echo "Use another experiment for testing or download Eaton 2010 data."; exit 1; }
echo "Script complete."
echo "Remember to use samtools to index the raw bam files."
echo "See file for the for loop"
#for bam in *_raw.bam; do samtools index "$bam" ; done
