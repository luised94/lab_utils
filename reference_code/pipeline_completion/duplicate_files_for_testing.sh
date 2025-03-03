#!/usr/bin/env bash
echo "Starting script for copying bam and fastq files for pipeline testing"

declare -A SAMPLES=(
    ['test']="$HOME/data/241010Bel/alignment/processed_245018_sequence_to_S288C_sorted.bam"
    ['input']="$HOME/data/241010Bel/alignment/processed_245003_sequence_to_S288C_sorted.bam"
    ['reference']="$HOME/data/100303Bel/alignment/processed_034475_sequence_to_S288C_sorted.bam"
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
        [[ ! -e $copy_path ]] && { cp $input $copy_path; }
    fi
done

(( missing_files )) && { echo "Aborting: Missing bam files"; echo "Use another experiment for testing or download Eaton 2010 data."; exit 1; }
echo "Script complete."

declare -A FASTQ_FILES=(
    ['test']="$HOME/data/241010Bel/fastq/consolidated_245018_sequence.fastq"
    ['input']="$HOME/data/241010Bel/fastq/consolidated_245003_sequence.fastq"
    ['reference']="$HOME/data/100303Bel/fastq/consolidated_034475_sequence.fastq"
)

declare -i missing_files=0
for sample_type in "${!FASTQ_FILES[@]}";
do
    input="${FASTQ_FILES[$sample_type]}"
    copy_path="$OUTDIR/fastq/${sample_type}_raw.fastq"
    echo "Processing $sample_type..."
    echo "Input: $input"
    echo "Copy: $copy_path"

    if [[ ! -e "$input" ]]; then
        echo "[ERROR] Missing $sample_type file: $input" 1>&2
        missing_files=1
    else
        echo "Will copy. Fastq file found: $input"
        [[ ! -e $copy_path ]] && { cp $input $copy_path; }
    fi
done

(( missing_files )) && { echo "Aborting: Missing fastq files"; echo "Use another experiment for testing or download Eaton 2010 data."; exit 1; }
echo "Script complete."
