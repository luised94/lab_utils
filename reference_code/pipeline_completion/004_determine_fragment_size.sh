#!/usr/bin/env bash
################################################################################
# Determine fragment size after deduplication and blacklist removal
################################################################################
# Purpose: Do preprocessing before determining fragment size
# Usage: Run as script.
# $./004_determine_fragment_size.sh
# DEPENDENCIES: bash, data from 250207Bel BMC experiment copied to preprocessing_test directory
# OUTPUT: Files
# NOTES: Updated the files being used in analysis, went from the data in the 241010Bel to the data in 250207Bel
# AUTHOR: LEMR
# DATE: 2025-04-24
################################################################################
SCRIPT_NAME=$(basename "$0")
CURRENT_TIMESTAMP=$(date +'%Y-%m-%d %H:%M:%S')
echo "=========================================="
echo "Script: $SCRIPT_NAME"
echo "Start Time: $CURRENT_TIMESTAMP"
echo "=========================================="

# Output config
OUTDIR="$HOME/data/preprocessing_test"
SUB_DIRS=("align" "predictd" "peaks" "coverage")
# Create directory structure
mkdir -p "${SUB_DIRS[@]/#/$OUTDIR/}"

# Genome data
GENOME_SIZE=12000000  # 1.2e7 in integer form

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
# Initialization of array with keys and file paths
# Store find results in a variable
mapfile -t FILES < <( find ~/data/preprocessing_test/align -maxdepth 1 -type f -name "*_blFiltered.bam" | sort )

if [[ "${#FILES}"  -eq 0 ]];
then
    echo "No *_blFiltered.bam files found in ~/data/preprocessing_test/align" >&2
    exit 1
fi

echo -e "\n=== Analyzing Fragment Sizes ==="
for filepath in "${FILES[@]}";
do
  filename=$(basename filepath)
  sample_name=${filename%.bam}
  sample_type=$(echo "$sample_name" | cut -d_ -f1-2 )
  log_file="$OUTDIR/predictd/${sample_type}_macs2.log"
  echo "--- Sample information ---"
  echo "Filename: $filename"
  echo "Sample name: $sample_name"
  echo "Sample type: $sample_type"
  echo "Log file: $filename"

  # Only run MACS2 if log is missing/incomplete
  if [[ ! -f "$log_file" ]] || ! grep -q 'predicted fragment length is' "$log_file"; then
    echo "Running fragment size prediction..."
    macs2 predictd \
        --ifile "$filepath" \
        --mfold 2 200 \
        --gsize "$GENOME_SIZE" \
        --outdir "$OUTDIR/predictd" 2> "$log_file"
  else
    echo "Using existing prediction results in: $log_file"
  fi

  # Validate and extract fragment size
  if ! frag_size=$(grep -oP 'predicted fragment length is \K\d+' "$log_file"); then
    echo -e "\nERROR: Fragment analysis failed for $sample_type"
    echo "Debug info:"
    echo "Log file: $log_file"
    [[ -f "$log_file" ]] && tail -n 20 "$log_file"
    exit 4
  fi
  echo "Fragment size: $frag_size"

done # end fragment size analysis for loop
