#!/usr/bin/env bash
################################################################################
# Shift bam files according to fragment size determined by MACS2
################################################################################
# Purpose: Use fragment size to predict fragment size
# Usage: Run as script.
# $./005_shift_bam_files.sh
# DEPENDENCIES: bash, chrom.sizes, shift_reads, samtools, 
# OUTPUT: chrom.sizes file
# NOTES:
# AUTHOR: LEMR
# DATE: 2025-04-24
# UPDATE:
################################################################################
SCRIPT_NAME=$(basename "$0")
CURRENT_TIMESTAMP=$(date +'%Y-%m-%d %H:%M:%S')
echo "=========================================="
echo "Script: $SCRIPT_NAME"
echo "Start Time: $CURRENT_TIMESTAMP"
echo "=========================================="
module load samtools
echo "Loaded samtools module..."
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
OUTDIR="$HOME/data/preprocessing_test"
SUB_DIRS=("align" "predictd" "peaks" "coverage")
# Create directory structure
mkdir -p "${SUB_DIRS[@]/#/$OUTDIR/}"
THREADS=$( nproc )

CHROM_SIZES_FILE="$OUTDIR/chrom.sizes"
AWK_SHIFT_FILE="$HOME/lab_utils/core_scripts/shift_reads.awk"

# Check for required files
if [[ ! -f "$CHROM_SIZES_FILE" ]]; then
  echo "[ERROR] ${CHROM_SIZES_FILE} not present."
  echo "Please run 002_determine_chromosome_sizes.sh"
  exit 2
fi

if [[ ! -f "$AWK_SHIFT_FILE" ]]; then
  echo "[ERROR] ${AWK_SHIFT_FILE} not present."
  echo "Download from lab utils repository"
  exit 2
fi

# Initialization of array with keys and file paths
# Create an array of files
mapfile -t FILES < <(find "$HOME/data/preprocessing_test/align" -maxdepth 1 -type f -name "*_blFiltered.bam" | sort)

# Count unique sample types (excluding inputs)
number_of_sample_types=$(basename -a "${FILES[@]}" | cut -d_ -f1-2 | sort -u | wc -l)

# Check if any files were found
if [ ${#FILES[@]} -eq 0 ]; then
  echo "No *.bam files found in ~/data/preprocessing_test/align" >&2
  exit 1
fi

# === Load fragment samples for sample types ===
# Read log file paths into array
mapfile -t FRAG_LOGS < <(find "$OUTDIR/predictd/" -type f -name "*macs2.log")

# Check if any logs were found
if [[ ${#FRAG_LOGS[@]} -eq 0 ]]; then
  echo "[ERROR] Failed to find fragment size logs in $OUTDIR/predictd/" >&2
  exit 1
fi

# Validate number of log files
if [[ ${#FRAG_LOGS[@]} -ne $number_of_sample_types ]]; then
  echo "[ERROR] Expected $number_of_sample_types fragment logs, found ${#FRAG_LOGS[@]}" >&2
  echo "Found logs:" >&2
  printf '%s\n' "${FRAG_LOGS[@]}" >&2
  exit 2
fi

# Initialize fragments associative array
declare -A FRAGMENT_SIZES=()
for log in "${FRAG_LOGS[@]}"; do
  # Extract sample type from path
  #sample_name=$(basename "$(dirname "$log")" | sed 's/_predictd//')
  sample_name=$(basename "$log" | cut -d_ -f1-2 )

  # Extract fragment size
  frag_size=$(grep -oP 'predicted fragment length is \K\d+' "$log")

  echo "Log file: $log"
  echo "Sample type: $sample_name"
  echo "Fragment size: $frag_size"

  # Check if extraction succeeded
  if [[ -z "$frag_size" ]]; then
    echo "[ERROR] Failed to extract fragment size from $log" >&2
    exit 3
  fi

  # Validate fragment size
  if ! [[ "$frag_size" =~ ^[0-9]+$ ]] || (( frag_size <= 0 )); then
    echo "[ERROR] Invalid fragment size ($frag_size) in $log" >&2
    exit 4
  fi

  # Store in associative array
  FRAGMENT_SIZES[$sample_name]=$frag_size
done

# === Step 4: Read Shifting ===
for filepath in "${FILES[@]}";
do
  filename=$(basename "$filepath")
  sample_name=${filename%.bam}
  sample_type=$( echo "$sample_name" | cut -d_ -f1-2 )
  shifted_bam="$OUTDIR/align/${sample_name}_shifted.bam"
  frag_size=${FRAGMENT_SIZES[$sample_type]}
  shift_size=$((frag_size / 2))

  # Debug
  echo -e "\n=== Processing $sample_type ==="
  echo "Filepath: $filepath"
  echo "Filename: $filename"
  echo "Sample name: $sample_name"
  echo "Sample type: $sample_type"
  echo "Shifted bam: $shifted_bam"
  echo "Fragment size: $frag_size"
  echo "Shift amount: $shift_size bp"

  # Sanity checks
  [[ -n "$frag_size" ]] || { echo "Fragment size missing for $sample_type"; exit 1; }
  [[ "$shifted_bam" =~ \.bam$ ]] || { echo "Output must be .bam: $shifted_bam" >&2; exit 1; }

  # Check if output already exists and is valid
  if [[ -f "$shifted_bam" ]] && samtools quickcheck -v "$shifted_bam" &> /dev/null; then
    echo "Skipping existing shifted BAM: $shifted_bam"
    continue
  fi

  if [[ "$sample_type" == "input" ]]; then
    echo "Skipping input shifting..."
    continue
  fi

  #samtools view -h "$filepath" | \
  #  awk -v shift="$shift_size" -v chrom_file="$CHROM_SIZES_FILE" -v log_file="$OUTDIR/shift_stats.log" -f "$AWK_SHIFT_FILE" | \
  #  samtools sort -@ $(( THREADS / 2 )) | \
  #  samtools view -bS - > "$shifted_bam" || \
  #  {
  #    echo "Shifting failed for $filepath" >&2
  #    rm -f "$shifted_bam"
  #    return 1
  #  }
  ## Cleanup and validation
  ##rm "$chrom_sizes_file"
  #samtools index "$shifted_bam"
  #echo "Shift validation:"
  ##samtools idxstats "$shifted_bam"

  #if [[ -f "$shifted_bam" ]]; then
  #  echo "Shift Complete: $(basename "$shifted_bam")"
  #else
  #  echo "[WARNING] Shift failed for $shifted_bam."
  #fi
done # end bam shift for loop
