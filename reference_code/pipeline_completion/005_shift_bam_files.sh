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
THREADS=8
CHROM_SIZES_FILE="$OUTDIR/chrom.sizes"
OUTDIR="$HOME/data/preprocessing_test"
SUB_DIRS=("align" "predictd" "peaks" "coverage")
# Create directory structure
mkdir -p "${SUB_DIRS[@]/#/$OUTDIR/}"

# Generate chrom.sizes if missing
if [[ ! -f "$CHROM_SIZES_FILE" ]]; then
  echo "[ERROR] ${CHROM_SIZES_FILE} not present."
  echo "Please run 002_determine_chromosome_sizes.sh"
else
  echo "Chromosome sizes present..."
fi

module load samtools

# Initialization of array with keys and file paths
# Create an array of files
mapfile -t FILES < <(find ~/data/preprocessing_test/align -maxdepth 1 -type f -name "*_blFiltered.bam" | sort)
# Count unique sample types (excluding inputs)
number_of_sample_types=$(basename -a "${FILES[@]}" | cut -d_ -f1-2 | grep -v "input" | sort -u | wc -l)

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

  # Check if extraction succeeded
  if [[ -z "$frag_size" ]]; then
    echo "[ERROR] Failed to extract fragment size from $log" >&2
    exit 3
  fi

  echo "Sample type: $sample_name"
  echo "Fragment size: $frag_size"

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
  sample_type=$( echo "$sample_name" | cut -d_ -f1 )
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

  samtools view -h "$filepath" | \
    awk -v shift="$shift_size" -v chrom_file="$OUTDIR/chrom.sizes" -v log_file="$OUTDIR/shift_stats.log" -f "$HOME/lab_utils/core_scripts/shift_reads.awk" | \
    samtools sort -@ $THREADS | \
    samtools view -bS - > "$shifted_bam" || {
      echo "Shifting failed for $filepath" >&2
      rm -f "$shifted_bam"
      return 1
    }
  # Cleanup and validation
  #rm "$chrom_sizes_file"
  samtools index "$shifted_bam"
  echo "Shift validation:"
  samtools idxstats "$shifted_bam"

  if [[ -f "$shifted_bam" ]]; then
    echo "Shift Complete: $(basename "$shifted_bam")"
  else
    echo "[WARNING] Shift failed for $shifted_bam."
  fi
done # end bam shift for loop
