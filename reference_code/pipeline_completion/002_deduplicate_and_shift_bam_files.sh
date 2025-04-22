#!/usr/bin/env bash
################################################################################
# Remove duplicates and shift reads for sample bam files
################################################################################
# Purpose: Generate bam files with duplicated reads removed and the bam files with shifted reads.
# Usage: Run as script.
# $./002_deduplicate_and_shift_bam_files.sh
# DEPENDENCIES: bash, data from 250207Bel BMC experiment
# OUTPUT: Duplicate files with renaming for easier identification and isolation from data directories
# NOTES: Updated the files being used in analysis, went from the data in the 241010Bel to the data in 250207Bel
# AUTHOR: LEMR
# DATE: 2025-04-17
# UPDATE: 2025-02-07
################################################################################
echo "Starting script"
#set -euo pipefail

# === Cluster Environment Setup ===
#if [[ "$(hostname)" != "luria" ]]; then
#    echo "Error: This script must be run on luria cluster" 1>&2
#    exit 1
#fi

# Load required modules
module load picard/2.18.26 java samtools
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
# Cluster config
THREADS=8
declare -A CHROM_SIZES=()

# File paths
# Store find results in a variable
# See ./lab_utils/reference_code/pipeline_completion/duplicate_files_for_testing.sh to see what samples should be initialized.
# Alternatively verify the directory with the find command.
FILES=$(find ~/data/preprocessing_test/align -maxdepth 1 -type f -name "*_raw.bam" | sort)
#mapfile -d '' FILES < <(find ~/data/preprocessing_test/align -maxdepth 1 -type f -name "*_raw.bam" -print0)

# Step 2: Check if any files were found
if [ -z "$FILES" ]; then
    echo "No *_raw.bam files found in ~/data/preprocessing_test/align" >&2
    exit 1
fi

# Process the files if found
declare -A SAMPLES
while IFS= read -r filepath; do
    filename=$(basename "$filepath")
    key=${filename%_raw.bam}
    SAMPLES["$key"]="$filepath"
done <<< "$FILES"

# Check if the SAMPLES array is not empty
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "SAMPLES array is empty!"
  echo "Check ~/data/preprocessing_test/align directory for '*_raw.bam' bam files."
else
  echo "SAMPLES array initialized successfully:"
  for key in "${!SAMPLES[@]}"; do
    echo "  Sample key: $key"
    echo "    Filepath: ${SAMPLES[$key]}"
  done
fi

#SELECTED_SAMPLE_KEYS=$(echo "${!SAMPLES[@]}" | sed 's/ /\n/g' | grep -v "input")
mapfile -t SELECTED_SAMPLE_KEYS < <(echo "${!SAMPLES[@]}" | tr ' ' '\n' | grep -E '^test|^reference')

# Reference data
GENOME_SIZE=12000000  # 1.2e7 in integer form
GENOME_FASTA="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C/SaccharomycescerevisiaeS288C_refgenome.fna"
if [[ ! -f $GENOME_FASTA ]]; then
  echo "$GENOME_FASTA does not exist."
  exit 1
fi

# Output config
OUTDIR="$HOME/data/preprocessing_test"
SUB_DIRS=("align" "predictd" "peaks" "coverage")

##############################################
# Initialization & Validation
##############################################
# Create directory structure
mkdir -p "${SUB_DIRS[@]/#/$OUTDIR/}"

# Display parameter summary
echo "=== Pipeline Configuration ==="
echo "Genome Size: $GENOME_SIZE"
echo "Output Directory: $OUTDIR"
echo "Threads: $THREADS"
echo "Samples:"
for stype in "${!SAMPLES[@]}"; do
  echo " - $stype: ${SAMPLES[$stype]}"
done
echo "=============================="

##############################################
# Preprocessing Workflow
##############################################
# === Step 1: Duplicate Removal ===
echo -e "\n=== Marking Duplicates ==="
for sample_type in "${!SAMPLES[@]}"; do
  raw_bam="${SAMPLES[$sample_type]}"
  deduplicated_bam="$OUTDIR/align/${sample_type}_deduped.bam"

  echo "Processing $sample_type..."
  echo "Input: $raw_bam"
  echo "Output: $deduplicated_bam"

  if [[ -f "$deduplicated_bam" ]]; then
    echo "Skipping existing: $deduplicated_bam"
  else
    java -jar $PICARD_JAR MarkDuplicates \
        I="$raw_bam" \
        O="$deduplicated_bam" \
        M="$OUTDIR/${sample_type}_dup_metrics.txt" \
        REMOVE_DUPLICATES=true
  fi
  # Validate output
  if [[ ! -s "$deduplicated_bam" ]]; then
    echo "ERROR: Failed to create $deduplicated_bam" >&2
    exit 2
  fi
  #echo "$sample_type Deduped Reads: "
  #samtools view -c "$output"
  #echo "Inspect first few reads"
  #samtools view "$output" | head
done # end deduplicate bam for loop

# === Step 1b: Index Deduplicated BAMs ===
echo -e "\n=== Indexing Deduplicated Files ==="
for sample_type in "${!SAMPLES[@]}"; do
  deduplicated_bam="$OUTDIR/align/${sample_type}_deduped.bam"
  index_file="${deduplicated_bam}.bai"
  echo "Indexing $sample_type deduped BAM..."

  if [[ -f "$index_file" ]]; then
    echo "Index exists: $index_file"
  else
    echo "Indexing: $deduplicated_bam"
    samtools index "$deduplicated_bam"
    exit_code=$?
    if [[ $exit_code -eq 0 ]] && [[ -f "$index_file" ]]; then
      echo "File $deduplicated_bam indexed succesfully"
    else
      echo "[ERROR] Indexing failed for ${deduplicated_bam} (exit $exit_code)" >&2
      exit 3
    fi
  fi
done # end indexing deduplicated bam file for loop

# === Step 2: Fragment Analysis ===
echo -e "\n=== Analyzing Fragment Sizes ==="
declare -A FRAGMENTS=()
for sample_type in "${SELECTED_SAMPLE_KEYS[@]}"; do
  echo "Processing $sample_type..."
  deduplicated_bam="$OUTDIR/align/${sample_type}_deduped.bam"
  outdir="$OUTDIR/predictd/${sample_type}_predictd"
  log_file="$outdir/macs2.log"

  mkdir -p "$outdir"

  # Only run MACS2 if log is missing/incomplete
  if [[ ! -f "$log_file" ]] || ! grep -q 'predicted fragment length is' "$log_file"; then
    echo "Running fragment size prediction..."
    macs2 predictd \
        --ifile "$deduplicated_bam" \
        --mfold 2 200 \
        --gsize "$GENOME_SIZE" \
        --outdir "$outdir" 2> "$log_file"
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

  FRAGMENTS[$sample_type]=$frag_size
  echo -e "Fragment Size Analysis Complete\n  Sample: $sample_type\n  Size: ${frag_size}bp\n  Log: ${log_file/$OUTDIR/\$OUTDIR}\n"
done # end fragment size analysis for loop

# === Step 3: Initializing chromosome lengths ===
# In main workflow, after variable initialization
echo -e "\n=== Chromosome Size Initialization ==="
CHROM_SIZES_FILE="$OUTDIR/chrom.sizes"
# Generate .fai if missing
if [[ ! -f "${GENOME_FASTA}.fai" ]]; then
  echo "Creating FASTA index for reference genome..."
  samtools faidx "$GENOME_FASTA" || {
    echo "[ERROR] Failed to index reference FASTA" >&2
    exit 10
  }
else
  echo "Genome fasta index present. Skipping..."
fi

# Generate chrom.sizes if missing
if [[ ! -f "$CHROM_SIZES_FILE" ]]; then
  echo -e "\nBuilding chromosome size file..."
  echo "Input FASTA: $(basename "$GENOME_FASTA")"
  echo "Output sizes: $CHROM_SIZES_FILE"

  awk '{print $1 "\t" $2}' "${GENOME_FASTA}.fai" > "$CHROM_SIZES_FILE" || {
    echo "[ERROR] Failed to create chrom.sizes" >&2
    exit 11
  }

  # Validate non-empty output
  [[ -s "$CHROM_SIZES_FILE" ]] || {
    echo "[ERROR] chrom.sizes file is empty" >&2
    exit 12
  }
else
  echo "Chrom sizes file present. Skipping..."
fi

# Load into associative array with debug output
echo -e "\nLoading chromosome sizes:"
while IFS=$'\t' read -r chrom size; do
  # Skip empty lines
  [[ -z "$chrom" ]] && continue

  # Trim additional fields after first tab
  chrom="${chrom%%[[:space:]]*}"

  echo " - ${chrom}: ${size}bp"
  CHROM_SIZES["$chrom"]="$size"
done < "$CHROM_SIZES_FILE"

echo -e "ChromDB status: Loaded ${#CHROM_SIZES[@]} chromosomes\n"
echo "======================================="

# Example chromosome checks
test_chr="chrI"
echo "Debug: ${test_chr} length = ${CHROM_SIZES[$test_chr]:-UNDEFINED}"

# === Step 4: Read Shifting ===
echo -e "\n=== Shifting Reads ==="
conda deactivate
echo "Deactivated macs2 conda environment"
if ! command -v macs2 &> /dev/null; then
  echo "MACS2 not found in environment. Deactivation confirmed"
fi

module load python/2.7.13 deeptools/3.0.1
echo "Activated python and deeptools"
# Apply to samples
for sample_type in "${SELECTED_SAMPLE_KEYS[@]}"; do
  deduplicated_bam="$OUTDIR/align/${sample_type}_deduped.bam"
  shifted_bam="$OUTDIR/align/${sample_type}_shifted.bam"
  frag_size=${FRAGMENTS[$sample_type]}
  shift_size=$((frag_size / 2))

  # Debug
  echo -e "\n=== Processing $sample_type ==="
  echo "Deduplicated Bam: $deduplicated_bam"
  echo "Fragment size: $frag_size"
  echo "Shift amount: $shift_size bp"
  echo "Chromosome count: ${#CHROM_SIZES[@]}"

  # Sanity checks
  [[ -n "$frag_size" ]] || { echo "Fragment size missing for $sample_type"; exit 1; }
  [[ "$deduplicated_bam" =~ \.bam$ ]] || { echo "deduplicated_bam must be .bam: $deduplicated_bam" >&2; exit 1; }
  [[ "$shifted_bam" =~ \.bam$ ]] || { echo "Output must be .bam: $shifted_bam" >&2; exit 1; }

  # Check if output already exists and is valid
  if [[ -f "$shifted_bam" ]] && samtools quickcheck -v "$shifted_bam" &> /dev/null; then
    echo "Skipping existing shifted BAM: $shifted_bam"
    continue
  fi
  samtools view -h "$deduplicated_bam" | \
    awk -v shift="$shift_size" -v chrom_file="$OUTDIR/chrom.sizes" -v log_file="$OUTDIR/shift_stats.log" -f "$HOME/lab_utils/core_scripts/shift_reads.awk" | \
    samtools sort -@ $THREADS | \
    samtools view -bS - > "$shifted_bam" || {
      echo "Shifting failed for $deduplicated_bam" >&2
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
    echo "Shift failed for $shifted_bam."
  fi
done # end bam shift for loop
