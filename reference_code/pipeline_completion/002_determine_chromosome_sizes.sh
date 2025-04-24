#!/usr/bin/env bash
################################################################################
# Determine chromosome sizes for shifting
################################################################################
# Purpose: Generate chromosome sizes for shifting, used by awk script
# Usage: Run as script.
# $./002_determine_chromosome_sizes.sh
# DEPENDENCIES: bash, fasta of s288c, see download script for specific version
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

# Genome data
GENOME_SIZE=12000000  # 1.2e7 in integer form
GENOME_FASTA="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C/SaccharomycescerevisiaeS288C_refgenome.fna"
if [[ ! -f $GENOME_FASTA ]]; then
  echo "$GENOME_FASTA does not exist."
  exit 1
fi

# Output config
OUTDIR="$HOME/data/preprocessing_test"
SUB_DIRS=("align" "predictd" "peaks" "coverage")

# Create directory structure
mkdir -p "${SUB_DIRS[@]/#/$OUTDIR/}"
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
# Reusable component start
echo -e "\nLoading chromosome sizes:"
declare -A CHROM_SIZES=()
while IFS=$'\t' read -r chrom size; do
  # Skip empty lines
  [[ -z "$chrom" ]] && continue

  # Trim additional fields after first tab
  chrom="${chrom%%[[:space:]]*}"

  echo " - ${chrom}: ${size}bp"
  CHROM_SIZES["$chrom"]="$size"
done < "$CHROM_SIZES_FILE"
# Reusable component end

# Example chromosome checks
echo -e "ChromDB status: Loaded ${#CHROM_SIZES[@]} chromosomes\n"
test_chr="chrI"
echo "Debug: ${test_chr} length = ${CHROM_SIZES[$test_chr]:-UNDEFINED}"
echo "======================================="
