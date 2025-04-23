#!/usr/bin/env bash
################################################################################
# Call peaks for raw, deduplicated and the shifted bam
################################################################################
# Purpose: Generate peak calling files for preprocessed bam files
# Usage: Run as script.
# $./003_deduplicate_and_shift_bam_files.sh
# DEPENDENCIES: bash, data from 250207Bel BMC experiment and scripts 001 and 002 from the pipeline script
# OUTPUT: Duplicate files with renaming for easier identification and isolation from data directories
# NOTES: Updated the files being used in analysis, went from the data in the 241010Bel to the data in 250207Bel
# AUTHOR: LEMR
# DATE: 2025-02-07
# UPDATE: 2025-04-22
################################################################################
#set -euo pipefail

# === Cluster Environment Setup ===
#if [[ "$(hostname)" != "luria" ]]; then
#    echo "Error: This script must be run on luria cluster" 1>&2
#    exit 1
#fi

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

# Initialization of array with keys and file paths
# Create an array of files
mapfile -t FILES < <(find ~/data/preprocessing_test/align -maxdepth 1 -type f -name "*.bam" | sort)
# Count unique sample types (excluding inputs)
number_of_sample_types=$(basename -a "${FILES[@]}" | cut -d_ -f1-2 | grep -v "input" | sort -u | wc -l)

# Check if any files were found
if [ ${#FILES[@]} -eq 0 ]; then
  echo "No *.bam files found in ~/data/preprocessing_test/align" >&2
  exit 1
fi

# Initialize Conda
CONDA_ROOT=~/miniforge3
if [[ -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
  . "$CONDA_ROOT/etc/profile.d/conda.sh"
else
  echo "ERROR: Conda not found at $CONDA_ROOT" >&2
  exit 1
fi

# Activate MACS2 environment
echo "Activating MACS2 environment..."
conda activate macs2_env 2> /dev/null

# Verify MACS2 installation
if ! command -v macs2 &> /dev/null; then
  echo "ERROR: MACS2 not installed. Install with:" >&2
  echo "Run: conda install -c bioconda macs2=2.2.7.1" >&2
  exit 2
fi

echo "MACS2 environment ready in $(which python)"
# === Load fragment samples for sample types ===
# Read log file paths into array
mapfile -t FRAG_LOGS < <(find "$OUTDIR/predictd/" -type f -name "macs2.log")

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
  sample_name=$(basename "$(dirname "$log")" | sed 's/_predictd//')

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

# Find input control
INPUT_CONTROL=""
for filepath in "${FILES[@]}"; do
  filename=$(basename "$filepath")
  if [[ "$filename" == input_*_deduped.bam ]]; then
    INPUT_CONTROL="$filepath"
    echo "Found input control file: $INPUT_CONTROL"
    break
  fi
done

# Check the input file was found
if [[ -z "$INPUT_CONTROL" ]]; then
  echo "[WARNING] No input control file found. Some analyses will be skipped."
fi

# Analysis variants with significance thresholds
declare -A PEAK_MODES=(
  ['narrow_p']="--nomodel --call-summits --pvalue 1e-6"
  ['narrow_q']="--nomodel --call-summits --qvalue 0.01"
  ['broad_p']="--broad --broad-cutoff 0.1 --nomodel --pvalue 1e-6"
  ['broad_q']="--broad --broad-cutoff 0.1 --nomodel --qvalue 0.01"
  ['auto_p']="--fix-bimodal --call-summits --pvalue 1e-6"
  ['auto_q']="--fix-bimodal --call-summits --qvalue 0.01"
)

# Process the files if found
for filepath in "${FILES[@]}"; do
  filename=$(basename "$filepath")
  key=${filename%.bam}
  sample_type=$(echo "$key" | cut -d'_' -f1)
  sample_name=$(echo "$key" | cut -d'_' -f1-2)
  bam_type=$(echo "$key" | cut -d'_' -f3 )
  fragment_size=${FRAGMENT_SIZES[$sample_name]}

  echo "---Sample information---"
  echo "  Filepath: $filepath"
  echo "  Filename: $filename"
  echo "  Sample_name: $sample_name"
  echo "  Sample type: $sample_type"
  echo "  Bam type: $bam_type"
  echo "  Fragment_size: $fragment_size"

  if [[ "$sample_type" == "input" ]]; then
    echo "Skipping input sample peak calling..."
    continue
  fi

  for peak_mode in "${!PEAK_MODES[@]}" ; do
    output_dir="$OUTDIR/peaks"

    # Always run without input control
    output_name="${key}_${peak_mode}"
    output_name_prefix="${output_name}_noInput"

    # Split the peak mode parameters into an array
    IFS=' ' read -r -a mode_params <<< "${PEAK_MODES[$peak_mode]}"

    # Core parameters used in all calls
    flags=(
      callpeak
      --treatment "${filepath}"
      --name "${output_name_prefix}"
      --outdir "$output_dir"
      --gsize "$GENOME_SIZE"
      --SPMR
    )

    # Add the peak mode-specific parameters
    flags+=("${mode_params[@]}")

    # Print the command that would be executed
    echo "  ---Peak mode information---"
    echo "    Output name: $output_name"
    echo -e "    COMMAND No input:\nmacs2 ${flags[*]}"

    if [[ -f "${output_name_prefix}.narrowPeak" || -f "${output_name_prefix}.broadPeak" ]]; then
      echo "  Skipping existing output: $(basename "$output_name_prefix")*Peak"
    else
      echo "Executing macs2 command for ${output_name_prefix}..."
      # Execute the command
      macs2 "${flags[@]}" > "${output_dir}/${output_name}.log" 2>&1 || {
        echo "[ERROR] MACS2 peak calling failed for $output_name" >&2
      }
    fi

    # === Second run WITH input control (if available) ===
    # Create the base flags array for input control run
    output_name_prefix="${output_name}_withInput"
    if [[ -f "${output_name_prefix}.narrowPeak" || -f "${output_name_prefix}.broadPeak" ]]; then
      echo "  Skipping existing output: $(basename "$output_name_prefix")*Peak"
    else
      echo "If input available, executing macs2 command for ${output_name_prefix} ..."
      if [[ -n "$INPUT_CONTROL" ]]; then
        flags=(
          callpeak
          --treatment "$filepath"
          --control "$INPUT_CONTROL"
          --name "${output_name_prefix}"
          --outdir "$output_dir"
          --gsize "$GENOME_SIZE"
          --SPMR
        )

        # Add the peak mode-specific parameters
        flags+=("${mode_params[@]}")

        # Print the command that would be executed
        echo -e "    COMMAND with input:\nmacs2 ${flags[*]}"

        # Execute the command
        macs2 "${flags[@]}" >> "${output_dir}/${output_name}.log" 2>&1 || {
          echo "[ERROR] MACS2 peak calling failed for $output_name" >&2
        }
      fi
    fi

  done # end of peak calling for loop
done # end of file for loop
echo "Peak calling done."
