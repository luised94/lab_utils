#!/usr/bin/env bash
################################################################################
# Call peaks for raw, deduplicated and the shifted bam
################################################################################
# Purpose: Generate peak calling files for preprocessed bam files
# Usage: Run as script.
# $./007_call_peaks_for_preprocessed_bam_files.sh
# DEPENDENCIES: bash, data from 250207Bel BMC experiment and scripts 001-005
# OUTPUT: A bunch of bed and peak files
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

# Genome data
GENOME_SIZE=12000000  # 1.2e7 in integer form

# Output config
OUTDIR="$HOME/data/preprocessing_test"
SUB_DIRS=("align" "predictd" "peaks" "coverage")
# Create directory structure
mkdir -p "${SUB_DIRS[@]/#/$OUTDIR/}"

# Initialization of array with keys and file paths
# Create an array of files
mapfile -t FILES < <(find ~/data/preprocessing_test/align -maxdepth 1 -type f -name "*.bam" | sort)
# Count unique sample types (excluding inputs)
#number_of_sample_types=$(basename -a "${FILES[@]}" | cut -d_ -f1-2 | grep -v "input" | sort -u | wc -l)

# Check if any files were found
if [ ${#FILES[@]} -eq 0 ]; then
  echo "No *.bam files found in ~/data/preprocessing_test/align" >&2
  exit 1
fi

# Find input control
INPUT_FILE_PATTERN="input_*_blFiltered.bam"
# Use find to get all matching files
mapfile -t matching_files < <(find "$OUTDIR/align/" -name "$INPUT_FILE_PATTERN" -type f)

# Check how many files were found
if [[ ${#matching_files[@]} -eq 0 ]]; then
    echo "[ERROR] No input control files matching pattern '$INPUT_FILE_PATTERN' found."
    exit 1
fi

if [[ ${#matching_files[@]} -gt 1 ]]; then
    echo "[ERROR] Multiple input control files found matching pattern '$INPUT_FILE_PATTERN':"
    printf "  %s\n" "${matching_files[@]}"
    echo "Please specify a more precise pattern."
    exit 1
fi

# Exactly one file found - success!
INPUT_CONTROL="${matching_files[0]}"
echo "Found input control file: $INPUT_CONTROL"

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
  basename=${filename%.bam}
  sample_type=$(echo "$basename" | cut -d'_' -f1)
  sample_name=$(echo "$basename" | cut -d'_' -f1-2)
  bam_type=$(echo "$basename" | cut -d'_' -f3 )
  fragment_size=${FRAGMENT_SIZES[$sample_name]}

  echo "---Sample information---"
  echo "  Filepath: $filepath"
  echo "  Filename: $filename"
  echo "  Sample_name: $sample_name"
  echo "  Sample type: $sample_type"
  echo "  Bam type: $bam_type"
  echo "  Fragment_size: $fragment_size"

  #if [[ "$sample_type" == "input" ]]; then
  #  echo "Skipping input sample peak calling..."
  #  continue
  #fi

  for peak_mode in "${!PEAK_MODES[@]}" ; do
    output_dir="$OUTDIR/peaks"

    # Always run without input control
    output_name="${basename}_${peak_mode}"
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
      #macs2 "${flags[@]}" > "${output_dir}/${output_name}.log" 2>&1 || {
      #  echo "[ERROR] MACS2 peak calling failed for $output_name" >&2
      #}
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
        #macs2 "${flags[@]}" >> "${output_dir}/${output_name}.log" 2>&1 || {
        #  echo "[ERROR] MACS2 peak calling failed for $output_name" >&2
        #}
      fi
    fi

  done # end of peak calling for loop
done # end of file for loop
# Create array of bigwig files
declare -a PEAK_FILES
mapfile -t PEAK_FILES < <(find "$OUTDIR/peaks" -name "*.xls" | sort | uniq )

# Validate array
if [[ ${#PEAK_FILES[@]} -eq 0 ]];
then
  echo "ERROR: No bigwig files found in $OUTDIR/peaks" >&2
  echo "Verify output directory or $0 " >&2
  #exit 3
fi

echo "Processed ${#FILES[@]} bam files..."
echo "Generated ${#PEAK_FILES[@]} peak files..."
echo "Peak calling done."
