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
OUTPUT_DIR="$OUTDIR/peaks"

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
#declare -A PEAK_MODES=(
  #['narrow_q']="--nomodel --call-summits --qvalue 0.01"
  #['broad_q']="--broad --broad-cutoff 0.1 --nomodel --qvalue 0.01"
  #['auto_q']="--fix-bimodal --call-summits --qvalue 0.01"
#  #['broad_p']="--broad --broad-cutoff 0.1 --nomodel --pvalue 1e-6"
#  #['narrow_p']="--nomodel --call-summits --pvalue 1e-6"
#  #['auto_p']="--fix-bimodal --call-summits --pvalue 1e-6"
#)

declare -a BASE_FLAGS=(
      callpeak
      --outdir "$OUTPUT_DIR"
      --gsize "$GENOME_SIZE"
      --SPMR
      --qvalue 0.01
)

declare -A INPUT_PARAMETERS=(
  ['withInput']="--control ${INPUT_CONTROL}"
  ['noInput']=""
)

declare -A PEAK_TYPES=(
  ['narrow']="--call-summits"
  ['broad']="--broad --broad-cutoff 0.1"
)

# Process the files if found
for filepath in "${FILES[@]}";
do

  filename=$(basename "$filepath")
  basename=${filename%.bam}
  sample_type=$(echo "$basename" | cut -d'_' -f1)
  sample_name=$(echo "$basename" | cut -d'_' -f1-2)
  bam_type=$(echo "$basename" | awk -F_ '{print $NF}' )

  echo "---Sample information---"
  echo "  Filepath: $filepath"
  echo "  Filename: $filename"
  echo "  Sample_name: $sample_name"
  echo "  Sample type: $sample_type"
  echo "  Bam type: $bam_type"


  # Assign mutually exclusive options based on bam types
  if [[ "$bam_type" == "shifted" ]]; then
    echo "Setting the flags for shifted bam file..."
    bam_type_flags=( "--nomodel" )
  else
    echo "Setting the flags for non-shifted bam file..."
    bam_type_flags=( "--fix-bimodal" )
  fi
  ##########

  for input_parameter in "${!INPUT_PARAMETERS[@]}";
  do

    prefix_with_input_parameter="${basename}_${input_parameter}"
    if [[ $input_parameter == withInput ]] && [[ -n "$INPUT_CONTROL" ]];
    then
      input_control_flag=("${INPUT_PARAMETERS[$input_parameter]}")
    fi

    for peak_type in "${!PEAK_TYPES[@]}" ;
    do

      final_output_name_prefix="${prefix_with_input_parameter}_${peak_type}"
      if [[ -f "${final_output_name_prefix}.narrowPeak" || -f "${final_output_name_prefix}.broadPeak" ]];
      then
        echo "Skipping ${final_output_name_prefix} file"
        continue
      fi

      IFS=' ' read -r -a mode_params <<< "${PEAK_TYPES[$peak_type]}"
      macs2_command_flags=(
        "${BASE_FLAGS[@]}"
        --name "$final_output_name_prefix"
        --treatment "${filepath}"
        "${bam_type_flags[@]}"
        "${input_control_flag[@]}"
        "${mode_params[@]}"
      )

      echo "  Input parameter: ${input_parameter}"
      echo "  Peak type: ${peak_type}"
      echo "  Prefix with input: ${prefix_with_input_parameter}"
      echo "  Output name: ${final_output_name_prefix}"
      echo -e "  COMMAND:\nmacs2 ${macs2_command_flags[*]}"
      #macs2 "${macs2_command_flags[@]}" > "${output_dir}/${basename}.log" 2>&1 || {
      #  echo "[ERROR] MACS2 peak calling failed for $output_name" >&2
      #}
    done

  done

done # end of peak calling for loop

# Create array of bigwig files
declare -a PEAK_FILES
mapfile -t PEAK_FILES < <(find "$OUTDIR/peaks" -name "*.xls" | sort | uniq )

# Validate array
if [[ ${#PEAK_FILES[@]} -eq 0 ]];
then
  echo "ERROR: No peak files found in $OUTDIR/peaks" >&2
  echo "Verify output directory or $0 " >&2
  #exit 3
fi

echo "Processed ${#FILES[@]} bam files..."
echo "Generated ${#PEAK_FILES[@]} peak files..."
echo "Peak calling done."
