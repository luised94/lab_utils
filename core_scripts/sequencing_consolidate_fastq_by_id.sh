#!/bin/bash
################################################################################
# Consolidate fastq files from different lanes.
# Author: Luis | Date: 2025-10-20 | Version: 2.0.0
################################################################################
# PURPOSE:
#   Consolidate fastq files from different lanes into a single large fastq file.
# USAGE:
#   From the command line
#   $ srun sequencing_consolidate_fastq_by_id.sh <EXPERIMENT_ID> [--active-run]
# DEPENDENCIES:
#   bash 4.2
#   Assumes fastq files are in the fastq directory and everything is removed. (cleanup script.)
# OUTPUTS:
#   Consolidated fastq files.
################################################################################
#============================== 
# Usage and help
#============================== 
show_usage() {
  cat << EOF
Usage: srun $(basename "$0") <fastq_directory> [-v]

Description:
  Consolidate fastq files from different lanes of the sequencer into a single file. Paired end files are consolidated in lane order.

Arguments:
  EXPERIMENT_ID    Experiment id of experiment (from BMC submission.)
                   (e.g., 250930Bel)

Options:
  -h, --help        Show this help message
  --active-run      Execute rsync command

Output:
  Directory with fastq files downloaded to the local directly in initial directory structure.

Example:
  $(basename "$0") 250930Bel
  $(basename "$0") 250930Bel --active-run
  $(basename "$0") -h

EOF
  exit 0
}

#============================== 
# Argument error handling
#============================== 
# Check for arguments
echo "Handling arguments..."
MIN_NUMBER_OF_ARGS=1
MAX_NUMBER_OF_ARGS=2
EXPECTED_EXPERIMENT_ID_PATTERN=^[0-9]{6}Bel$ # Do not quote regular expression.
# Check for help flag
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
  show_usage
fi

# Verify number of arguments
if [[ $# -lt $MIN_NUMBER_OF_ARGS ]]; then
    echo "Error: Missing required argument EXPERIMENT_ID." >&2
    show_usage
    exit 1
:w
    echo "Error: Too many arguments provided ($#)." >&2
    show_usage
    exit 1
fi

# Handle first argument: Remove trailing slash and validate pattern
EXPERIMENT_ID=${1%/} # Ensure argument does not have trailing slashes.
echo "Running error handling..."
if [[ ! $EXPERIMENT_ID =~ $EXPECTED_EXPERIMENT_ID_PATTERN ]]; then
  echo "Error: EXPERIMENT_ID does not match expected pattern." >&2
  echo "Please adjust EXPERIMENT_ID accordingly." >&2
  echo "EXPERIMENT ID PATTERN: $EXPECTED_EXPERIMENT_ID_PATTERN" >&2
  echo "EXPERIMENT_ID: $EXPERIMENT_ID" >&2
  exit 1
fi

# Handle second argument: Set dry-run mode.
DRY_RUN=true
if [[ $# -eq $MAX_NUMBER_OF_ARGS ]]; then
    if [[ "$2" != "--active-run" ]]; then
        echo "Error: Unknown option '$2'" >&2
        echo "Use --active-run to perform actual sync." >&2
        exit 1
    fi
    DRY_RUN=false
fi

# Ensure script is run inside a Slurm allocation
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    cat >&2 <<EOF
Error: This script must be run within a Slurm job.

To run interactively:
    srun $0 <EXPERIMENT_ID> [--active-run]

To submit as a batch job:
    echo "$0 <EXPERIMENT_ID>" [--active-run] | sbatch

EOF
    exit 1
fi

#============================== 
# Configuration
#============================== 
echo "Setting configuration..."
# Filename parsing indices (split on _ and -)
# Example: 250930Bel_D25-12496-2_1_sequence.fastq
# Parts: [250930Bel, D25, 12496, 2, 1, sequence.fastq]
SAMPLE_ID_START_IDX=1   # "D25"
SAMPLE_ID_END_IDX=2     # "12496"
LANE_IDX=3              # "2"
READ_PAIR_IDX=4         # "1" or "2"
NCORES=$(nproc)
OUTPUT_PREFIX="consolidated"
OUTPUT_SUFFIX=".fastq"
#FILETYPE_TO_KEEP="*.fastq"
#EXCLUDING_PATTERN="*unmapped*"

#============================== 
# Setup and preprocessing
#============================== 
EXPERIMENT_DIR="$HOME/data/${EXPERIMENT_ID}"
FASTQ_DIRECTORY="$EXPERIMENT_DIR/fastq/"
DOCUMENTATION_DIR="$(dirname "$FASTQ_DIRECTORY")/documentation"
MANIFEST_FILE="$DOCUMENTATION_DIR/paired_reads_manifest.tsv"

#============================== 
# Error handling
#============================== 

if [[ ! -d "$FASTQ_DIRECTORY" ]]; then
  echo "Error: FASTQ_DIRECTORY does not exist. Please verify experiment id." >&2
  echo "FASTQ_DIRECTORY: $FASTQ_DIRECTORY" >&2
  echo "Run $0 -h for additional help." >&2
  exit 1

fi

if [[ ! -d "$DOCUMENTATION_DIR" ]]; then
  echo "Error: DOCUMENTATION_DIR does not exist. Please verify experiment id." >&2
  echo "DOCUMENTATION_DIR: $DOCUMENTATION_DIR" >&2
  echo "Run $0 -h for additional help." >&2
  exit 1

fi

# 1. Check for any subdirectory (mindepth 1)
# Use quit when you need to know if there is just more thna zero.
if [[ -n "$(find "$FASTQ_DIRECTORY" -mindepth 1 -type d -print -quit 2>/dev/null)" ]]; then
    echo "ERROR: Subdirectories still exist in: $FASTQ_DIRECTORY" >&2
    exit 1
fi

# 2. Check for non-FASTQ files in top level
if [[ -n "$(find "$FASTQ_DIRECTORY" -maxdepth 1 -type f ! -name "*.fastq" -print -quit 2>/dev/null)" ]]; then
    echo "ERROR: Non-FASTQ files found in: $FASTQ_DIRECTORY" >&2
    exit 1
fi

# 3. Check that at least one FASTQ file exists in top level
if [[ -z "$(find "$FASTQ_DIRECTORY" -maxdepth 1 -type f -name "*.fastq" -print -quit 2>/dev/null)" ]]; then
    echo "ERROR: No FASTQ files found in: $FASTQ_DIRECTORY" >&2
    exit 1
fi

echo "Confirmed directory has been cleaned from other bmc files."
# ############################################
# Main logic
# ############################################
echo "Using FASTQ directory: ${FASTQ_DIRECTORY}"
echo "Manifest file: $MANIFEST_FILE"

# Validate manifest exists
if [[ -f "$MANIFEST_FILE" ]]; then
    echo "Warning: Manifest file already exists."

fi

# ============================================
# Discover all fastq files
# ============================================
echo "Looking for fastq files in: $FASTQ_DIRECTORY"

mapfile -t all_fastq_files < <(find "$FASTQ_DIRECTORY" -type f -name "*.fastq")
echo "Total fastq files found: ${#all_fastq_files[@]}"

# Check if any files found
if [[ ${#all_fastq_files[@]} -eq 0 ]]; then
  echo "ERROR: No .fastq files found in directory"
  exit 1
fi

# ============================================
# Extract sample IDs, lanes, and read types in one pass
# ============================================
echo "Analyzing FASTQ files..."
# For detection (uniqueness)
declare -A sample_id_map
declare -A lane_map
declare -A pair_indicator_map

# For consolidation (file grouping)
declare -A sample_r1_files
declare -A sample_r2_files
declare -A sample_se_files

# Loop overall files and extract based on indices.
for fastq_file in "${all_fastq_files[@]}"; do
  filename=$(basename "$fastq_file")

  if [[ "$VERBOSE" == true ]]; then
    echo "Filename: $filename"
  fi

  if [[ "$filename" =~ unmapped ]]; then
    if [[ "$VERBOSE" == true ]]; then
      echo "Skipping unmapped fastq file"
    fi
    continue
  fi

  # Split on _ and - to get components
  IFS='_-' read -ra parts <<< "$filename"

  # --- Extract the metadata ---
  sample_id="${parts[$SAMPLE_ID_START_IDX]}-${parts[$SAMPLE_ID_END_IDX]}"
  lane_number="${parts[$LANE_IDX]}"
  read_indicator="${parts[$READ_PAIR_IDX]}"

  # --- Add the values to the associative array ---
  sample_id_map["$sample_id"]=1
  lane_map["$lane_number"]=1
  pair_indicator_map["$read_indicator"]=1

  # Store for consolidation (grouped by sample, with lane prefix for sorting)
  if [[ "$read_indicator" == "1" ]]; then
    sample_r1_files["$sample_id"]+="$lane_number:$fastq_file"$'\n'

  elif [[ "$read_indicator" == "2" ]]; then
    sample_r2_files["$sample_id"]+="$lane_number:$fastq_file"$'\n'

  elif [[ "$read_indicator" == "NA" ]]; then
    sample_se_files["$sample_id"]+="$lane_number:$fastq_file"$'\n'

  else
    echo "WARNING: Unknown read indicator '$read_indicator' in file: $filename"

  fi

done

echo "Extracting unique metadata values..."
mapfile -t unique_sample_ids < <(printf '%s\n' "${!sample_id_map[@]}" | sort)
mapfile -t detected_lanes < <(printf '%s\n' "${!lane_map[@]}" | sort -n)
mapfile -t unique_pair_indicator < <(printf '%s\n' "${!pair_indicator_map[@]}" | sort)
LANES_PER_SAMPLE=${#detected_lanes[@]}

# ============================================
# Error handling for metadata extraction
# ============================================
if [[ ${#unique_sample_ids[@]} -eq 0 ]]; then
  echo "ERROR: No valid sample IDs extracted from filenames" >&2
  echo "Ensure the unique ids logic is correct and no updates have occured to names." >&2
  exit 1

fi

if [[ ${#detected_lanes[@]} -eq 0 ]]; then
  echo "ERROR: No lanes detected after extraction."
  exit 1

fi

if [[ ${#unique_pair_indicator[@]} -eq 0 ]]; then
  echo "ERROR: No lanes detected after extraction."
  exit 1

fi

# Validate lane numbers are reasonable (1-4 typical)
for lane in "${detected_lanes[@]}"; do
  if [[ ! "$lane" =~ ^[1-4]$ ]]; then
    echo "WARNING: Unusual lane number detected: $lane"
  fi
done

# Determine read type and validate
has_paired_indicators=false
has_single_indicator=false
has_invalid=false

for indicator in "${unique_pair_indicator[@]}"; do
  if [[ "$indicator" == "1" || "$indicator" == "2" ]]; then
    has_paired_indicators=true
  elif [[ "$indicator" == "NA" ]]; then
    has_single_indicator=true
  else
    has_invalid=true
    echo "  ERROR: Invalid read indicator found: '$indicator'"
  fi
done

# Check for errors
if [[ "$has_invalid" == true ]]; then
  echo "  ERROR: Unexpected read indicators found"
  echo "  Expected: '1', '2' (paired-end) or 'NA' (single-end)"
  exit 1
fi

if [[ "$has_paired_indicators" == true && "$has_single_indicator" == true ]]; then
  echo "  ERROR: Mixed read types detected (both paired and single-end)"
  echo "  Found indicators: ${unique_pair_indicator[*]}"
  exit 1
fi

# Set configuration
if [[ "$has_paired_indicators" == true ]]; then
  IS_PAIRED_END=true
  EXPECTED_READ_PAIRS_PER_LANE=2
elif [[ "$has_single_indicator" == true ]]; then
  IS_PAIRED_END=false
  EXPECTED_READ_PAIRS_PER_LANE=1
else
  echo "  ERROR: No valid read indicators found"
  exit 1
fi

# Check read type detection succeeded
if [[ -z "$IS_PAIRED_END" ]]; then
  echo "ERROR: Failed to detect read type"
  exit 1
fi

# Expected number of samples per file.
EXPECTED_FILES_PER_SAMPLE=$((LANES_PER_SAMPLE * EXPECTED_READ_PAIRS_PER_LANE))

echo "Processing ${#unique_sample_ids[@]} sample IDs"
echo "----------------"
echo "Unique ids:"
printf '%s\n' "${unique_sample_ids[@]}" | xargs -n6 | sed 's/^/  /' | column -t
echo "Lanes found: ${detected_lanes[*]}"
echo "Found read indicators: ${unique_pair_indicator[*]}"
echo "IS_PAIRED_END: $IS_PAIRED_END"
echo "Expected files per sample: $EXPECTED_FILES_PER_SAMPLE"
echo "----------------"

# ============================================
# Consolidation planning and execution
# ============================================
echo ""
echo "============================================"
echo "Starting consolidation process..."
echo "============================================"
echo ""

# Process each sample
samples_processed=0
samples_skipped=0

for sample_id in "${unique_sample_ids[@]}"; do
  echo "--------------------"
  echo "Sample: $sample_id"
  
  if [[ "$IS_PAIRED_END" == true ]]; then
    # Extract and sort R1 files
    r1_files_list=$(printf "%s" "${sample_r1_files[$sample_id]}" | sort -t: -k1,1n | cut -d: -f2-)
    # Extract and sort R2 files
    r2_files_list=$(printf "%s" "${sample_r2_files[$sample_id]}" | sort -t: -k1,1n | cut -d: -f2-)
    
    # Convert to arrays for counting
    mapfile -t r1_files_array <<< "$r1_files_list"
    mapfile -t r2_files_array <<< "$r2_files_list"
    
    echo "  Read type: Paired-end"
    echo "  R1 files found: ${#r1_files_array[@]}"
    echo "  R2 files found: ${#r2_files_array[@]}"
    
    # Display files in lane order
    echo "  R1 files (in lane order):"
    printf '    %s\n' "${r1_files_array[@]}"
    echo "  R2 files (in lane order):"
    printf '    %s\n' "${r2_files_array[@]}"
    
    # Define output filenames
    output_r1="${FASTQ_DIRECTORY}${OUTPUT_PREFIX}_${sample_id}_R1${OUTPUT_SUFFIX}"
    output_r2="${FASTQ_DIRECTORY}${OUTPUT_PREFIX}_${sample_id}_R2${OUTPUT_SUFFIX}"
    
    echo "  Would create:"
    echo "    $output_r1"
    echo "    $output_r2"
    
  else
    # Extract and sort single-end files
    se_files_list=$(printf "%s" "${sample_se_files[$sample_id]}" | sort -t: -k1,1n | cut -d: -f2-)
    
    # Convert to array for counting
    mapfile -t se_files_array <<< "$se_files_list"
    
    echo "  Read type: Single-end"
    echo "  Files found: ${#se_files_array[@]}"
    
    # Display files in lane order
    echo "  Files (in lane order):"
    printf '    %s\n' "${se_files_array[@]}"
    
    # Define output filename
    output_se="${FASTQ_DIRECTORY}${OUTPUT_PREFIX}_${sample_id}_NA${OUTPUT_SUFFIX}"
    
    echo "  Would create:"
    echo "    $output_se"
  fi
  
  echo "--------------------"
  echo ""
done

echo "============================================"
echo "Consolidation planning complete"
echo "Samples to process: ${#unique_sample_ids[@]}"
echo "============================================"

echo "Consolidation complete..."
#read -rp "Proceed with job submission? (y/n): " confirm
#confirm=$(echo "$confirm" | tr '[:upper:]' '[:lower:]')
#
#if [[ "$confirm" != "y" ]]; then
#  echo "Job submission cancelled"
#  exit 4
#fi

# Process each unique ID
#for id in "${unique_sample_ids[@]}"; do
#  echo "--------------------"
#  echo "Processing ID: $id"
#  # Find files using array and glob pattern
#  files=( *"${id}"*_sequence.fastq )
#
#  if [ ${#files[@]} -eq 0 ]; then
#    echo "ERROR: No files found for ID: $id"
#    continue
#  fi
#
#  if [ ! ${#files[@]} -eq 2 ]; then
#    echo -e "[ERROR]: Found ${#files[@]} for $id\nExpected 2"
#    continue
#  fi
#
#  # Validate all files before processing
#  for file in "${files[@]}"; do
#    if ! [ -f "$file" ]; then
#      echo "Error: $file not found"
#      exit 1
#    fi
#    if ! [[ "$file" =~ \.fastq$ ]]; then
#      echo "Error: $file is not a FASTQ file"
#      exit 1
#    fi
#    echo "Validated: $file"
#  done
#
#  output_file="consolidated_${id}_sequence.fastq"
#  tmp_file="${output_file}.tmp"
#
#  if [ -f "$output_file" ]; then
#    echo "$output_file already exists. Skipping..."
#    continue
#  fi
#
#  # Consolidate files with atomic write
#  if cat -- "${files[@]}" > "$tmp_file"; then
#    mv "$tmp_file" "$output_file"
#    echo "Successfully created $output_file"
#
#    # Verify file content
#    if [ -s "$output_file" ]; then
#      # Accurate size calculation
#      original_size=$(wc -c "${files[@]}" | awk '/total/ {print $1}')
#      new_size=$(wc -c < "$output_file")
#
#      echo "Original files total size: $original_size bytes"
#      echo "New file size: $new_size bytes"
#
#      # Calculate original data checksum
#      orig_checksum=$(cat "${files[@]}" | md5sum | cut -d' ' -f1)
#      new_checksum=$(md5sum "$output_file" | cut -d' ' -f1)
#
#      if [[ "$orig_checksum" != "$new_checksum" ]]; then
#        echo "Error: Consolidated file checksum mismatch!" >&2
#        echo "Expected: $orig_checksum" >&2
#        echo "Actual:   $new_checksum" >&2
#        exit 5
#      fi
#
#      # Only remove if verification passed
#      echo "Removing original files..."
#      rm -f -- "${files[@]}"
#      echo "Original files removed"
#    else
#      echo "Error: Consolidated file is empty"
#      rm -f "$output_file"
#      exit 5
#    fi
#  else
#    echo "Error during consolidation"
#    rm -f "$tmp_file" "$output_file"
#    exit 5
#  fi
#  echo "--------------------"
#done
#
#echo "Consolidation completed successfully in ${target_dir}"
## Return to original directory
#cd "$original_dir"
