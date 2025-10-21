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
if [[ $# -lt $MIN_NUMBER_OF_ARGS ]] || [[ $# -gt $MAX_NUMBER_OF_ARGS ]]; then
  echo "Error: No argument provided." >&2
  show_usage
fi

# Check for help flag
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
  show_usage
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
EXPECTED_EXPERIMENT_ID_PATTERN=^[0-9]{8}Bel$ # Do not quote regular expression.
NCORES=$(nproc)
#FILETYPE_TO_KEEP="*.fastq"
#EXCLUDING_PATTERN="*unmapped*"

#============================== 
# Setup and preprocessing
#============================== 
EXPERIMENT_ID=${1%/} # Ensure argument does not have trailing slashes.
EXPERIMENT_DIR="$HOME/data/${EXPERIMENT_ID}"
FASTQ_DIRECTORY="$EXPERIMENT_DIR/fastq/"
DOCUMENTATION_DIR="$(dirname "$FASTQ_DIRECTORY")/documentation"
MANIFEST_FILE="$DOCUMENTATION_DIR/paired_reads_manifest.tsv"

#============================== 
# Error handling
#============================== 
echo "Running error handling..."
if [[ ! $EXPERIMENT_ID =~ $EXPECTED_EXPERIMENT_ID_PATTERN ]]; then
  echo "Error: EXPERIMENT_ID does not match expected pattern." >&2
  echo "Please adjust EXPERIMENT_ID accordingly." >&2
  echo "EXPERIMENT ID PATTERN: $EXPECTED_EXPERIMENT_ID_PATTERN" >&2
  echo "EXPERIMENT_ID: $EXPERIMENT_ID" >&2
  exit 1
fi

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
declare -A sample_id_map  # Use associative array for uniqueness
declare -A lane_map
declare -A pair_indicator_map

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

  # --- Extract sample ID ---
  sample_id="${parts[$SAMPLE_ID_START_IDX]}-${parts[$SAMPLE_ID_END_IDX]}"
  sample_id_map["$sample_id"]=1

  # --- Extract lane number ---
  lane_number="${parts[$LANE_IDX]}"
  lane_map["$lane_number"]=1

  # --- Extract read pair indicator ---
  read_indicator="${parts[$READ_PAIR_IDX]}"
  pair_indicator_map["$read_indicator"]=1

done

echo "Extracting unique metadata values..."
mapfile -t unique_sample_ids < <(printf '%s\n' "${!sample_id_map[@]}" | sort)
mapfile -t detected_lanes < <(printf '%s\n' "${!lane_map[@]}" | sort -n)
mapfile -t unique_pair_indicator < <(printf '%s\n' "${!pair_indicator_map[@]}" | sort)
EXPECTED_LANES_PER_SAMPLE=${#detected_lanes[@]}

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

echo "Processing ${#unique_sample_ids[@]} sample IDs"
echo "----------------"
echo "Unique ids:"
printf '%s\n' "${unique_sample_ids[@]}" | xargs -n6 | sed 's/^/  /' | column -t
echo "  Lanes found: ${detected_lanes[*]}"
echo "  Found read indicators: ${unique_pair_indicator[*]}"
echo "----------------"

#read -rp "Proceed with job submission? (y/n): " confirm
#confirm=$(echo "$confirm" | tr '[:upper:]' '[:lower:]')
#
#if [[ "$confirm" != "y" ]]; then
#  echo "Job submission cancelled"
#  exit 4
#fi

echo "Job confirmed. Proceed with consolidation..."

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
