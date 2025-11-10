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
fi

if [[ $# -gt $MAX_NUMBER_OF_ARGS ]]; then
    echo "Error: Too many arguments provided." >&2
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

# Filename parsing configuration
FILENAME_DELIMITERS='_-'  # Delimiters used to split filename components
OUTPUT_PREFIX="consolidated"
OUTPUT_SUFFIX=".fastq"
#NCORES=$(nproc)

# Manifest output configuration
MANIFEST_FILENAME="consolidated_reads_manifest.tsv"

#============================== 
# Setup and preprocessing
#============================== 
EXPERIMENT_DIR="$HOME/data/${EXPERIMENT_ID}"
FASTQ_DIRECTORY="$EXPERIMENT_DIR/fastq/"
DOCUMENTATION_DIR="${EXPERIMENT_DIR}/documentation"
MANIFEST_FILEPATH="$DOCUMENTATION_DIR/$MANIFEST_FILENAME"

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

# ============================================
# Check if consolidation has already occurred
# ============================================
consolidated_count=$(find "$FASTQ_DIRECTORY" -maxdepth 1 -name "${OUTPUT_PREFIX}_*${OUTPUT_SUFFIX}" 2>/dev/null | wc -l)

if [[ $consolidated_count -gt 0 ]]; then
  echo "INFO: Found $consolidated_count consolidated files"
  echo "INFO: Consolidation appears complete"

  # Check if manifest exists
  MANIFEST_FILEPATH="$DOCUMENTATION_DIR/$MANIFEST_FILENAME"

  if [[ -f "$MANIFEST_FILEPATH" ]]; then
    echo "INFO: Manifest file already exists: $MANIFEST_FILEPATH"
    echo "INFO: Nothing to do - exiting"
    exit 0

  else
    echo "WARNING: Manifest file not found: $MANIFEST_FILEPATH"

    if [[ "$DRY_RUN" == true ]]; then
      echo "INFO: Would regenerate manifest (use --active-run to create it)"
      exit 0

    fi

    echo "INFO: Regenerating manifest from consolidated files..."

    # Detect read type from consolidated filenames
    r1_count=$(find "$FASTQ_DIRECTORY" -maxdepth 1 -name "${OUTPUT_PREFIX}_*_R1${OUTPUT_SUFFIX}" 2>/dev/null | wc -l)
    r2_count=$(find "$FASTQ_DIRECTORY" -maxdepth 1 -name "${OUTPUT_PREFIX}_*_R2${OUTPUT_SUFFIX}" 2>/dev/null | wc -l)
    na_count=$(find "$FASTQ_DIRECTORY" -maxdepth 1 -name "${OUTPUT_PREFIX}_*_NA${OUTPUT_SUFFIX}" 2>/dev/null | wc -l)

    if [[ $r1_count -gt 0 && $r2_count -gt 0 ]]; then
      echo "Detected paired-end data (R1: $r1_count, R2: $r2_count)"

      # Create paired-end manifest
      echo -e "sample_id\tread1_path\tread2_path" > "$MANIFEST_FILEPATH"

      # Find all R1 files and extract sample IDs
      for r1_file in "$FASTQ_DIRECTORY"/${OUTPUT_PREFIX}_*_R1${OUTPUT_SUFFIX}; do
        # Extract sample_id from filename: consolidated_SAMPLEID_R1.fastq
        filename=$(basename "$r1_file")
        sample_id=$(echo "$filename" | sed "s/^${OUTPUT_PREFIX}_\(.*\)_R1${OUTPUT_SUFFIX}$/\1/")

        r2_file="${FASTQ_DIRECTORY}${OUTPUT_PREFIX}_${sample_id}_R2${OUTPUT_SUFFIX}"

        if [[ -f "$r2_file" ]]; then
          echo -e "${sample_id}\t${r1_file}\t${r2_file}" >> "$MANIFEST_FILEPATH"

        else
          echo "WARNING: Missing R2 file for sample: $sample_id"

        fi

      done

    elif [[ $na_count -gt 0 ]]; then
      echo "Detected single-end data (files: $na_count)"

      # Create single-end manifest
      echo -e "sample_id\tread_path" > "$MANIFEST_FILEPATH"

      # Find all NA files and extract sample IDs
      for na_file in "$FASTQ_DIRECTORY"/${OUTPUT_PREFIX}_*_NA${OUTPUT_SUFFIX}; do
        # Extract sample_id from filename: consolidated_SAMPLEID_NA.fastq
        filename=$(basename "$na_file")
        sample_id=$(echo "$filename" | sed "s/^${OUTPUT_PREFIX}_\(.*\)_NA${OUTPUT_SUFFIX}$/\1/")

        echo -e "${sample_id}\t${na_file}" >> "$MANIFEST_FILEPATH"
      done

    else
      echo "ERROR: Could not determine read type from consolidated files"
      exit 1

    fi

    manifest_entries=$(tail -n +2 "$MANIFEST_FILEPATH" | wc -l)
    echo "[OK] Manifest regenerated: $MANIFEST_FILEPATH"
    echo "  Entries: $manifest_entries"
    exit 0

  fi

fi

# ############################################
# Main logic
# ############################################
echo "Using FASTQ directory: ${FASTQ_DIRECTORY}"
echo "Manifest file: $MANIFEST_FILEPATH"

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

  # --- Skip unmapped and consolidated files if present ---
  if [[ "$filename" =~ unmapped ]]; then
    if [[ "$VERBOSE" == true ]]; then
      echo "Skipping unmapped fastq file"

    fi
    continue

  fi

  if [[ "$filename" =~ ^${OUTPUT_PREFIX}_ ]]; then
    continue

  fi

  # Split on _ and - to get components
  IFS="$FILENAME_DELIMITERS" read -ra parts <<< "$filename"

  # Inside your loop, after splitting the filename into 'parts'
  num_parts=${#parts[@]}

  # Check the number of parts to decide the structure
  if [[ $num_parts -eq 6 ]]; then
      # Format: 250930Bel_D25-12496-2_1_sequence.fastq
      # Parts: [250930Bel, D25, 12496, 2, 1, sequence.fastq]
      SAMPLE_ID_START_IDX=1
      SAMPLE_ID_END_IDX=2
      LANE_IDX=3
      READ_PAIR_IDX=4
  elif [[ $num_parts -eq 5 ]]; then
      # Format: 250930Bel_D25-12496_1_sequence.fastq
      # Parts: [250930Bel, D25, 12496, 1, sequence.fastq]
      SAMPLE_ID_START_IDX=1
      SAMPLE_ID_END_IDX=2
      LANE_IDX=-1 # Use an invalid index to indicate no lane
      READ_PAIR_IDX=3
  else
      echo "ERROR: Unexpected number of parts in filename: $filename"
      continue # Skip to the next file
  fi


  # --- Extract sample ID ---
  sample_id="${parts[$SAMPLE_ID_START_IDX]}-${parts[$SAMPLE_ID_END_IDX]}"
  sample_id_map["$sample_id"]=1

  # --- Extract lane number (conditionally) ---
  if [[ $LANE_IDX -ne -1 ]]; then
      lane_number="${parts[$LANE_IDX]}"
  else
      # If there's no lane, assign "NA" to both the variable and the map
      lane_number="NA"
  fi
  lane_map["$lane_number"]=1 # Add to map after setting the variable

  # --- Extract read pair indicator ---
  read_indicator="${parts[$READ_PAIR_IDX]}"
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

# Determine the filename convention based on detected lanes
LANE_FORMAT=""
if [[ ${#detected_lanes[@]} -eq 1 && "${detected_lanes[0]}" == "NA" ]]; then
    # This is a valid case: only files with no lane identifiers were found.
    echo "INFO: No lane identifiers were found in filenames. Setting format to 'none'."
    LANE_FORMAT="none"

else
    # This block executes if actual lane numbers were found.
    echo "INFO: Detected lanes: ${detected_lanes[*]}. Setting format to 'numeric'."
    LANE_FORMAT="numeric"

fi

for lane in "${detected_lanes[@]}"; do
  if [[ "$lane" == "NA" ]]; then
    continue

  fi

  if [[ ! "$lane" =~ ^[1-4]$ ]]; then
    echo "WARNING: Unusual lane number detected: $lane"

  fi

done

if [[ ${#unique_pair_indicator[@]} -eq 0 ]]; then
  echo "ERROR: No read pairs detected after extraction."
  exit 1

fi

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

  # --- Process a paired end read files ---
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

    # Validate file counts match expected lanes
    if [[ ${#r1_files_array[@]} -ne $LANES_PER_SAMPLE ]]; then
      echo "  ERROR: R1 file count (${#r1_files_array[@]}) does not match expected lanes ($LANES_PER_SAMPLE)"
      ((samples_skipped++))
      continue

    fi

    if [[ ${#r2_files_array[@]} -ne $LANES_PER_SAMPLE ]]; then
      echo "  ERROR: R2 file count (${#r2_files_array[@]}) does not match expected lanes ($LANES_PER_SAMPLE)"
      ((samples_skipped++))
      continue

    fi

    # Validate R1 and R2 counts match
    if [[ ${#r1_files_array[@]} -ne ${#r2_files_array[@]} ]]; then
      echo "  ERROR: R1 count (${#r1_files_array[@]}) does not match R2 count (${#r2_files_array[@]})"
      ((samples_skipped++))
      continue

    fi

    # Define output filenames
    output_r1="${FASTQ_DIRECTORY}${OUTPUT_PREFIX}_${sample_id}_R1${OUTPUT_SUFFIX}"
    output_r2="${FASTQ_DIRECTORY}${OUTPUT_PREFIX}_${sample_id}_R2${OUTPUT_SUFFIX}"

    # Check if output already exists
    if [[ -f "$output_r1" && -f "$output_r2" ]]; then
      echo "  Output files already exist - skipping"
      echo "    $output_r1"
      echo "    $output_r2"
      ((samples_skipped++))
      continue

    fi

    # Display files in lane order
    echo "  R1 files (in lane order):"
    printf '    %s\n' "${r1_files_array[@]}"
    echo "  R2 files (in lane order):"
    printf '    %s\n' "${r2_files_array[@]}"

    echo "  Will create:"
    echo "    $output_r1"
    echo "    $output_r2"

    # Only consolidate if not in dry-run mode
    if [[ "$DRY_RUN" == false ]]; then
      # Consolidate and verify R1 files
      echo "  Consolidating and verifying R1 files..."
      tmp_r1="${output_r1}.tmp"

      # Use tee to write to the file and pipe to md5sum simultaneously
      # This captures the checksum of the stream being written.
      stream_checksum_r1=$(cat "${r1_files_array[@]}" | tee "$tmp_r1" | md5sum | cut -d' ' -f1)

      # Now verify the file on disk matches the stream it came from.
      file_checksum_r1=$(md5sum "$tmp_r1" | cut -d' ' -f1)

      if [[ "$stream_checksum_r1" != "$file_checksum_r1" ]]; then
        echo "  ERROR: R1 checksum mismatch during consolidation." >&2
        rm -f "$tmp_r1"
        ((samples_skipped++))
        continue
      else
        # Verification passed, atomically move the file.
        mv "$tmp_r1" "$output_r1"
        echo "  [OK] R1 consolidation and verification passed."
      fi

      # Consolidate and verify R2 files
      echo "  Consolidating and verifying R2 files..."
      tmp_r2="${output_r2}.tmp"

      stream_checksum_r2=$(cat "${r2_files_array[@]}" | tee "$tmp_r2" | md5sum | cut -d' ' -f1)
      file_checksum_r2=$(md5sum "$tmp_r2" | cut -d' ' -f1)

      if [[ "$stream_checksum_r2" != "$file_checksum_r2" ]]; then
        echo "  ERROR: R2 checksum mismatch during consolidation." >&2
        rm -f "$tmp_r2" "$output_r1" # Clean up successfully created R1 file
        ((samples_skipped++))
        continue
      else
        mv "$tmp_r2" "$output_r2"
        echo "  [OK] R2 consolidation and verification passed."
      fi

      # Both R1 and R2 were successful, now remove originals
      echo "  Removing original files..."
      rm -f "${r1_files_array[@]}" "${r2_files_array[@]}"
      echo "  [OK] Original files removed"

      echo "  [OK] Sample consolidated successfully"
      ((samples_processed++))

    else
      echo "  [DRY-RUN] Would consolidate ${#r1_files_array[@]} R1 files and ${#r2_files_array[@]} R2 files"
      ((samples_processed++))

    fi

  # --- Process a single end read file ---
  else
    # Extract and sort single-end files

    se_files_list=$(printf "%s" "${sample_se_files[$sample_id]}" | sort -t: -k1,1n | cut -d: -f2-)

    # Convert to array for counting
    mapfile -t se_files_array <<< "$se_files_list"

    echo "  Read type: Single-end"
    echo "  Files found: ${#se_files_array[@]}"

    # Validate file count matches expected lanes
    if [[ ${#se_files_array[@]} -ne $LANES_PER_SAMPLE ]]; then
      echo "  ERROR: File count (${#se_files_array[@]}) does not match expected lanes ($LANES_PER_SAMPLE)"
      ((samples_skipped++))
      continue

    fi

    # Define output filename
    output_se="${FASTQ_DIRECTORY}${OUTPUT_PREFIX}_${sample_id}_NA${OUTPUT_SUFFIX}"

    # Check if output already exists
    if [[ -f "$output_se" ]]; then
      echo "  Output file already exists - skipping"
      echo "    $output_se"
      ((samples_skipped++))
      continue

    fi

    # Display files in lane order
    echo "  Files (in lane order):"
    printf '    %s\n' "${se_files_array[@]}"

    echo "  Will create:"
    echo "    $output_se"


    # Only consolidate if not in dry-run mode
    if [[ "$DRY_RUN" == false ]]; then
      # Consolidate and verify single-end files in one step
      echo "  Consolidating and verifying files..."
      tmp_se="${output_se}.tmp"

      # Use tee to write to the file and pipe to md5sum simultaneously
      stream_checksum_se=$(cat "${se_files_array[@]}" | tee "$tmp_se" | md5sum | cut -d' ' -f1)

      # Verify the file on disk matches the stream it came from.
      file_checksum_se=$(md5sum "$tmp_se" | cut -d' ' -f1)

      if [[ "$stream_checksum_se" != "$file_checksum_se" ]]; then
        echo "  ERROR: Checksum mismatch during consolidation." >&2
        rm -f "$tmp_se"
        ((samples_skipped++))
        continue
      else
        # Verification passed, atomically move the file.
        mv "$tmp_se" "$output_se"
        echo "  [OK] Consolidation and verification passed."
      fi

      # Verification is complete, now remove originals
      echo "  Removing original files..."
      rm -f "${se_files_array[@]}"
      echo "  [OK] Original files removed"
      echo "  [OK] Sample consolidated successfully"
      ((samples_processed++))

    else
      echo "  [DRY-RUN] Would consolidate ${#se_files_array[@]} files"
      ((samples_processed++))

    fi
  fi

  echo "--------------------"
  echo ""
done

echo "============================================"
if [[ "$DRY_RUN" == true ]]; then
  echo "Consolidation dry-run complete (no files created)"
else
  echo "Consolidation complete"
fi
echo "Total samples: ${#unique_sample_ids[@]}"
echo "Samples processed: $samples_processed"
echo "Samples skipped: $samples_skipped"
echo "============================================"

# ============================================
# Generate manifest file for successfully processed samples
# ============================================
if [[ "$DRY_RUN" == false && $samples_processed -gt 0 ]]; then
  echo "Generating manifest file..."

  MANIFEST_PATH="$DOCUMENTATION_DIR/$MANIFEST_FILENAME"

  # Create manifest with header
  if [[ "$IS_PAIRED_END" == true ]]; then
    echo -e "sample_id\tread1_path\tread2_path" > "$MANIFEST_PATH"
  else
    echo -e "sample_id\tread_path" > "$MANIFEST_PATH"
  fi

  # Add entries for each successfully consolidated sample
  manifest_entries=0
  for sample_id in "${unique_sample_ids[@]}"; do
    if [[ "$IS_PAIRED_END" == true ]]; then
      output_r1="${FASTQ_DIRECTORY}${OUTPUT_PREFIX}_${sample_id}_R1${OUTPUT_SUFFIX}"
      output_r2="${FASTQ_DIRECTORY}${OUTPUT_PREFIX}_${sample_id}_R2${OUTPUT_SUFFIX}"

      # Only add to manifest if both files exist
      if [[ -f "$output_r1" && -f "$output_r2" ]]; then
        echo -e "${sample_id}\t${output_r1}\t${output_r2}" >> "$MANIFEST_PATH"
        ((manifest_entries++))
      fi
    else
      output_se="${FASTQ_DIRECTORY}${OUTPUT_PREFIX}_${sample_id}_NA${OUTPUT_SUFFIX}"

      # Only add to manifest if file exists
      if [[ -f "$output_se" ]]; then
        echo -e "${sample_id}\t${output_se}" >> "$MANIFEST_PATH"
        ((manifest_entries++))
      fi
    fi
  done

  echo "[OK] Manifest created: $MANIFEST_PATH"
  echo "  Entries in manifest: $manifest_entries"

  # Verify manifest entry count matches processed samples
  if [[ $manifest_entries -ne $samples_processed ]]; then
    echo "  WARNING: Manifest entries ($manifest_entries) does not match samples processed ($samples_processed)"
  fi

elif [[ "$DRY_RUN" == true ]]; then
  echo "Manifest generation skipped (dry-run mode)"
  echo "Would create: $DOCUMENTATION_DIR/$MANIFEST_FILENAME"
else
  echo "No samples processed - manifest not created"
fi

echo ""
echo "============================================"
echo "Script complete"
echo "============================================"
