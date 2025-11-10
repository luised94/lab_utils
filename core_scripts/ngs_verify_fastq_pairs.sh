#!/bin/bash
################################################################################
# Verify paired end fastq files and output manifest
# Author: Luis | Date: 2025-10-20 | Version: 1.0.0
################################################################################
# PURPOSE:
#   For a fastq directory, find paired end files according to BMC naming convention and verify that paired end files have same number of reads (simplest check)
#   Output manifest of paired files for subsequent scripts.
# USAGE:
#   From the command line
#   $ ./ngs_verify_fastq_pairs.sh /path/to/fastq/directory
# DEPENDENCIES:
#   bash 4.2
# OUTPUTS:
#   No outputs produced by file or script.
################################################################################
# ============================================
# Usage and Help
# ============================================
show_usage() {
  cat << EOF
Usage: srun $(basename "$0") <fastq_directory> [-v]

Description:
  Verifies FASTQ file structure, detects single/paired-end reads,
  validates read counts, and generates manifest for downstream processing.

Arguments:
  EXPERIMENT_ID    Experiment ID for BMC sequencing project
                     (e.g., 250930Bel)

Options:
  -h, --help        Show this help message
  -v                Verbose output (show detailed processing steps)

Output:
  Creates paired_reads_manifest.tsv in <experiment>/documentation/
  Skips generation if manifest already exists (delete to regenerate)

Example:
  $(basename "$0") ~/data/250930Bel/fastq
  $(basename "$0") ~/data/250930Bel/fastq -v

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
elif [[ $# -gt $MAX_NUMBER_OF_ARGS ]]; then
    echo "Error: Too many arguments provided ($#)." >&2
    show_usage
    exit 1
fi

# Handle first argument: Remove trailing slash and validate pattern
EXPERIMENT_ID=${1%/}
echo "Running error handling..."
if [[ ! $EXPERIMENT_ID =~ $EXPECTED_EXPERIMENT_ID_PATTERN ]]; then
  echo "Error: EXPERIMENT_ID does not match expected pattern." >&2
  echo "Please adjust EXPERIMENT_ID accordingly." >&2
  echo "EXPERIMENT ID PATTERN: $EXPECTED_EXPERIMENT_ID_PATTERN" >&2
  echo "EXPERIMENT_ID: $EXPERIMENT_ID" >&2
  exit 1
fi

VERBOSE=false
if [[ $# -eq $MAX_NUMBER_OF_ARGS ]]; then
  if [[ "$2" != "-v" ]]; then
    echo "Error: Invalid second argument: '$2'" >&2
    echo "Use -v to perform output verbose output." >&2
    echo "Run '$0 -h' for usage." >&2
    exit 1
  fi
  VERBOSE=true
fi

# Ensure script is run inside a Slurm allocation
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    cat >&2 <<EOF
Error: This script must be run within a Slurm job.

To run interactively:
    srun $0 <EXPERIMENT_ID> [-v]

To submit as a batch job:
    echo "$0 <EXPERIMENT_ID>" [-v] | sbatch

EOF
  exit 1
fi

#============================== 
# Configuration
#============================== 
echo "-------Start $0-------"
echo "Setting up configuration..."
NCORES=$(nproc)

#============================== 
# Setup and preprocessing
#============================== 
# Get directory from first argument
EXPERIMENT_ID=${1%/}
EXPERIMENT_DIR="$HOME/data/${EXPERIMENT_ID}"
FASTQ_DIRECTORY="$EXPERIMENT_DIR/fastq/"
DOCUMENTATION_DIR="$(dirname "$FASTQ_DIRECTORY")/documentation"
MANIFEST_FILE="$DOCUMENTATION_DIR/paired_reads_manifest.tsv"

# Ensure DOCUMENTATION_DIR exists
mkdir -p "$DOCUMENTATION_DIR"



# ============================================
# Validation and Error Checking
# ============================================
echo "-------Start $(basename "$0")-------"
echo "Validating inputs..."

# Check directory creation succeeded
if [[ ! -d "$DOCUMENTATION_DIR" ]]; then
  echo "ERROR: Failed to create documentation directory: $DOCUMENTATION_DIR"
  exit 1

fi

# Check directory exists
if [[ ! -e "$FASTQ_DIRECTORY" ]]; then
  echo "ERROR: Path does not exist: $FASTQ_DIRECTORY"
  exit 1
fi

# Check it's actually a directory
if [[ ! -d "$FASTQ_DIRECTORY" ]]; then
  echo "ERROR: Path is not a directory: $FASTQ_DIRECTORY"
  exit 1
fi

# Check directory is readable
if [[ ! -r "$FASTQ_DIRECTORY" ]]; then
  echo "ERROR: Directory is not readable: $FASTQ_DIRECTORY"
  exit 1
fi

# ############################################
# Main logic
# ############################################
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
  #
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
      lane_map["$lane_number"]=1

  else
      # If there's no lane, you can assign a default value
      lane_map["NA"]=1

  fi


  # --- Extract read pair indicator ---
  read_indicator="${parts[$READ_PAIR_IDX]}"
  pair_indicator_map["$read_indicator"]=1

done # end read metadata for loop

echo "Extracting unique metadata values..."
mapfile -t unique_sample_ids < <(printf '%s\n' "${!sample_id_map[@]}" | sort)
mapfile -t detected_lanes < <(printf '%s\n' "${!lane_map[@]}" | sort -n)
mapfile -t unique_pair_indicator < <(printf '%s\n' "${!pair_indicator_map[@]}" | sort)
EXPECTED_LANES_PER_SAMPLE=${#detected_lanes[@]}

# ============================================
# Error handling for metadata extraction
# ============================================
if [[ ${#unique_sample_ids[@]} -eq 0 ]]; then
  echo "ERROR: No valid sample IDs extracted from filenames"
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
  echo "ERROR: No lanes detected after extraction."
  exit 1

fi

if [[ "$VERBOSE" == true ]]; then
  echo "Unique samples found: ${#unique_sample_ids[@]}"
  echo "Sample IDs:"
  printf '  %s\n' "${unique_sample_ids[@]}"
  echo "  Lanes found: ${detected_lanes[*]}"
  echo "  Found read indicators: ${unique_pair_indicator[*]}"
  echo "  Expected lanes per sample: $EXPECTED_LANES_PER_SAMPLE"

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
  echo "  Detected: PAIRED-END"
elif [[ "$has_single_indicator" == true ]]; then
  IS_PAIRED_END=false
  EXPECTED_READ_PAIRS_PER_LANE=1
  echo "  Detected: SINGLE-END"
else
  echo "  ERROR: No valid read indicators found"
  exit 1
fi

# Check read type detection succeeded
if [[ -z "$IS_PAIRED_END" ]]; then
  echo "ERROR: Failed to detect read type"
  exit 1
fi

echo "Verifying file counts per sample..."
EXPECTED_FILES_PER_SAMPLE=$((EXPECTED_LANES_PER_SAMPLE * EXPECTED_READ_PAIRS_PER_LANE))
echo "  Expected files per sample: $EXPECTED_FILES_PER_SAMPLE"

# ============================================
# Verify how many files there are per sample
# ============================================
declare -A sample_file_lists  # Will store space-separated file lists

for sample_id in "${unique_sample_ids[@]}"; do
  # Find all files for this sample
  mapfile -t files_for_sample < <(
    printf '%s\n' "${all_fastq_files[@]}" | grep "$sample_id" | grep -v "unmapped"
  )

  file_count=${#files_for_sample[@]}

  if [[ $file_count -ne $EXPECTED_FILES_PER_SAMPLE ]]; then
    echo "  [WARNING] $sample_id: found $file_count files, expected $EXPECTED_FILES_PER_SAMPLE"

  else

    if [[ "$VERBOSE" == true ]]; then
      echo "  [OK] $sample_id: $file_count files"

    fi

  fi

  # Store for later use (joining array into string)
  sample_file_lists["$sample_id"]="${files_for_sample[*]}"
done

# ============================================
# Process based on read type
# ============================================
# ============================================
# Process based on read type
# ============================================
if [[ "$IS_PAIRED_END" == true ]]; then
  echo ""
  echo "Processing paired-end reads..."

  # ============================================
  # Generate Paired Reads Manifest and Prepare for Verification
  # ============================================
  echo "Pairing files and generating manifest: $MANIFEST_FILE"

  if [[ -f "$MANIFEST_FILE" ]]; then
    echo "Manifest already exists: $MANIFEST_FILE"
    echo "  To regenerate, delete the file and rerun this script"
    echo "  rm \"$MANIFEST_FILE\""
  else
    # Write header for the manifest file
    echo -e "sample_id\tlane\tread1_path\tread2_path" > "$MANIFEST_FILE"
    
    # Create a temporary file to hold pairs for the read count verification step
    temp_pairs_file=$(mktemp)

    for sample_id in "${unique_sample_ids[@]}"; do
      echo "  Processing sample: $sample_id"
      IFS=' ' read -ra sample_files <<< "${sample_file_lists[$sample_id]}"

      for lane_num in "${detected_lanes[@]}"; do
        r1="" r2=""
        
        # Find the R1/R2 pair based on the detected naming format
        if [[ "$LANE_FORMAT" == "numeric" ]]; then
          # Format includes lane numbers, so we match on the lane
          for file in "${sample_files[@]}"; do
            IFS='_-' read -ra parts <<< "$(basename "$file")"
            if [[ "${parts[$LANE_IDX]}" == "$lane_num" ]]; then
              [[ "${parts[$READ_PAIR_IDX]}" == "1" ]] && r1="$file"
              [[ "${parts[$READ_PAIR_IDX]}" == "2" ]] && r2="$file"
            fi
          done
        else # LANE_FORMAT is "none"
          # Format does not include lanes, so we only look for R1/R2
          for file in "${sample_files[@]}"; do
            IFS='_-' read -ra parts <<< "$(basename "$file")"
            [[ "${parts[$READ_PAIR_IDX]}" == "1" ]] && r1="$file"
            [[ "${parts[$READ_PAIR_IDX]}" == "2" ]] && r2="$file"
          done
        fi

        # If a valid pair is found, write to manifest and verification file
        if [[ -n "$r1" && -n "$r2" ]]; then
          if [[ "$VERBOSE" == true ]]; then
            echo "    Found pair for lane $lane_num: $(basename $r1) <-> $(basename $r2)"
          fi
          # Write the correct data to the manifest
          echo -e "$sample_id\t$lane_num\t$r1\t$r2" >> "$MANIFEST_FILE"
          # Write data for the verification step
          echo "$sample_id|$r1|$r2" >> "$temp_pairs_file"
        elif [[ -n "$r1" || -n "$r2" ]]; then
          echo "    [WARNING] Unmatched file for sample $sample_id in lane $lane_num."
        fi

      done
    done
    
    echo "Manifest complete: $(($(wc -l < "$MANIFEST_FILE") - 1)) pairs written"

    # ============================================
    # Verify read counts match
    # ============================================
    echo ""
    echo "Verifying read counts in pairs..."
    
    # Use xargs to run in parallel, messages are not sorted.
    cat "$temp_pairs_file" | xargs -P "$NCORES" -I {} bash -c '
      IFS="|" read -r sample r1 r2 <<< "{}"
      r1_lines=$(wc -l < "$r1")
      r2_lines=$(wc -l < "$r2")
      r1_reads=$((r1_lines / 4))
      r2_reads=$((r2_lines / 4))
      if [[ $r1_reads -eq $r2_reads ]]; then
        echo "OK|$sample|$(basename "$r1")|$(basename "$r2")|$r1_reads"
      else
        echo "ERROR|$sample|$r1_reads|$r2_reads"
      fi
    ' | while IFS="|" read -r status sample_or_r1 file1_or_r2 file2_or_count rest; do
        if [[ "$status" == "OK" ]]; then
          [[ "$VERBOSE" == true ]] && echo "  [OK] $sample_or_r1: $file1_or_r2 <-> $file2_or_count: $rest reads"
        else
          echo "  [ERROR] Read count mismatch for $sample_or_r1: R1=$file1_or_r2, R2=$file2_or_count"
        fi
      done
    
    # Clean up the temporary file
    rm "$temp_pairs_file"
  fi

else
  echo ""
  echo "Single-end mode - skipping pairing and read count verification"
  echo "Note: Single-end processing not yet implemented"

fi

# --- Final summary ---
echo ""
echo "============================================"
echo "Summary"
echo "============================================"
echo "  Samples processed: ${#unique_sample_ids[@]}"
echo "  Read type: $([ "$IS_PAIRED_END" == true ] && echo "PAIRED-END" || echo "SINGLE-END")"
echo "  Lanes per sample: ${#detected_lanes[@]}"
echo "  Manifest: $MANIFEST_FILE"
echo "  Ncores: $NCORES"
echo "For complete output, pass -v option as second argument position."
echo "-------End $(basename "$0")-------"
