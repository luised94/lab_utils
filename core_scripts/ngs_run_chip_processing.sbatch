#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --nice=10000
#SBATCH --exclude=c[5-22]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luised94@mit.edu
################################################################################
# ngs_run_chip_processing.sbatch
# Run CHIP fastq processing pipeline
# Author: Luis | Date: 2025-10-23 | Version: 2.0.0
################################################################################
# PURPOSE:
#   Filter, align and calculate bigwig coverage for fastq files from a given experiment id.
# USAGE:
#   From the command line, run via the slurm submission script. Although can submit a particular slurm task directly
#   $ ngs_submit_chip_processing.sh <EXPERIMENT_ID>
#   $ srun ngs_run_chip_processing.sh <EXPERIMENT_ID>
#   $ sbatch --array=1-N%16 ngs_run_chip_processing.sh <EXPERIMENT_ID>
# DEPENDENCIES:
#   bash 4.2, bowtie2, fastp, bedtools, samtools
#   Assumes fastq files are in the fastq directory and everything is removed. (cleanup script.)
#   Assumes manifest file has been generated.
# OUTPUTS:
#   Sorted and indexed bam files, bigwig files for visualization.
#   fastp files generates quality control files for filtering process. To sync html and json files locally, run 'rsync -nav username@luria.mit.edu:~/logs/* ~/logs/'.
################################################################################
show_usage() {
  cat << EOF
Usage: 
2) ngs_submit_chip_processing.sh <EXPERIMENT_ID>
1) srun $0 <EXPERIMENT_ID>
2) sbatch $0 <EXPERIMENT_ID>

Process consolidated FASTQ files via the CHIP processing pipeline.

Arguments:
  EXPERIMENT_ID    Full path to experiment directory
                         Example: \$HOME/data/250930Bel
Options:
  -h, --help        Show this help message

Requirements:
  - Must run on luria cluster as sbatch job
  - Must be in git repository (lab_utils)
  - Consolidated FASTQ files must exist
  - Manifest file must exist: documentation/consolidated_reads_manifest.tsv

The script will:
  1. Read manifest to determine samples and read type
  2. Ensure file specified by slurm task id exists.
  3. Determines command options.
  4. Process fastq files with fastp, bowtie2, bamcoverage and bedtools
EOF
  exit 0

}

# Validate SLURM_ARRAY_TASK_ID
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  echo "Error: This script must be run as a SLURM array job"
  echo "Use: sbatch --array=1-N%16 $0 <experiment_directory>"
  exit 1

fi


#==============================
# Argument error handling
#==============================
# Check for arguments
echo "Handling arguments..."
MIN_NUMBER_OF_ARGS=1
MAX_NUMBER_OF_ARGS=1
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
#DRY_RUN=true
#if [[ $# -eq $MAX_NUMBER_OF_ARGS ]]; then
#    if [[ "$2" != "--active-run" ]]; then
#        echo "Error: Unknown option '$2'" >&2
#        echo "Use --active-run to perform actual sync." >&2
#        exit 1
#    fi
#    DRY_RUN=false
#fi

# Capture initial environment variables for logging
# Any variables at this time will not be logged.
# Everything after this mapfile command will be!
mapfile -t initial_vars_declarations < <(compgen -A variable | sort)

#==============================
# Configuration
#==============================
echo "Setting configuration..."
CPU_THREADS="$SLURM_CPUS_PER_TASK"
THREADS=$(nproc)
# SLURM job configuration
#REPO_DIRECTORY_PATH=$(git rev-parse --show-toplevel)
#PROCESSING_SCRIPT_TO_SUBMIT="${REPO_DIRECTORY_PATH}/core_scripts/ngs_run_chip_processing.sbatch"
#NGS_FILE_PATTERNS=("consolidated*.fastq" "*.bam" "*.bw")

# --- Command configuration ---
# Any settings not explictly defined can be assumed to be default.
# fastp parameters
MAXIMUM_N_BASE_COUNT=5
MAXIMUM_UNQUALIFIED_BASE_PERCENT=50
COMPLEXITY_WINDOW_SIZE=4
OVERREPRESENTATION_SAMPLING=50
MINIMUM_BASE_QUALITY=20
MINIMUM_READ_LENGTH=25

# bowtie2 parameters
ALIGNMENTS_TO_REPORT=1
EXTENSIONS_TO_TRY=15
SETS_OF_SEEDS=2
MAX_PENALTY_MISMATCH=4
NON_NUCLEOTIDE_PENALTY=1

# bamCoverage parameters
BIN_SIZE=10
EFFECTIVE_GENOME_SIZE=12157105
MINIMUM_MAPPING_QUALITY=20
NORM_METHOD="CPM"

# Genome reference files (Saccharomyces cerevisiae S288C)
GENOME_DIR="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C"
GENOME_INDEX="$GENOME_DIR/SaccharomycescerevisiaeS288C_index"
BLACKLIST_BED_FILE="$HOME/data/feature_files/20250423_merged_saccharomyces_cerevisiae_s288c_blacklist.bed"
FASTQ_FILEPATH_PATTERN="consolidated*.fastq"
#SAMPLE_ID_BASH_REGEX='^consolidated_(.*)_sequence\.fastq$'
#SAMPLE_ID_SED_REGEX= "consolidated_\(.*\)_sequence\.fastq"

# Manifest output configuration
MANIFEST_FILENAME="consolidated_reads_manifest.tsv"

# Directory structure
EXPERIMENT_DIR="$HOME/data/${EXPERIMENT_ID}"
FASTQ_DIR="${EXPERIMENT_DIR}/fastq"
ALIGNMENT_DIR="${EXPERIMENT_DIR}/alignment"
COVERAGE_DIR="${EXPERIMENT_DIR}/coverage"
DOCUMENTATION_DIR="${EXPERIMENT_DIR}/documentation"
MANIFEST_FILEPATH="${DOCUMENTATION_DIR}/$MANIFEST_FILENAME"
QUALITY_CONTROL_DIR="${EXPERIMENT_DIR}/quality_control"
THREADS=$( nproc )

#=======================
# Setup logging
#=======================
# Read in logging functions and then setup
TOOL_NAME="chip_processing"
source "$HOME/lab_utils/core_scripts/functions_for_logging.sh"
# Sets log file in the logs directory. Outputs LOG_DIR.
eval "$(setup_logging "${EXPERIMENT_ID}" ${TOOL_NAME} ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID} )"

# @QUES: I dont think this is necessary.
#log_message "INFO" "=== Starting ChIP-seq processing pipeline ==="
#log_message "INFO" "Experiment: $EXPERIMENT_ID"
#log_message "INFO" "Array task ID: $SLURM_ARRAY_TASK_ID"

# --- Check REQUIRED FILES ---
# Ensure genome index exists.
if [[ ! -f "${GENOME_INDEX}.1.bt2" ]]; then
  log_message "ERROR" "Genome index not found: $GENOME_INDEX"
  exit 1

fi

# Check if blacklist file exists
if [[ ! -f "$BLACKLIST_BED_FILE" ]]; then
  log_message "ERROR" "Blacklist file not found: $BLACKLIST_BED_FILE. Blacklist runs will be skipped."
  exit 1

fi

log_message "INFO" "Validation complete - all required files present"

# --- Ensure directories exist ---
mkdir -p "$FASTQ_DIR"
mkdir -p "$ALIGNMENT_DIR"
mkdir -p "$QUALITY_CONTROL_DIR"
mkdir -p "$COVERAGE_DIR"
#============================== 
# READ MANIFEST AND EXTRACT SAMPLE INFORMATION
#============================== 
# Validate manifest exists
if [[ ! -f "$MANIFEST_FILEPATH" ]]; then
  log_message "ERROR" "Manifest file not found: $MANIFEST_FILEPATH"
  exit 1

fi

# Get manifest row for this array task (add 1 to skip header)
ROW_NUMBER=$((SLURM_ARRAY_TASK_ID + 1))
MANIFEST_ROW=$(sed -n "${ROW_NUMBER}p" "$MANIFEST_FILEPATH")
if [[ -z "$MANIFEST_ROW" ]]; then
  log_message "ERROR" "No data found for task ID ${SLURM_ARRAY_TASK_ID} in manifest"
  log_message "ERROR" "Manifest file: $MANIFEST_FILEPATH"
  exit 1
fi

# Parse manifest columns
IFS=$'\t' read -r SAMPLE_ID read1_path read2_path <<< "$MANIFEST_ROW"

# Validate minimum fields
if [[ -z "$SAMPLE_ID" ]] || [[ -z "$read1_path" ]]; then
  log_message "ERROR" "Manifest row must have at least 2 fields: $MANIFEST_ROW"
  exit 1
fi

log_message "INFO" "Sample ID: $SAMPLE_ID"
log_message "INFO" "Read 1: $read1_path"
if [[ -n "$read2_path" ]]; then
  log_message "INFO" "Read 2: $read2_path"
fi

#========================
# Detect read type and build tool parameters
#========================

if [[ -z "$read2_path" ]]; then
  # Single-end configuration
  read_type="single-end"
  log_message "INFO" "Processing as single-end reads"

  # Validate input exists
  if [[ ! -f "$read1_path" ]]; then
      log_message "ERROR" "Read1 file not found: $read1_path"
      exit 1
  fi

  # Define output paths
  fastp_out1="${fastq_dir}/fastpfiltered_${SAMPLE_ID}_R1.fastq"

  # Build tool parameters
  fastp_input_params="--in1 ${read1_path}"
  fastp_output_params="--out1 ${fastp_out1}"
  bowtie2_input_params="-U ${fastp_out1}"

else
  # Paired-end configuration
  read_type="paired-end"
  log_message "INFO" "Processing as paired-end reads"

  # Validate inputs exist
  if [[ ! -f "$read1_path" ]]; then
      log_message "ERROR" "Read1 file not found: $read1_path"
      exit 1
  fi

  if [[ ! -f "$read2_path" ]]; then
      log_message "ERROR" "Read2 file not found: $read2_path"
      exit 1
  fi

  # Define output paths
  fastp_out1="${fastq_dir}/fastpfiltered_${SAMPLE_ID}_R1.fastq"
  fastp_out2="${fastq_dir}/fastpfiltered_${SAMPLE_ID}_R2.fastq"

  # Build tool parameters
  fastp_input_params="--in1 ${read1_path} --in2 ${read2_path}"
  fastp_output_params="--out1 ${fastp_out1} --out2 ${fastp_out2}"
  bowtie2_input_params="-1 ${fastp_out1} -2 ${fastp_out2}"

fi

#===================================================
# DEFINE OUTPUT FILE PATHS (COMMON FOR BOTH READ TYPES)
#==========================================
# Alignment outputs
SORTED_BAM_FILE="${ALIGNMENT_DIR}/${SAMPLE_ID}_to_S288C_sorted.bam"
BLACKLIST_FILTERED_BAM_PATH="${ALIGNMENT_DIR}/${SAMPLE_ID}_blFiltered.bam"
# Coverage output
BIGWIG_PATH="${COVERAGE_DIR}/${SAMPLE_ID}_${NORM_METHOD}.bw"

# QC outputs
FASTP_JSON="${QUALITY_CONTROL_DIR}/${SAMPLE_ID}_fastp.json"
FASTP_HTML="${QUALITY_CONTROL_DIR}/${SAMPLE_ID}_fastp.html"

# I dont think these are necessary
#log_message "INFO" "Output paths configured"
#log_message "INFO" "  Sorted BAM: $(basename "$SORTED_BAM_FILE")"
#log_message "INFO" "  Filtered BAM: $(basename "$BLACKLIST_FILTERED_BAM_PATH")"
#log_message "INFO" "  BigWig: $(basename "$BIGWIG_PATH")"

#===================================================
# Log configuration variables
#==========================================
log_message "CONFIG" "=== Script Configuration ==="
mapfile -t variables_after_config < <(compgen -A variable | sort)

while IFS= read -r variable_name; do
  # Skip empty lines
  [[ -z "$variable_name" ]] && continue

  # Skip internal loop variables and arrays
  case "$variable_name" in
      initial_vars_declarations|variables_after_config|variable_name|variable_value|formatted_log_line)
          continue
          ;;
  esac

  # Get variable value
  variable_value="${!variable_name}"

  # Format and log
  printf -v formatted_log_line "%s=%q" "$variable_name" "$variable_value"
  log_message "CONFIG" "$formatted_log_line"

done < <(comm -13 <(printf "%s\n" "${initial_vars_declarations[@]}") \
                 <(printf "%s\n" "${variables_after_config[@]}"))

log_message "CONFIG" "=== End Configuration ==="


BAMCOVERAGE_COMMON_PARAMS="--bam ${BLFILTERED_BAM_PATH} \
  --outFileName ${BAMCOVERAGE_PATH} \
  --binSize ${BIN_SIZE} \
  --ignoreDuplicates \
  --minMappingQuality ${MINIMUM_MAPPING_QUALITY} \
  --normalizeUsing ${NORM_METHOD} \
  --numberOfProcessors ${SLURM_CPUS_PER_TASK}"

# Add RPGC-specific parameters
if [ "${NORM_METHOD}" == "RPGC" ]; then
  BAMCOVERAGE_COMMON_PARAMS+=" --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE}"
  log_message "INFO" "Added effectiveGenomeSize parameter for RPGC normalization"
fi


# Get current fastq file
mapfile -t PREFILTERED_FASTQ_FILEPATHS < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "$FASTQ_FILEPATH_PATTERN" | sort)
# Get half of the number of files, then use that to get the rest
FASTQ_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
PREFILTERED_FASTQ_PATH="${PREFILTERED_FASTQ_FILEPATHS[$FASTQ_INDEX]}"
FASTQ_BASENAME=$(basename "$PREFILTERED_FASTQ_PATH")

# Grab sample id from fastq filepath
#SAMPLE_ID=$(echo "$FASTQ_BASENAME" | sed -n "s/$SAMPLE_ID_SED_REGEX/\\1/p")
if [[ "$FASTQ_BASENAME" =~ $SAMPLE_ID_BASH_REGEX ]]; then
  SAMPLE_ID="${BASH_REMATCH[1]}"
else 
  log_message "ERROR" "Invalid input filename format: ${FASTQ_BASENAME}"
  log_message "ERROR" "Expected format: consolidate_<ID>_sequence.fastq"
  log_message "ERROR" "SAMPLE ID not extracted. Check fastq basename or SAMPLE_ID_BASH_REGEX"
  exit 1

fi

# Set the output path files
# Generate output name

log_message "CONFIG" "=== Script Configuration ==="
mapfile -t variables_after_config < <(compgen -A variable | sort)

# We pipe this into a while loop to process each new variable
while IFS= read -r variable_name; do
  # Skip empty lines
  [[ -z "$variable_name" ]] && continue
  # We don't want to log these as they're not part of your pipeline config
  case "$variable_name" in
  initial_vars_declarations|variables_after_config|variable_name|variable_value|formatted_log_line)
    continue
    ;;
  esac

  # If we set a prefix requirement, check if this variable matches
  # If VAR_PREFIX is "PIPELINE_", this only logs variables starting with PIPELINE_
  #if [[ -n "$VAR_PREFIX" ]]; then
  #  [[ "$variable_name" =~ ^${VAR_PREFIX} ]] || continue
  #fi

  # Get the actual value of this variable using indirect expansion
  variable_value="${!variable_name}"

  # Format the variable for logging: use %q to safely quote
  printf -v formatted_log_line "%s=%q" "$variable_name" "$variable_value"

  # Write to both console and log file
  log_message "CONFIG" "$formatted_log_line"

  # The comm command compares two sorted files line by line
  # -1 suppresses lines unique to first file (initial_vars_declarations)
  # -3 suppresses lines common to both files
done < <(comm -13 <(printf "%s\n" "${initial_vars_declarations[@]}") \
                 <(printf "%s\n" "${variables_after_config[@]}"))
log_message "CONFIG" "=== End Configuration ==="

# Setup additional parameters based on samples
# Error handling
TOTAL_FILES=${#PREFILTERED_FASTQ_FILEPATHS[@]}
if [ "$TOTAL_FILES" -eq 0 ]; then
  log_message "ERROR" "No fastq files found in ${FASTQ_DIR}"
  log_message "ERROR" "Check the following:"
  log_message "ERROR" " - Files were consolidated using script 'bmc_consolidate_fastq_by_id.sh'"
  log_message "ERROR" " - FInd name expression: $FASTQ_FILEPATH_PATTERN"
  log_message "ERROR" "---------------------"
  exit 1

fi

if [[ ! -f "$FASTQ_PATH" ]]; then
  log_message "ERROR" "No fastq file found for index $FASTQ_INDEX"
  exit 1

fi

if [[ ! "$SAMPLE_ID" =~ ^[0-9]{1,6}$ ]]; then
  log_message "ERROR" "Invalid or missing sample ID from filename: ${FASTQ_BASENAME}"
  log_message "ERROR" "SAMPLE_ID: ${SAMPLE_ID}"
  exit 1

fi

log_message "INFO" "If reached this step, verify the configuration variables before continuing..."
exit 0

# Load required modules
module purge
module load fastp/0.20.0

# Fastp filtering
log_message "INFO" "Starting fastp filtering"
# Options are not available in fastp 0.20.0 version available in the linux cluster.
#--dedup \
#--dup_calc_accuracy "$DUPLICATION_CALC_ACCURACY" \
#--compression "$COMPRESSION_LEVEL" \
# Add --in2 --out2 when I figure out how the paired end reads are delivered.
if measure_performance "fastp_filtering" \
  fastp \
      --in1 "$FASTQ_PATH" \
      --out1 "$FASTP_FILTERED_FASTQ_PATH" \
      --adapter_sequence auto \
      --cut_window_size "$COMPLEXITY_WINDOW_SIZE" \
      --cut_mean_quality "$MINIMUM_BASE_QUALITY" \
      --cut_front \
      --cut_tail \
      --cut_right \
      --n_base_limit "$MAXIMUM_N_BASE_COUNT" \
      --average_qual "$MINIMUM_BASE_QUALITY" \
      --qualified_quality_phred "$MINIMUM_BASE_QUALITY" \
      --unqualified_percent_limit "$MAXIMUM_UNQUALIFIED_BASE_PERCENT" \
      --length_required "$MINIMUM_READ_LENGTH" \
      --thread "$CPU_THREADS" \
      --overrepresentation_analysis \
      --overrepresentation_sampling "$OVERREPRESENTATION_SAMPLING" \
      --json "$QUALITY_CONTROL_DIR/${SAMPLE_ID}_fastp.json" \
      --html "$QUALITY_CONTROL_DIR/${SAMPLE_ID}_fastp.html"; then
  log_message "INFO" "Performed fastp filtering."

else
  log_message "ERROR" "Fastp filtering failed for ${SAMPLE_ID}"
  exit 1

fi
log_message "INFO" "Fastp complete..."

# Load required modules
module purge
module load bowtie2
module load samtools

#if [ -f "${OUTPUT_BAM}" ]; then
#  log_message "SKIP" "File already exists: ${OUTPUT_BAM}"
#  log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
#  exit 0
#fi

# Alignment and sorting
log_message "INFO" "Starting alignment and sorting"
# Error handling
if measure_performance "alignment_and_sorting" \
  bowtie2 -x "$GENOME_INDEX" \
          -U "$FASTP_FILTERED_FASTQ_PATH" \
          --mp $MAX_PENALTY_MISMATCH \
          --np $NON_NUCLEOTIDE_PENALTY \
          --sensitive \
          -D $EXTENSIONS_TO_TRY \
          -R $SETS_OF_SEEDS \
          -k $ALIGNMENTS_TO_REPORT \
          -p "$CPU_THREADS" \
          2>> "${ERROR_LOG}" | \
  samtools view -@ "$SLURM_CPUS_PER_TASK" -bS - 2>> "${ERROR_LOG}" | \
  samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$PREBLACKLIST_FILTERED_BAM_PATH" - 2>> "${ERROR_LOG}"; then
  log_message "INFO" "Starting BAM indexing"
  if measure_performance "indexing" samtools index "$PREBLACKLIST_FILTERED_BAM_PATH"; then
      log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
  else
      log_message "ERROR" "BAM indexing failed for ${SAMPLE_NAME}"
      exit 1
  fi
else
  log_message "ERROR" "Alignment/sorting failed for ${SAMPLE_NAME}"
  exit 1
fi

log_message "INFO" "Bowtie2 complete..."

# Start the
# Load required modules
module purge
module load python/2.7.13
module load deeptools/3.0.1
module load samtools

# Execute blacklist filtering ------
# Check if BAM file already exists
#if [ -f "$OUTPUT_BAM" ]; then
#  log_message "SKIP" "File already exists: ${OUTPUT_BAM}"
#
#  if [ -f "${OUTPUT_BAM}.bai" ]; then
#    log_message "SKIP" "File already exists: ${OUTPUT_BAM}.bai"
#    exit 0
#  fi
#
#  if ! measure_performance "index" samtools index "$OUTPUT_BAM"; then
#    log_message "ERROR" "BAM indexing failed for ${SAMPLE_NAME}"
#    exit 1
#  fi
#
#  log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
#  exit 0
#fi

# Perform blacklist filtering
log_message "INFO" "Starting blacklist filtering processing"
if ! measure_performance "blacklist_filtering" \
  alignmentSieve --bam "${BAM_PATH}" \
               --blackListFileName "$BLACKLIST_BED_FILE" \
               --outFile "$OUTPUT_BAM" \
               --numberOfProcessors $(( THREADS / 2 )); then
  log_message "ERROR" "Blacklist processing failed for ${SAMPLE_NAME}"
  exit 1
fi
log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"

# Index the resulting BAM file
if ! measure_performance "index" samtools index "$OUTPUT_BAM"; then
  log_message "ERROR" "Blacklist BAM indexing failed for ${SAMPLE_NAME}"
  exit 1
fi
log_message "INFO" "Successfully indexed for ${SAMPLE_NAME}"

# Log completion
log_message "INFO" "Completed alignementSieve blacklist filtering..."

# Load required modules
module purge
module load python/2.7.13
module load deeptools/3.0.1

#if [ -f "$OUTPUT_BIGWIG" ]; then
#  log_message "SKIP" "File already exists: ${OUTPUT_BIGWIG}"
#  log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
#  exit 0
#fi

# Execute bamCoverage
log_message "INFO" "Starting bamCoverage processing"
# shellcheck disable=SC2086 # Intentional: $BAMCOVERAGE_COMMON_PARAMS should be word-split
if ! measure_performance "bamcoverage" bamCoverage $BAMCOVERAGE_COMMON_PARAMS; then
  log_message "bamCoverage processing failed for ${SAMPLE_NAME}"
  exit 1

fi

# Log completion
log_message "INFO" "bamCoverage command complete..."
log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
