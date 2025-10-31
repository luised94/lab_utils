#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --nice=10000
#SBATCH --exclude=c[5-22]
#SBATCH --output=/home/%u/data/logs/slurm_%A_%a.out
#SBATCH --error=/home/%u/data/logs/slurm_%A_%a.out
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

# Capture initial environment variables for logging
# Any variables at this time will not be logged.
# Everything after this mapfile command will be!
mapfile -t initial_vars_declarations < <(compgen -A variable | sort)

#==============================
# Configuration
#==============================
echo "Setting configuration..."
CPU_THREADS="$SLURM_CPUS_PER_TASK"

# --- Command configuration ---
# Any settings not explictly defined can be assumed to be default.
# fastp parameters
MAXIMUM_N_BASE_COUNT=5
MAXIMUM_UNQUALIFIED_BASE_PERCENT=50
COMPLEXITY_WINDOW_SIZE=4
OVERREPRESENTATION_SAMPLING=50
MINIMUM_BASE_QUALITY=20
MINIMUM_READ_LENGTH=25
FASTP_OUT_PREFIX="fastpfiltered"

# bowtie2 parameters
ALIGNMENTS_TO_REPORT=1
EXTENSIONS_TO_TRY=15
SETS_OF_SEEDS=2
MAX_PENALTY_MISMATCH=4
NON_NUCLEOTIDE_PENALTY=1
BAM_OUT_SUFFIX="S288C_sorted"
BLFILTERED_OUT_SUFFIX="blFiltered"

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

# Directory structure
EXPERIMENT_DIR="$HOME/data/${EXPERIMENT_ID}"
FASTQ_DIR="${EXPERIMENT_DIR}/fastq"
ALIGNMENT_DIR="${EXPERIMENT_DIR}/alignment"
COVERAGE_DIR="${EXPERIMENT_DIR}/coverage"
DOCUMENTATION_DIR="${EXPERIMENT_DIR}/documentation"
QUALITY_CONTROL_DIR="${EXPERIMENT_DIR}/quality_control"

# Manifest output configuration
MANIFEST_FILENAME="consolidated_reads_manifest.tsv"
MANIFEST_FILEPATH="${DOCUMENTATION_DIR}/$MANIFEST_FILENAME"

#=======================
# Setup logging
#=======================
# Read in logging functions and then setup
TOOL_NAME="chip_processing"
source "$HOME/lab_utils/core_scripts/functions_for_logging.sh"
# Sets log file in the logs directory. Outputs LOG_DIR.
eval "$(setup_logging "${EXPERIMENT_ID}" ${TOOL_NAME} ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID} )"
log_message "INFO" "=== Starting ChIP-seq processing pipeline ==="

# @QUES: I dont think this is necessary.
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
IFS=$'\t' read -r SAMPLE_ID READ1_PATH READ2_PATH <<< "$MANIFEST_ROW"

# Validate minimum fields
if [[ -z "$SAMPLE_ID" ]] || [[ -z "$READ1_PATH" ]]; then
  log_message "ERROR" "Manifest row must have at least 2 fields: $MANIFEST_ROW"
  exit 1

fi

#log_message "INFO" "Sample ID: $SAMPLE_ID"
#log_message "INFO" "Read 1: $READ1_PATH"
#if [[ -n "$READ2_PATH" ]]; then
  #log_message "INFO" "Read 2: $READ2_PATH"
#fi

#========================
# Detect read type and build tool parameters
#========================

if [[ -z "$READ2_PATH" ]]; then
  # Single-end configuration
  READ_TYPE="single-end"
  log_message "INFO" "Processing as single-end reads"

  # Validate input exists
  if [[ ! -f "$READ1_PATH" ]]; then
      log_message "ERROR" "Read1 file not found: $READ1_PATH"
      exit 1
  fi

  # Define output paths
  FASTP_OUT1="${FASTQ_DIR}/fastpfiltered_${SAMPLE_ID}_R1.fastq"

  # Build tool parameters
  FASTP_INPUT_PARAMS="--in1 ${READ1_PATH}"
  FASTP_OUTPUT_PARAMS="--out1 ${FASTP_OUT1}"
  BOWTIE2_INPUT_PARAMS="-U ${FASTP_OUT1}"

else
  # Paired-end configuration
  READ_TYPE="paired-end"
  log_message "INFO" "Processing as paired-end reads"

  # Validate inputs exist
  if [[ ! -f "$READ1_PATH" ]]; then
      log_message "ERROR" "Read1 file not found: $READ1_PATH"
      exit 1
  fi

  if [[ ! -f "$READ2_PATH" ]]; then
      log_message "ERROR" "Read2 file not found: $READ2_PATH"
      exit 1
  fi

  # Define output paths
  FASTP_OUT1="${FASTQ_DIR}/${FASTP_OUT_PREFIX}_${SAMPLE_ID}_R1.fastq"
  FASTP_OUT2="${FASTQ_DIR}/${FASTP_OUT_PREFIX}_${SAMPLE_ID}_R2.fastq"

  # Build tool parameters
  FASTP_INPUT_PARAMS="--in1 ${READ1_PATH} --in2 ${READ2_PATH}"
  FASTP_OUTPUT_PARAMS="--out1 ${FASTP_OUT1} --out2 ${FASTP_OUT2}"
  BOWTIE2_INPUT_PARAMS="-1 ${FASTP_OUT1} -2 ${FASTP_OUT2}"

fi

#===================================================
# DEFINE OUTPUT FILE PATHS (COMMON FOR BOTH READ TYPES)
#==========================================
# Alignment outputs
SORTED_BAM_OUT_FILE="${ALIGNMENT_DIR}/${SAMPLE_ID}_${BAM_OUT_SUFFIX}.bam"
BLACKLIST_FILTERED_OUT_FILE="${ALIGNMENT_DIR}/${SAMPLE_ID}_${BLFILTERED_OUT_SUFFIX}.bam"
# Coverage output
BIGWIG_OUT_PATH="${COVERAGE_DIR}/${SAMPLE_ID}_${NORM_METHOD}.bw"

# QC outputs
FASTP_JSON="${QUALITY_CONTROL_DIR}/${SAMPLE_ID}_fastp.json"
FASTP_HTML="${QUALITY_CONTROL_DIR}/${SAMPLE_ID}_fastp.html"

# I dont think these are necessary
#log_message "INFO" "Output paths configured"
#log_message "INFO" "  Sorted BAM: $(basename "$SORTED_BAM_OUT_FILE")"
#log_message "INFO" "  Filtered BAM: $(basename "$BLACKLIST_FILTERED_OUT_FILE")"
#log_message "INFO" "  BigWig: $(basename "$BIGWIG_OUT_PATH")"

#===================================================
# Log configuration variables
#==========================================
log_message "CONFIG" "=== Script Configuration ==="
mapfile -t variables_after_config < <(compgen -A variable | sort)

mapfile -t new_variables < <(
    comm -13 <(printf "%s\n" "${initial_vars_declarations[@]}" | sort) \
             <(printf "%s\n" "${variables_after_config[@]}" | sort)
)

# Log new user-defined variables
for variable_name in "${new_variables[@]}"; do
    [[ -z "$variable_name" ]] && continue

    # Skip internal/script variables
    case "$variable_name" in
        initial_vars_declarations|variables_after_config|variable_name|variable_value|formatted_log_line)
        continue
        ;;
    esac

     # Optional: only log ALL_CAPS config-style vars
     if [[ "$variable_name" =~ ^[A-Z_][A-Z0-9_]*$ ]]; then
        printf -v formatted_log_line "%s=%q" "$variable_name" "${!variable_name}"
        log_message "CONFIG" "$formatted_log_line"
     fi

done

log_message "CONFIG" "=== End Configuration ==="

# ============================================================================
# PIPELINE STEP 1: FASTP QUALITY FILTERING
# ============================================================================
log_message "INFO" "=== Step 1: fastp quality filtering ==="

# Check if output already exists
if [[ -f "$FASTP_OUT1" ]] && [[ -f "$FASTP_JSON" ]]; then
    if [[ "$READ_TYPE" == "paired-end" ]] && [[ ! -f "$FASTP_OUT2" ]]; then
        log_message "WARN" "Paired-end run but FASTP_OUT2 missing - will reprocess"

    else
        log_message "SKIP" "fastp outputs already exist - skipping"
        log_message "INFO" "  Output: $FASTP_OUT1"

        if [[ "$READ_TYPE" == "paired-end" ]]; then
            log_message "INFO" "  Output: $FASTP_OUT2"

        fi

    fi

else
    # Load fastp module
    module purge
    module load fastp/0.20.0

    log_message "INFO" "Starting fastp filtering for ${READ_TYPE}"

    # Execute fastp with dynamically built parameters
    # shellcheck disable=SC2086
    if measure_performance "fastp_filtering" \
        fastp \
            ${FASTP_INPUT_PARAMS} \
            ${FASTP_OUTPUT_PARAMS} \
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
            --json "$FASTP_JSON" \
            --html "$FASTP_HTML"; then
        log_message "INFO" "fastp filtering complete"

    else
        log_message "ERROR" "fastp filtering failed for ${SAMPLE_ID}"
        exit 1
    fi
fi

# ============================================================================
# PIPELINE STEP 2: BOWTIE2 ALIGNMENT AND SORTING
# ============================================================================

log_message "INFO" "=== Step 2: bowtie2 alignment and sorting ==="

# Check if output already exists
if [[ -f "$SORTED_BAM_OUT_FILE}.bai" ]]; then
    log_message "SKIP" "Sorted BAM index already exist - skipping"
    log_message "INFO" "  Output: $SORTED_BAM_OUT_FILE"

else
    # Load required modules
    module purge
    module load bowtie2
    module load samtools

    log_message "INFO" "Starting alignment for ${READ_TYPE}"

    # Execute bowtie2 alignment with piped sorting
    # shellcheck disable=SC2086
    if measure_performance "alignment_and_sorting" \
        bowtie2 \
            -x "$GENOME_INDEX" \
            ${BOWTIE2_INPUT_PARAMS} \
            --mp "$MAX_PENALTY_MISMATCH" \
            --np "$NON_NUCLEOTIDE_PENALTY" \
            --sensitive \
            -D "$EXTENSIONS_TO_TRY" \
            -R "$SETS_OF_SEEDS" \
            -k "$ALIGNMENTS_TO_REPORT" \
            -p "$CPU_THREADS" \
            2>> "${ERROR_LOG}" | \
        samtools view -@ "$CPU_THREADS" -bS - 2>> "${ERROR_LOG}" | \
        samtools sort -@ "$CPU_THREADS" -o "$SORTED_BAM_OUT_FILE" - 2>> "${ERROR_LOG}"; then

        log_message "INFO" "Alignment and sorting complete"
        log_message "INFO" "Starting BAM indexing"

        if measure_performance "indexing" samtools index "$SORTED_BAM_OUT_FILE"; then
            log_message "INFO" "BAM indexing complete"

        else
            log_message "ERROR" "BAM indexing failed for ${SAMPLE_ID}"
            exit 1

        fi
    else
        log_message "ERROR" "Alignment/sorting failed for ${SAMPLE_ID}"
        exit 1

    fi

fi


# ============================================================================
# PIPELINE STEP 3: ALIGNMENTSIEVE BLACKLIST FILTERING
# ============================================================================

log_message "INFO" "=== Step 3: alignmentSieve blacklist filtering ==="

# Check if output already exists
if [[ -f "$BLACKLIST_FILTERED_OUT_FILE" ]] && [[ -f "${BLACKLIST_FILTERED_OUT_FILE}.bai" ]]; then
    log_message "SKIP" "Blacklist-filtered BAM and index already exist - skipping"
    log_message "INFO" "  Output: $BLACKLIST_FILTERED_OUT_FILE"

else
    # Load required modules
    module purge
    module load python/2.7.13
    module load deeptools/3.0.1
    module load samtools

    log_message "INFO" "Starting blacklist filtering"

    # Execute alignmentSieve
    if measure_performance "blacklist_filtering" \
        alignmentSieve \
            --bam "$SORTED_BAM_OUT_FILE" \
            --blackListFileName "$BLACKLIST_BED_FILE" \
            --outFile "$BLACKLIST_FILTERED_OUT_FILE" \
            --numberOfProcessors "$((CPU_THREADS / 2))"; then

        log_message "INFO" "Blacklist filtering complete"
        log_message "INFO" "Starting BAM indexing"

        if measure_performance "index" samtools index "$BLACKLIST_FILTERED_OUT_FILE"; then
            log_message "INFO" "Blacklist-filtered BAM indexing complete"

        else
            log_message "ERROR" "Blacklist-filtered BAM indexing failed for ${SAMPLE_ID}"
            exit 1

        fi

    else
        log_message "ERROR" "Blacklist filtering failed for ${SAMPLE_ID}"
        exit 1

    fi

fi

# ============================================================================
# PIPELINE STEP 4: BAMCOVERAGE VISUALIZATION
# ============================================================================

log_message "INFO" "=== Step 4: bamCoverage visualization ==="

# Check if output already exists
if [[ -f "$BIGWIG_OUT_PATH" ]]; then
    log_message "SKIP" "BigWig file already exists - skipping"
    log_message "INFO" "  Output: $BIGWIG_OUT_PATH"

else
    # Load required modules
    module purge
    module load python/2.7.13
    module load deeptools/3.0.1

    log_message "INFO" "Starting bamCoverage"

    # Build bamCoverage parameters
    bamcoverage_params=(
        --bam "$BLACKLIST_FILTERED_OUT_FILE"
        --outFileName "$BIGWIG_OUT_PATH"
        --binSize "$BIN_SIZE"
        --ignoreDuplicates
        --minMappingQuality "$MINIMUM_MAPPING_QUALITY"
        --normalizeUsing "$NORM_METHOD"
        --numberOfProcessors "$CPU_THREADS"
    )

    # Add RPGC-specific parameter if needed
    if [[ "${NORM_METHOD}" == "RPGC" ]]; then
      bamcoverage_params+=(--effectiveGenomeSize "${EFFECTIVE_GENOME_SIZE}")
        log_message "INFO" "Using RPGC normalization with effective genome size: ${EFFECTIVE_GENOME_SIZE}"

    fi

    # Execute bamCoverage
    # shellcheck disable=SC2086
    if measure_performance "bamcoverage" bamCoverage "${bamcoverage_params[@]}"; then
        log_message "INFO" "bamCoverage complete"

    else
        log_message "ERROR" "bamCoverage failed for ${SAMPLE_ID}"
        exit 1

    fi
fi

# ============================================================================
# COMPLETION
# ============================================================================
log_message "INFO" "=== Pipeline complete for ${SAMPLE_ID} ==="
log_message "INFO" "All steps finished successfully"
log_message "INFO" "Final outputs:"
log_message "INFO" "  - Trimmed FASTQ: ${FASTQ_DIR}/fastpfiltered_${SAMPLE_ID}_R*.fastq"
log_message "INFO" "  - Sorted BAM: $SORTED_BAM_FILE"
log_message "INFO" "  - Filtered BAM: $BLACKLIST_FILTERED_OUT_FILE"
log_message "INFO" "  - BigWig: $BIGWIG_PATH"
log_message "INFO" "  - QC reports: ${QUALITY_CONTROL_DIR}/${SAMPLE_ID}_fastp.{json,html}"
