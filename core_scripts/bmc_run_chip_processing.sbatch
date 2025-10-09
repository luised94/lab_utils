#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --nice=10000
#SBATCH --exclude=c[5-22]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luised94@mit.edu
# Script: bmc_run_chip_processing.sbatch
# Purpose: Executes processing of fastq files. fastp, bowtie2, bam coverage
# Usage: sbatch --array=1-N%16 $HOME/lab_utils/core_scripts/bmc_run_chip_processing.sbatch <experiment_directory>
# OR RUN via bmc_submit_chip_processing.sh
# Author: Luis
# Date: 2025-10-09

# Validate input arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch --array=1-N%16 $0 <experiment_directory>"
    echo "Example: sbatch --array=1-10%16 $0 /home/user/data/240304Bel"
    echo "Note: Array range should not exceed the number of fastq files"
    exit 1
fi

# Validate SLURM_ARRAY_TASK_ID
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: This script must be run as a SLURM array job"
    echo "Use: sbatch --array=1-N%16 $0 <experiment_directory>"
    exit 1
fi

# Constants
GENOME_DIR="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C"
GENOME_INDEX="$GENOME_DIR/SaccharomycescerevisiaeS288C_index"
# coverage parameters - hardcoded values
BIN_SIZE=10
EFFECTIVE_GENOME_SIZE=12157105
MIN_MAPPING_QUALITY=20
# Normalization methods array
#declare -a NORM_METHODS=("RPKM" "CPM" "BPM" "RPGC")
declare -a NORM_METHODS=("CPM")

BLACKLIST_BED_FILE="$HOME/data/feature_files/20250423_merged_saccharomyces_cerevisiae_s288c_blacklist.bed"

# Build bamCoverage command with specific parameters for each method
COMMON_PARAMS="--bam ${BAM_PATH} \
    --outFileName ${OUTPUT_BIGWIG} \
    --binSize ${BIN_SIZE} \
    --minMappingQuality ${MIN_MAPPING_QUALITY} \
    --ignoreDuplicates \
    --normalizeUsing ${NORM_METHOD} \
    --numberOfProcessors ${SLURM_CPUS_PER_TASK}"

# Ensure genome index exists.
if [ ! -f "${GENOME_INDEX}.1.bt2" ]; then
    log_message "Error: Genome index not found: $GENOME_INDEX"
    exit 1
fi

# Check if blacklist file exists
if [[ ! -f "$BLACKLIST_BED_FILE" ]];
then
    echo "Blacklist file not found: $BLACKLIST_BED_FILE. Blacklist runs will be skipped."
    exit 1
fi

# Parse arguments
EXPERIMENT_DIR="$1"
EXPERIMENT_ID=$( basename "$EXPERIMENT_DIR" )
FASTQ_DIR="${EXPERIMENT_DIR}/fastq"
BAM_DIR="${EXPERIMENT_DIR}/alignment"
THREADS=$( nproc )

# Get current fastq file
FASTQ_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
FASTQ_PATH="${FASTQ_FILES[$FASTQ_INDEX]}"

# Get current BAM file
# Calculate array indices
BAM_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

# Get current BAM file
# Calculate array indices
BAM_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / ${#NORM_METHODS[@]} ))
NORM_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % ${#NORM_METHODS[@]} ))

# Get current BAM file and normalization method
BAM_PATH="${BAM_FILES[$BAM_INDEX]}"
# Get current BAM file and normalization method
BAM_PATH="${BAM_FILES[$BAM_INDEX]}"
NORM_METHOD="${NORM_METHODS[$NORM_INDEX]}"

FASTQ_BASENAME=$(basename "$FASTQ_PATH")
SAMPLE_ID=$(echo "$FASTQ_BASENAME" | sed -n 's/consolidated_\(.*\)_sequence\.fastq/\1/p')
# Generate output name
SAMPLE_NAME=$(basename --suffix=.fastq "$FASTQ_PATH" )
OUTPUT_BAM="${EXPERIMENT_DIR}/alignment/${SAMPLE_NAME}_to_S288C_sorted.bam"
# Set output name
SAMPLE_NAME=$(basename --suffix=_sorted.bam "$BAM_PATH" )
OUTPUT_BAM="${EXPERIMENT_DIR}/alignment/${SAMPLE_NAME}_blFiltered.bam"

# Construct output path
OUTPUT_FASTQ="${EXPERIMENT_DIR}/fastq/processed_${SAMPLE_ID}_sequence.fastq"
OUTPUT_QC="${EXPERIMENT_DIR}/quality_control"
# Set output name
SAMPLE_NAME=$(basename --suffix=.bam "$BAM_PATH" )
OUTPUT_BIGWIG="${EXPERIMENT_DIR}/coverage/${SAMPLE_NAME}_${NORM_METHOD}.bw"

mapfile -t FASTQ_FILES < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "consolidated*.fastq" | sort)
TOTAL_FILES=${#FASTQ_FILES[@]}

## Validate array range
#if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_FILES" ]; then
#    log_message "WARNING" "Task ID ${SLURM_ARRAY_TASK_ID} exceeds number of jobs ${TOTAL_FILES}"
#    exit 1
#fi

if [ ! -f "$FASTQ_PATH" ]; then
    log_message "Error: No fastq file found for index $FASTQ_INDEX"
    log_message "ERROR" "---------------------"
    exit 1
fi

if [[ ! "$FASTQ_BASENAME" =~ ^consolidated_.*_sequence\.fastq$ ]]; then
    log_message "ERROR" "Invalid input filename format: ${FASTQ_BASENAME}"
    log_message "ERROR" "Expected format: consolidate_<ID>_sequence.fastq"
    exit 1
fi

if [ ! -f "$BAM_PATH" ]; then
    log_message "WARNING" "Task ID ${SLURM_ARRAY_TASK_ID} bam path does not exist."
    log_message "WARNING" "Input: ${BAM_PATH}"
    exit 1
fi

if [[ ! "$SAMPLE_ID" =~ ^[0-9]{1,6}$ ]]; then
    log_message "ERROR" "Invalid or missing sample ID from filename: ${FASTQ_BASENAME}"
    exit 1
fi

if [ "$TOTAL_FILES" -eq 0 ]; then
    log_message "ERROR" "No fastq files found in ${FASTQ_DIR}"
    log_message "ERROR" "Check the following:"
    log_message "ERROR" " - Files were consolidated using script 'bmc_consolidate_fastq_by_id.sh'"
    log_message "ERROR" "---------------------"
    exit 1
fi

if [ ! -f "$BAM_PATH" ]; then
    log_message "WARNING" "Task ID ${SLURM_ARRAY_TASK_ID} bam path does not exist."
    log_message "WARNING" "Input: ${BAM_PATH}"
    exit 1
fi
# Read in setup_logging function
source "$HOME/lab_utils/core_scripts/functions_for_logging.sh"
readonly TOOL_NAME="chip_processing"
# Sets log file in the logs directory. Outputs LOG_DIR.
eval "$(setup_logging "${EXPERIMENT_ID}" ${TOOL_NAME} ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID} )"

# Validate output path doesn't already exist (optional)
if [ -f "$OUTPUT_FASTQ" ]; then
    log_message "WARNING" "Output file already exists: ${OUTPUT_FASTQ}"
    log_message "INFO" "Task completed successfully"
    exit 0
fi

mkdir -p "${EXPERIMENT_DIR}/alignment"
mkdir -p "$OUTPUT_QC"
mkdir -p "${EXPERIMENT_DIR}/coverage"

# Load required modules
module purge
module load fastp/0.20.0

# Set quality parameters based on experiment ID
if [[ "$EXPERIMENT_ID" == "100303Bel" ]]; then
    # Legacy experiment from 2010
    MINIMUM_BASE_QUALITY=20
    MINIMUM_READ_LENGTH=20
    log_message "INFO" "Using legacy parameters for 2010 experiment: Q${MINIMUM_BASE_QUALITY}, L${MINIMUM_READ_LENGTH}"
else
    # Modern experiments
    MINIMUM_BASE_QUALITY=20
    MINIMUM_READ_LENGTH=50
    log_message "INFO" "Using standard parameters: Q${MINIMUM_BASE_QUALITY}, L${MINIMUM_READ_LENGTH}"
fi

# Quality Filtering Parameters
readonly MAXIMUM_UNQUALIFIED_BASE_PERCENT=50
readonly MAXIMUM_N_BASE_COUNT=5

# Read Length Parameters
readonly MAXIMUM_READ_LENGTH=150

# Complexity Filtering
readonly COMPLEXITY_WINDOW_SIZE=4
#readonly COMPLEXITY_THRESHOLD=30
readonly OVERREPRESENTATION_SAMPLING=50
readonly DUPLICATION_CALC_ACCURACY=3

# Performance Parameters
#readonly COMPRESSION_LEVEL=0
readonly CPU_THREADS="$SLURM_CPUS_PER_TASK"

# Fastp filtering
log_message "INFO" "Starting fastp filtering"
# Options are not available in fastp 0.20.0 version available in the linux cluster.
#--dedup \
#--dup_calc_accuracy "$DUPLICATION_CALC_ACCURACY" \
#--compression "$COMPRESSION_LEVEL" \
if measure_performance "fastp_filtering" \
    fastp \
        --in1 "$FASTQ_PATH" \
        --out1 "$OUTPUT_FASTQ" \
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
        --json "$OUTPUT_QC/${SAMPLE_ID}_fastp.json" \
        --html "$OUTPUT_QC/${SAMPLE_ID}_fastp.html"; then

    log_message "INFO" "Performed fastp filtering."
else
    log_message "ERROR" "Fastp filtering failed for ${SAMPLE_ID}"
    exit 1
fi
# Run rsync -nav username@luria.mit.edu:~/logs/* ~/logs/ to grab the html files and json.

# Load required modules
module purge
module load bowtie2
module load samtools

if [ -f "${OUTPUT_BAM}" ]; then
  log_message "SKIP" "File already exists: ${OUTPUT_BAM}"
  log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
  exit 0
fi

# Alignment and sorting
# Most sensitive alignment mode
# Report only the best alignment
# Limit seed attempts
# Match bonus
# Maximum penalty for mismatch
# Penalty for ambiguous bases
# Read gap open, extend penalties
# Reference gap open, extend penalties
# Alignment score threshold
#--rdg 5,1 \
#--rfg 5,1 \
#--score-min L,0,-0.3 \
#--ma 2 \
log_message "INFO" "Starting alignment and sorting"
if measure_performance "alignment_and_sorting" \
    bowtie2 -x "$GENOME_INDEX" \
            -U "$FASTQ_PATH" \
            -p "$SLURM_CPUS_PER_TASK" \
            --sensitive \
            -k 1 \
            -D 15 \
            -R 2 \
            --mp 4 \
            --np 1 \
            2>> "${ERROR_LOG}" | \
    samtools view -@ "$SLURM_CPUS_PER_TASK" -bS - 2>> "${ERROR_LOG}" | \
    samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$OUTPUT_BAM" - 2>> "${ERROR_LOG}"; then
    log_message "INFO" "Starting BAM indexing"
    if measure_performance "indexing" samtools index "$OUTPUT_BAM"; then
        log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
    else
        log_message "ERROR" "BAM indexing failed for ${SAMPLE_NAME}"
        exit 1
    fi
else
    log_message "ERROR" "Alignment/sorting failed for ${SAMPLE_NAME}"
    exit 1
fi

# Log completion
log_message "INFO" "Task completed successfully"

# Load required modules
module purge
module load python/2.7.13
module load deeptools/3.0.1
module load samtools

# Execute blacklist filtering ------
# Check if BAM file already exists
if [ -f "$OUTPUT_BAM" ]; then
  log_message "SKIP" "File already exists: ${OUTPUT_BAM}"

  if [ -f "${OUTPUT_BAM}.bai" ]; then
    log_message "SKIP" "File already exists: ${OUTPUT_BAM}.bai"
    exit 0
  fi

  if ! measure_performance "index" samtools index "$OUTPUT_BAM"; then
    log_message "ERROR" "BAM indexing failed for ${SAMPLE_NAME}"
    exit 1
  fi

  log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
  exit 0
fi

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
  log_message "ERROR" "BAM indexing failed for ${SAMPLE_NAME}"
  exit 1
fi
log_message "INFO" "Successfully indexed for ${SAMPLE_NAME}"

# Log completion
log_message "INFO" "Task completed successfully"

# Load required modules
module purge
module load python/2.7.13
module load deeptools/3.0.1

mapfile -t BAM_FILES < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*blFiltered.bam" | sort)
TOTAL_FILES=${#BAM_FILES[@]}
TOTAL_JOBS=$((TOTAL_FILES * ${#NORM_METHODS[@]}))
if [ "$TOTAL_FILES" -eq 0 ]; then
    log_message "ERROR" "No BAM files found in ${BAM_DIR}"
    exit 1
fi

# Add RPGC-specific parameters
if [ "${NORM_METHOD}" == "RPGC" ]; then
    COMMON_PARAMS+=" --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE}"
    log_message "INFO" "Added effectiveGenomeSize parameter for RPGC normalization"
fi
log_message "INFO" "Parameters used:\n"
log_message "INFO" "$COMMON_PARAMS"

if [ -f "$OUTPUT_BIGWIG" ]; then
  log_message "SKIP" "File already exists: ${OUTPUT_BIGWIG}"
  log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
  exit 0
fi

# Execute bamCoverage
log_message "INFO" "Starting bamCoverage processing"
# shellcheck disable=SC2086 # Intentional: $COMMON_PARAMS should be word-split
if ! measure_performance "bamcoverage" bamCoverage $COMMON_PARAMS; then
    log_error "bamCoverage processing failed for ${SAMPLE_NAME}"
    exit 1
fi

log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
# Log completion
log_message "INFO" "Task completed successfully"
