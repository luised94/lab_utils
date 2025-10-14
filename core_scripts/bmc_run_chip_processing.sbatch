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
# From single fastq to bigwig files
# Usage: sbatch --array=1-N%16 $HOME/lab_utils/core_scripts/bmc_run_chip_processing.sbatch <experiment_directory>
# OR RUN via bmc_submit_chip_processing.sh
# To sync html and json files locally, run 'rsync -nav username@luria.mit.edu:~/logs/* ~/logs/'.
# Author: Luis
# Date: 2025-10-09

# Validate input arguments
# Script that submits this script via srun checks that the directory exists.
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


# Capture variables in environment to filter afterwards
mapfile -t initial_vars_declarations < <(compgen -A variable | sort)

# Parse arguments
# EXPERIMENT_DIR is verified during slurm submission. No error handling required.
EXPERIMENT_DIR="$1"
EXPERIMENT_ID=$( basename "$EXPERIMENT_DIR" )
FASTQ_DIR="${EXPERIMENT_DIR}/fastq"
BAM_DIR="${EXPERIMENT_DIR}/alignment"
COVERAGE_DIR="${EXPERIMENT_DIR}/coverage"
QC_DIR="${EXPERIMENT_DIR}/quality_control"
THREADS=$( nproc )

# Read in logging functions and then setup
source "$HOME/lab_utils/core_scripts/functions_for_logging.sh"
TOOL_NAME="chip_processing"
# Sets log file in the logs directory. Outputs LOG_DIR.
eval "$(setup_logging "${EXPERIMENT_ID}" ${TOOL_NAME} ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID} )"

# Configurations
FASTQ_FILEPATH_PATTERN="consolidated*.fastq"
SAMPLE_ID_BASH_REGEX='^consolidated_(.*)_sequence\.fastq$'

# Additional required files
GENOME_DIR="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C"
GENOME_INDEX="$GENOME_DIR/SaccharomycescerevisiaeS288C_index"
BLACKLIST_BED_FILE="$HOME/data/feature_files/20250423_merged_saccharomyces_cerevisiae_s288c_blacklist.bed"

#SAMPLE_ID_SED_REGEX= "consolidated_\(.*\)_sequence\.fastq"

# Fastp configuration variables
MAXIMUM_UNQUALIFIED_BASE_PERCENT=50
MAXIMUM_N_BASE_COUNT=5
COMPLEXITY_WINDOW_SIZE=4
OVERREPRESENTATION_SAMPLING=50
CPU_THREADS="$SLURM_CPUS_PER_TASK"
MINIMUM_BASE_QUALITY=20
MINIMUM_READ_LENGTH=25

# bowtie2 configuration parameters
# All other settings are default.
# Add <m1>, <m2> and --align-paired-reads
ALIGNMENTS_TO_REPORT=1
EXTENSIONS_TO_TRY=15
SETS_OF_SEEDS=2
MAX_PENALTY_MISMATCH=4
NON_NUCLEOTIDE_PENALTY=1

# bamCoverage configuration properites
# All other settings are default.
BIN_SIZE=10
EFFECTIVE_GENOME_SIZE=12157105
MIN_MAPPING_QUALITY=20
NORM_METHOD="CPM"

BAMCOVERAGE_COMMON_PARAMS="--bam ${BLFILTERED_BAM_PATH} \
    --outFileName ${BAMCOVERAGE_PATH} \
    --binSize ${BIN_SIZE} \
    --ignoreDuplicates \
    --minMappingQuality ${MIN_MAPPING_QUALITY} \
    --normalizeUsing ${NORM_METHOD} \
    --numberOfProcessors ${SLURM_CPUS_PER_TASK}"

# Add RPGC-specific parameters
if [ "${NORM_METHOD}" == "RPGC" ]; then
    BAMCOVERAGE_COMMON_PARAMS+=" --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE}"
    log_message "INFO" "Added effectiveGenomeSize parameter for RPGC normalization"
fi

# Ensure directories exist
mkdir -p "$FASTQ_DIR"
mkdir -p "$BAM_DIR"
mkdir -p "$QC_DIR"
mkdir -p "$COVERAGE_DIR"

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
FASTP_FILTERED_FASTQ_PATH="${FASTQ_DIR}/fastpfiltered_${SAMPLE_ID}_sequence.fastq"

FASTP_FILTERED_FASTQ_FILENAME=$(basename --suffix=.fastq "$FASTP_FILTERED_FASTQ_PATH" )
PREBLACKLIST_FILTERED_BAM_PATH="${BAM_DIR}/${FASTP_FILTERED_FASTQ_FILENAME}_to_S288C_sorted.bam"

PREBLFILT_BAM_FILENAME=$(basename --suffix=_sorted.bam "$PREBLACKLIST_FILTERED_BAM_PATH" )
BLFILTERED_BAM_PATH="${BAM_DIR}/${PREBLFILT_BAM_FILENAME}_blFiltered.bam"

BLFILTERED_BAM_NAME=$(basename --suffix=.bam "$BLFILTERED_BAM_PATH" )
BAMCOVERAGE_PATH="${EXPERIMENT_DIR}/coverage/${BLFILTERED_BAM_NAME}_${NORM_METHOD}.bw"

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
    log_message "ERROR" "---------------------"
    exit 1

fi

# REQUIRED FILES
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

if [[ ! -f "$FASTQ_PATH" ]]; then
    log_message "ERROR" "No fastq file found for index $FASTQ_INDEX"
    exit 1

fi

if [[ ! "$SAMPLE_ID" =~ ^[0-9]{1,6}$ ]]; then
    log_message "ERROR" "Invalid or missing sample ID from filename: ${FASTQ_BASENAME}"
    log_message "ERROR" "SAMPLE_ID: ${SAMPLE_ID}"
    exit 1

fi

# Output configuration and variables set
log_message "INFO" "Parameters used:\n"
log_message "INFO" "$BAMCOVERAGE_COMMON_PARAMS"

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
        --json "$QC_DIR/${SAMPLE_ID}_fastp.json" \
        --html "$QC_DIR/${SAMPLE_ID}_fastp.html"; then

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
