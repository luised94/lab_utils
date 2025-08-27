#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --nice=10000
#SBATCH --exclude=c[5-22]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luised94@mit.edu
# run_fastqc_array.sbatch
# Purpose: Execute FastQC as a SLURM array job
# Version: 1.1.0
# Compatibility: Bash 4.2+, SLURM

# Strict error handling
set -euo pipefail

# Validate input arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <experiment_directory> <fastq_subdirectory>"
    exit 1
fi

# Parse arguments
readonly EXPERIMENT_DIR="$1"
EXPERIMENT_ID=$( basename "$EXPERIMENT_DIR" )
readonly FASTQ_SUBDIR="$2"

# Validate SLURM array job context
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "Error: This script must be run as a SLURM array job"
    exit 1
fi

# Logging functions
source $HOME/lab_utils/core_scripts/functions_for_logging.sh
readonly TOOL_NAME="fastqc"
eval "$(setup_logging "${EXPERIMENT_ID}" ${TOOL_NAME} ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID} )"

readonly QUALITY_CONTROL_DIR="${EXPERIMENT_DIR}/quality_control"
mkdir -p "${QUALITY_CONTROL_DIR}"

# Find FASTQ files
FASTQ_DIR="${EXPERIMENT_DIR}/${FASTQ_SUBDIR}"
mapfile -t FASTQ_FILES < <(find "$FASTQ_DIR" -maxdepth 1 \( -name "*.fastq" -o -name "*.fq" \))
TOTAL_FILES=${#FASTQ_FILES[@]}

# Validate array task
if [[ ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL_FILES} ]]; then
    log_message "ERROR" "Task ID exceeds number of files"
    exit 1
fi

# Select current file
FASTQ_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
FASTQ_PATH="${FASTQ_FILES[${FASTQ_INDEX}]}"
SAMPLE_NAME=$(basename "${FASTQ_PATH}" | sed -E 's/\.(fastq|fq)$//')

# Load required modules
module purge
module load fastqc
module load java

# Execute FastQC with performance measurement
log_message "INFO" "Processing sample: ${SAMPLE_NAME}"
log_message "INFO" "Input file: ${FASTQ_PATH}"

if measure_performance "fastqc" fastqc \
    --outdir "${QUALITY_CONTROL_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK:-1}" \
    "${FASTQ_PATH}"; then
    
    log_message "INFO" "FastQC completed successfully for ${SAMPLE_NAME}"
else
    log_message "ERROR" "FastQC processing failed for ${SAMPLE_NAME}"
    exit 1
fi

# Log completion
log_message "INFO" "Task completed successfully"
