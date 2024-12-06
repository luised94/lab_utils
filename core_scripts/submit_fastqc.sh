#!/bin/bash
# submit_fastqc.sh
# Purpose: Submit FastQC processing as a SLURM array job
# Version: 1.0.0
# Compatibility: Bash 4.2+, SLURM

# Strict error handling
set -euo pipefail

# Export scripts directory
export PIPELINE_SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Usage function
usage() {
    echo "Usage: $0 <experiment_directory> [fastq_subdirectory]"
    echo "Example: $0 /path/to/experiment raw_data"
    exit 1
}

# Validate input
EXPERIMENT_DIR="${1:?Error: Experiment directory must be provided}"
FASTQ_SUBDIR="${2:-fastq}"  # Default to 'fastq' if not specified

# Validate experiment directory
if [[ ! -d "${EXPERIMENT_DIR}" ]]; then
    echo "Error: Experiment directory does not exist: ${EXPERIMENT_DIR}"
    exit 1
fi

# Locate FASTQ files
FASTQ_DIR="${EXPERIMENT_DIR}/${FASTQ_SUBDIR}"
if [[ ! -d "${FASTQ_DIR}" ]]; then
    echo "Error: FASTQ directory not found: ${FASTQ_DIR}"
    exit 1
fi

# Count FASTQ files
FASTQ_COUNT=$(find "${FASTQ_DIR}" -maxdepth 1 \( -name "*.fastq" -o -name "*.fq" \) | wc -l)

if [[ ${FASTQ_COUNT} -eq 0 ]]; then
    echo "Error: No FASTQ files found in ${FASTQ_DIR}"
    exit 1
fi

echo "Found ${FASTQ_COUNT} FASTQ files for processing"

# Submit SLURM array job
# Limit to 16 concurrent jobs, adjust as needed
sbatch --array=1-${FASTQ_COUNT}%16 \
    "${PIPELINE_SCRIPTS_DIR}/run_fastqc_array.sh" \
    "${EXPERIMENT_DIR}" \
    "${FASTQ_SUBDIR}"
