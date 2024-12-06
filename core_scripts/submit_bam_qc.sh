#!/bin/bash
# submit_bam_qc.sh
# Purpose: Submit BAM quality control processing as SLURM array job
# Version: 1.0.0
# Compatibility: Bash 4.2+, SLURM

# Strict error handling
set -euo pipefail

# Export scripts directory
export PIPELINE_SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Usage function
usage() {
    echo "Usage: $0 <experiment_directory> [bam_subdirectory]"
    echo "Example: $0 /path/to/experiment alignment"
    exit 1
}

# Validate input
EXPERIMENT_DIR="${1:?Error: Experiment directory must be provided}"
BAM_SUBDIR="${2:-alignment}"  # Default to 'alignment' if not specified

# Validate experiment directory
if [[ ! -d "${EXPERIMENT_DIR}" ]]; then
    echo "Error: Experiment directory does not exist: ${EXPERIMENT_DIR}"
    exit 1
fi

# Locate BAM files
BAM_DIR="${EXPERIMENT_DIR}/${BAM_SUBDIR}"
if [[ ! -d "${BAM_DIR}" ]]; then
    echo "Error: BAM directory not found: ${BAM_DIR}"
    exit 1
fi

# Count BAM files (sorted and indexed)
BAM_COUNT=$(find "${BAM_DIR}" -maxdepth 1 \( -name "*.sorted.bam" -o -name "*.bam" \) | wc -l)

if [[ ${BAM_COUNT} -eq 0 ]]; then
    echo "Error: No BAM files found in ${BAM_DIR}"
    exit 1
fi

echo "Found ${BAM_COUNT} BAM files for quality control"

# Submit SLURM array job
# Limit to 16 concurrent jobs, adjust as needed
sbatch --array=1-${BAM_COUNT}%16 \
    "${PIPELINE_SCRIPTS_DIR}/run_bam_qc_array.sh" \
    "${EXPERIMENT_DIR}" \
    "${BAM_SUBDIR}"
