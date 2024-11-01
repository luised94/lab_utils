#!/bin/bash
# bash/scripts/test_slurm_settings.sh

source "$HOME/lab_utils/bash/functions/logging_utils.sh"

#' Test SLURM Settings
#' @param experiment_dir Character Experiment directory
#' @return Integer 0 if successful
test_slurm_settings_main() {
    if [[ $# -ne 1 ]]; then
        echo "Usage: $0 <experiment_dir>"
        return 1
    }
    
    local experiment_dir="$1"
    local log_file
    log_file=$(initialize_logging "test_slurm")
    
    # Print SLURM Environment
    cat << EOF
=== SLURM Job Information ===
Job ID: ${SLURM_JOB_ID:-Not Set}
Array Job ID: ${SLURM_ARRAY_JOB_ID:-Not Set}
Array Task ID: ${SLURM_ARRAY_TASK_ID:-Not Set}
Job Name: ${SLURM_JOB_NAME:-Not Set}

=== SLURM Resource Allocation ===
Nodes: ${SLURM_JOB_NODELIST:-Not Set}
Node Count: ${SLURM_NNODES:-Not Set}
CPUs per Task: ${SLURM_CPUS_PER_TASK:-Not Set}
Tasks per Node: ${SLURM_NTASKS_PER_NODE:-Not Set}
Memory per Node: ${SLURM_MEM_PER_NODE:-Not Set}

=== Directory Information ===
Working Directory: ${SLURM_SUBMIT_DIR:-Not Set}
Experiment Directory: $experiment_dir

=== Module Status ===
$(module list 2>&1)

=== System Information ===
Hostname: $(hostname)
Current Directory: $(pwd)
User: $USER
Date: $(date)
EOF

    return 0
}

# Execute if run as script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    test_slurm_settings_main "$@"
fi
