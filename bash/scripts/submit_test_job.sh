# bash/scripts/submit_test_job.sh
#!/bin/bash

source "$HOME/lab_utils/bash/functions/slurm_wrapper.sh"

# Calculate array size based on experiment
experiment_dir="241028Bel"
array_range="1-4"  # For testing

# Submit test job
submit_slurm_job_main \
    "$array_range" \
    "test_slurm_settings.sh" \
    "$experiment_dir"
