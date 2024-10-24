#!/bin/bash
# functions moved

source "../functions/experiment_functions.sh"

main() {
    log_info "Starting experiment creation process"
    create_new_experiment
    log_info "Experiment creation completed"
}

main
