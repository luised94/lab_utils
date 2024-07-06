#!/bin/bash

#determine_experiment_id() {}

create_new_experiment() {
    local date_of_creation=$(date +%Y%m%d)
    local filename="${date_of_creation}_experiment.md"
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
    template_file="${SCRIPT_DIR}/templates/experiment_template.md" 
    sed "s/{{EXPERIMENT_NAME}}/trial/g; s/{{DATE}}/${date_of_creation}/g" "$template_file" > "$filename"
    nvim $filename
}
create_new_experiment
