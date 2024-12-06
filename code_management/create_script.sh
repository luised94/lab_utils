#STATUS: REMOVE.
#/bin/bash

determine_next_experiment_id() {
    #Find all experiment*.md files, extract second element using _ delimiter. If it is three digit number,
    # select the maximum. if none is found, outputs 001.
    echo $(find . -maxdepth 1 -type f -name "*.sh" | sed 's/\.\///g' | awk -F'_' '
        $1 ~ /^[0-9]{3}$/ {
            if ($1 > max) max = $1
        }
        END {
        printf "%03d", (max == "" ? 1 : max + 1)
        }
    ')
}
# Create the name of the file, find the template location (./templates/ relative to the script location)
# Used sed to replace the tags on the template file with appropriate values. 
#TODO Have to add descriptive name section that reads in name.
#TODO HAve to figure out how to deal with connected files, files with multiple stages, and use template for more specific experiments.
create_new_experiment() {
    local date_of_creation=$(date +%Y%m%d)
    experiment_index=$(determine_next_experiment_id)
    local filename="${date_of_creation}_${experiment_index}_experiment.md"
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

    #template_file="${SCRIPT_DIR}/templates/experiment_template.md" 
    #sed "s/{{EXPERIMENT_NAME}}/trial/g; s/{{DATE}}/${date_of_creation}/g" "$template_file" > "$filename"
    #nvim $filename
    echo $experiment_index
}

create_new_experiment
