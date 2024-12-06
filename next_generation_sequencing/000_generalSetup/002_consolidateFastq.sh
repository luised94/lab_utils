#STATUS: REMOVE.
#!/bin/bash
#DESCRIPTION: Consolidate multiple fastq files into one per sample.
#USAGE: ./000_consolidateFastq.sh <experiment_name>
#NOTE: Run inside a node. srun --pty bash

#set -euo pipefail
set -eo pipefail
EXPERIMENT_DIR=""
OUTPUT_DIR=""
FASTQ_PATHS=()
UNIQUE_IDS=()

#Useful error message.
print_usage() {
    echo "Printing usage statement:"
    echo "$0 consolidates fastq files for a given experiment."
    echo "Usage: $0 <experiment_name>"
    echo "Example: $0 \"240808Bel\""
    echo "Dont include trailing slash."
    exit 1
}

#Check that one argument was provided with no trailing slash and see that the argument is a directory.
validate_input(){
    local experiment_name="$1"
    echo "Starting input validation"
    if [ $# -ne 1 ]; then
        echo "Error: Experiment name is required." >&2
        print_usage
    fi
    if [[ "$1" == */ ]]; then
        echo "Error: Please provide the experiment name without a trailing slash." >&2
        echo "Example: 'test' instead of 'test/'" >&2
        exit 1
    fi
    EXPERIMENT_DIR="$HOME/data/$1"
    if [ ! -d "$EXPERIMENT_DIR" ]; then
        echo "Error: Experiment directory not found: $EXPERIMENT_DIR" >&2
        exit 1
    fi
    echo "Input validation complete."
}

#Create output_dir variable.
setup_directories() {
    echo "Setting up directories."
    OUTPUT_DIR="${EXPERIMENT_DIR}/fastq/"
    mkdir -p "$OUTPUT_DIR"
    echo "Directories set up."
}

#Find fastq files to create array. Do not include files with processed or unmapped strings in their name.
get_fastq_files() {
    echo "Gettitng fastq files."
    mapfile -t FASTQ_PATHS < <(find "$EXPERIMENT_DIR" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \) | sort )
    echo "Found ${#FASTQ_PATHS[@]} fastq files."
}

#Use awk to extract unique IDs from the fastq array.
get_unique_ids() {
    echo "Extracting unique IDs."
    #See lmarena thread awkForUniqueIDs to update the awk script.
    UNIQUE_IDS=($(printf '%s\n' "${FASTQ_PATHS[@]}" | awk -F'[_-]' '{print $3}' | sort -u))
    echo "Found ${#UNIQUE_IDS[@]} unique IDs."
}

#Create output file. Then for fastq file, determine if it has the unique ID and cat to the output file.
process_single_id() {
    local unique_id=$1
    local output_file="${OUTPUT_DIR}D24-${unique_id}_NA_sequence.fastq"
    echo "Processing ID: ${unique_id}, Output: ${output_file}"
    for fastq_path in "${FASTQ_PATHS[@]}"; do 
        if [[ $fastq_path =~ $unique_id ]]; then
            echo "Append and delete $fastq_path to $output_file"
            if cat "${fastq_path}" >> "$output_file"; then
                rm "$fastq_path"
                echo "Appended and deleted $fastq_path"
            else 
                echo "Error appending $fastq_path. Not deleted." >&2
                return 1
            fi
        fi
    done
}

#For each unique ID extracted by awk, process the id. See above.
process_fastq_files() {
    echo "Processing fastq files by unique ids."
    for unique_id in "${UNIQUE_IDS[@]}"; do
        process_single_id "$unique_id"
    done
    echo "All unique_ids processed."
}

#Output numbers that will let me know that processing went well.
# @todo: Check documentation folder to output the lines in the tsv file. That will let me confirm that everything went well.
print_summary() {
    echo "Files processed: ${#FASTQ_PATHS[@]}"
    echo "Number of Unique IDs: ${#UNIQUE_IDS[@]}"
    echo "Files in directory: $( find "$EXPERIMENT_DIR" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \) | sort | wc -l )"
    #echo "Number of samples in experiment: $(wc -l < "${EXPERIMENT_DIR}/documentation/${1}_bmc_table.tsv")"
}
main() {
    validate_input "$@"
    setup_directories
    get_fastq_files
    get_unique_ids
    process_fastq_files
    print_summary
    echo "Number of lines in experiment: $( wc -l < ${EXPERIMENT_DIR}/documentation/${1}_bmc_table.tsv)"
    echo "Substract one to account to account for header."
    #echo "Unique ID extracted by awk: ${UNIQUE_IDS[1]}"
}

main "$@"
