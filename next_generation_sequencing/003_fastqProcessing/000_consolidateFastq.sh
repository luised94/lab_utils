#!/bin/bash
#DESCRIPTION: Consolidate multiple fastq files into one per sample.
#USAGE: ./000_consolidateFastq.sh <experiment_name>

set -euo pipefail
EXPERIMENT_DIR=""
OUTPUT_DIR=""
FASTQ_PATHS=()
UNIQUE_IDS=()

print_usage() {
    echo "Printing usage statement:"
    echo "$0 consolidates fastq files for a given experiment."
    echo "Usage: $0 <experiment_name>"
    echo "Example: $0 \"240808Bel\""
    exit 1
}

validate_input(){
    local experiment_name=$1
    echo "Starting input validation"
    if [ $# -ne 1 ]; then
        echo "Error: Experiment name is required." >&2
        print_usage
    fi
    EXPERIMENT_DIR="$HOME/data/$1"
    if [ ! -d "$EXPERIMENT_DIR" ]; then
        echo "Error: Experiment directory not found: $EXPERIMENT_DIR" >&2
        exit 1
    fi
    echo "Input validation complete."
}

setup_directories() {
    echo "Setting up directories."
    OUTPUT_DIR="${EXPERIMENT_DIR}/fastq/"
    mkdir -p "$OUTPUT_DIR"
    echo "Directories set up."
}

get_fastq_files() {
    echo "Gettitng fastq files."
    mapfile -t FASTQ_PATHS < <(find "$EXPERIMENT_DIR" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \) | sort )
    echo "Found ${#FASTQ_PATHS[@]} fastq files."
}

get_unique_ids() {
    echo "Extracting unique IDs."
    UNIQUE_IDS=($(printf '%s\n' "${FASTQ_PATHS[@]}" | awk -F'-' '{print $2}' | sort -u))
    #UNIQUE_IDS=($(printf '%s\n' "${FASTQ_PATHS[@]}" | awk -F'[_-]' '{print $2}' | sort -u))
    echo "Found ${#UNIQUE_IDS[@]} unique IDs."
}

process_single_id() {
    local unique_id=$1
    local output_file="${EXPERIMENT_DIR}D24-${unique_id}_NA_sequence.fastq"
    echo "Processing ID: ${unique_id}, Output: ${output_file}"
    for fastq_path in "${FASTQ_PATHS[@]}"; do 
        echo "Appending $fastq_path"
        #if [[ $fastq_path =~ $unique_id ]]; then
        #    if cat "${fastq_path}" >> "$output_file"; then
        #        rm "$fastq_path"
        #        echo "Appended and deleted $fastq_path"
        #    else 
        #        echo "Error appending $fastq_path. Not deleted." >&2
        #        return 1
        #    fi
        #fi
    done
}

process_fastq_files() {
    echo "Processing fastq files."
    for unique_id in "${UNIQUE_IDS[@]}"; do
        process_single_id "$unique_id"
    done
    echo "All fastq files processed."
}
print_summary() {
    echo "Files processed: ${#FASTQ_PATHS[@]}"
    echo "Number of Unique IDs: ${#UNIQUE_IDS[@]}"
    echo "Files in directory: $( find "$EXPERIMENT_DIR" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \) | sort | wc -l )"
}
main() {
    validate_input "$@"
    setup_directories
    get_fastq_files
    get_unique_ids
    process_fastq_files
    print_summary
}

main "$@"

#ABSOLUTE_PATH_OF_DIR="$HOME/data/$1"
#OUTPUT_DIR="${ABSOLUTE_PATH_OF_DIR}fastq/"
#mkdir -p "$OUTPUT_DIR"
#
## INITIALIZE_ARRAY
#mapfile -t FASTQ_PATHS < <(find "${ABSOLUTE_PATH_OF_DIR}" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \) | sort )
#
##Could substitute cut with awk -F '-' '{print $2}'
#UNIQUE_IDS=($(printf '%s\n' "${FASTQ_PATHS[@]}" | cut -d- -f2 | uniq ))
#
## Concatenate FASTQ files for each unique ID, filtering paths containing the ID
## and writing to output files named D24-{UNIQUE_ID}_NA_sequence.fastq
#for UNIQUE_ID in ${UNIQUE_IDS[@]}; do
#    OUTPUT_FILE="${OUTPUT_DIR}D24-${UNIQUE_ID}_NA_sequence.fastq"
#    echo "Processing ID: ${UNIQUE_ID}, Output: ${OUTPUT_FILE}"
#    for FASTQ_PATH in "${FASTQ_PATHS[@]}"; do
#        if [[ $FASTQ_PATH =~ $UNIQUE_ID ]]; then 
#            cat "$FASTQ_PATH" >> $OUTPUT_FILE
#            if [ $? -eq 0 ]; then
#                rm "$FASTQ_PATH"
#                echo "Appended and deleted $FASTQ_PATH to $OUTPUT_FILE"
#            else
#                echo "Error appending $FASTQ_PATH. Not deleted"
#            fi
#        fi
#    done
#done
#
#echo "Files processed: $(echo ${#FASTQ_PATHS[@]})"
#echo "Number of Unique IDS: $(echo ${#UNIQUE_IDS[@]})"
#echo "Files in directory: $(find "${ABSOLUTE_PATH_OF_DIR}" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \) | sort | wc -l)"
