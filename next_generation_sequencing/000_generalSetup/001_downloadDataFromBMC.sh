#!/bin/bash
#DESCRIPTION: Download the fastq data from BMC servers. 
#USAGE: ./001_downloadDataFromBMC.sh <bmc_server> <experiment_directory>
#NOTES: Each function contains some form of validation. 
validate_input() {
    echo "Executing validate_input"
    if [ $# -ne 2 ]; then 
        echo "Usage: $0 <bmc_server> <experiment_directory>"
        echo "Example: $0 bmc_pub17 240808Bel"
        exit 1
    fi
}

download_files() {
    local bmc_server=$1
    local bell_lab_directory=$2
    local target_dir="$HOME/data/${bell_lab_directory}/fastq/"

    rsync_command="srun rsync -av --include '*/' --include '*.fastq' --exclude '*' /net/${bmc_server}/data/bmc/public/Bell/${bell_lab_directory}/ ${target_dir}"
    echo "Executing: ${rsync_command}"
    if ! eval ${rsync_command}; then 
        echo "Rysnc failed. Check your parameters or connection"
        exit 1
    fi
}

organize_files() {
    echo "Executing organize_files"
    local target_dir="$HOME/data/${bell_lab_directory}/fastq/"

    if [ ! -d "$target_dir" ]; then
        echo "$target_dir"
        echo "Target directory doesnt exist."
        echo "Run 000"
        exit 1
    fi
    cd "$target_dir" || exit 1
    find . -type f -name "*.fastq" -exec mv {} . \;
}

safe_remove_dirs() {
    echo "Executing safe_remove_dirs"
    local pattern=$1
    local dirs_to_remove=$(find . -type d -name "$pattern")
    if [ -n "$dirs_to_remove" ]; then
        echo "The following directories mathc the pattern ${pattern}:"
        echo "$dirs_to_remove"
        read -p "DELETE these directories? (y/n) " -n 1 -r 
        echo
        if [[ $REPLY =~ ^[YY]$ ]]; then
            local staging_dir="tmp/staging_deleter_$$"
            mkdir -p "$staging_dir"
            echo "Moving directories to staging area"
            echo "$dirs_to_remove" | xargs -I {} mv {} "$staging_dir/"
            echo "Removing directories from staging area"
            rm -rf "$staging_dir"
            echo "Deletion complete"
        else
            echo "Deletion aborted"
        fi
    else
        echo "No directories found matching pattern ${pattern}"
    fi
}

safe_remove_files() {
    echo "Executing safe_remove_files"
    local pattern=$1
    local files_to_remove=$(ls $pattern 2>/dev/null)
    if [ -n "$files_to_remove" ]; then
        echo "The following files match pattern ${pattern}: "
        echo "$files_to_remove"
        read -p "DELETE these files? (y/n) " -n 1 -r
        echo 
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            local staging_dir="tmp/staging_deleter_$$"
            mkdir -p "$staging_dir"
            echo "Moving files to staging area"
            echo "$files_to_remove" | xargs -I {} mv {} "$staging_dir/"
            echo "Removing files from staging area"
            rm -rf "$staging_dir"
            echo "Deletion complete"
        else
            echo "Deleteion aborted"
        fi
    else
        echo "No files found matching ${pattern}"
    fi
}
cleanup() {
    echo "Executing cleanup"
    safe_remove_dirs "*D24*"
    safe_remove_dirs "infosite*"
    safe_remove_files "*unmapped*.fastq"
}

main() {
    validate_params "$@"
    local bmc_server=$1
    local bell_lab_directory=$2
    download_files "$bmc_server" "$bell_lab_directory"
    organize_files
    cleanup
}

main "$@"

