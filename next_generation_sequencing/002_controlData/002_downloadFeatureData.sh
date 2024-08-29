#!/bin/bash
#DESCRIPTION: Download feature data from Rossi 2021 paper that will be used to categorical analysis and plot tracking.
#USAGE: $./002_sh_downloadFeatureData.sh

usage() {
    echo "Usage: $0 [download_dir]"
    echo "Downloads Rossi 2021 data from GitHub."
    echo "If download_dir is not specified, defaults to $HOME/data/feature_files"
}

main() {
    local download_dir="${1:-$HOME/data/feature_files}"
    
    echo "Starting Rossi 2021 data download..."
    
    if [ ! -d "$download_dir" ]; then
        echo "Creating directory: $download_dir"
        mkdir -p "$download_dir"
    fi
    
    echo "Cloning Rossi 2021 repository..."
    if git clone --depth=1 https://github.com/CEGRcode/2021-Rossi_Nature.git "$download_dir/rossi_2021"; then
        echo "Successfully cloned Rossi 2021 data to $download_dir/rossi_2021"
    else
        echo "Error: Failed to clone Rossi 2021 repository" >&2
        exit 1
    fi
    
    echo "Rossi 2021 data download complete."
}

# If user provides the -h or --help option, display the usage function output.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
        usage
        exit 0
    fi
    main "$@"
fi
