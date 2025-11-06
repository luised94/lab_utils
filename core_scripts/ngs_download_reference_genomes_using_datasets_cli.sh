#!/bin/bash
set -euo pipefail

# Show usage information
if [[ "${1:-}" == "-h" ]] || [[ "${1:-}" == "--help" ]]; then
    cat << 'EOF'
USAGE: ./ngs_download_reference_genomes_using_datasets_cli.sh

DESCRIPTION:
    Downloads reference genomes from NCBI and organizes them into a flat
    directory structure with standardized naming.

EXAMPLES:
    ./download_reference_genomes.sh

REQUIREMENTS:
    - NCBI datasets CLI tool (datasets command)
    - Sufficient disk space (estimate 10GB per large genome)
    - Internet connection

OUTPUT:
    Files organized in: ~/data/reference_genomes/
    Naming: Organism_Strain_AccessionSanitized_filetype.ext

EOF
    exit 0
fi
