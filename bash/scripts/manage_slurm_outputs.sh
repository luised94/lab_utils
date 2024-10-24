#!/bin/bash
# functions moved

source "../functions/slurm_output_manager.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]
Manage SLURM output files

Options:
    -d, --directory DIR    Search in specific directory
    -o, --organize        Organize files into dated structure
    -c, --cleanup        Remove old log files
    -l, --large          Check for large files
    -h, --help           Show this help message
EOF
}

function main() {
    local search_dir="."
    local do_organize=false
    local do_cleanup=false
    local check_large=false

    while [[ $# -gt 0 ]]; do
        case $1 in
            -d|--directory) search_dir="$2"; shift 2 ;;
            -o|--organize) do_organize=true; shift ;;
            -c|--cleanup) do_cleanup=true; shift ;;
            -l|--large) check_large=true; shift ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done

    if $do_organize; then
        organize_slurm_files "$search_dir"
    fi

    if $do_cleanup; then
        cleanup_old_logs
    fi

    if $check_large; then
        check_large_files "$search_dir"
    fi

    if ! $do_organize && ! $do_cleanup && ! $check_large; then
        find_slurm_files "$search_dir"
    fi
}

main "$@"
