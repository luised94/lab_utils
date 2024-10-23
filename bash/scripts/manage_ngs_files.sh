#!/usr/bin/env bash
# functions moved

set -o errexit
set -o nounset
set -o pipefail

source "../functions/ngs_file_manager.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]
Manage NGS data files

Options:
    -t, --target DIR     Target directory for file movement
    -s, --source DIR     Source directory to search
    -d, --depth NUM      Maximum search depth (default: ${DEFAULTS[MAX_DEPTH]})
    -b, --batch NUM      Batch size for processing (default: ${DEFAULTS[BATCH_SIZE]})
    -a, --analyze        Only analyze file distribution
    -h, --help          Show this help message
EOF
}

function main() {
    local target_dir=""
    local search_dir="."
    local do_analyze=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -t|--target) target_dir="$2"; shift 2 ;;
            -s|--source) search_dir="$2"; shift 2 ;;
            -d|--depth) DEFAULTS[MAX_DEPTH]="$2"; shift 2 ;;
            -b|--batch) DEFAULTS[BATCH_SIZE]="$2"; shift 2 ;;
            -a|--analyze) do_analyze=true; shift ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done
    
    if $do_analyze; then
        analyze_file_distribution "$search_dir"
    else
        if [[ -z "$target_dir" ]]; then
            log_error "Target directory must be specified"
            show_usage
            exit 1
        fi
        move_ngs_files "$target_dir" "$search_dir"
    fi
}

main "$@"
