#!/bin/bash
# functions moved

set -euo pipefail

source "../functions/fasta_processor.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]
Reformats S288C genome headers to UCSC format

Options:
    -m, --method METHOD  Use specific method (awk|while)
    -n, --no-backup     Skip backup creation
    -r, --restore       Restore from backup
    -h, --help          Show this help message
EOF
}

function find_s288c_genome() {
    local base_dir="$1"
    local pattern="${GENOME_FILES[S288C_PATTERN]}"
    
    log_info "Searching for S288C genome"
    
    local genome_path=$(find "$base_dir" -type f -name "$pattern")
    
    if [ -z "$genome_path" ]; then
        log_error "S288C genome not found"
        return 1
    }
    
    echo "$genome_path"
}

function restore_from_backup() {
    local genome_path="$1"
    local backup_file="${genome_path%_refgenome.fna}${GENOME_FILES[BACKUP_SUFFIX]}"
    
    if [ ! -f "$backup_file" ]; then
        log_error "Backup file not found: $backup_file"
        return 1
    }
    
    log_info "Restoring from backup"
    cp "$backup_file" "$genome_path"
}

function main() {
    local method="awk"
    local do_backup=${FILE_OPERATIONS[BACKUP]}
    local restore=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -m|--method) method="$2"; shift 2 ;;
            -n|--no-backup) do_backup=false; shift ;;
            -r|--restore) restore=true; shift ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done
    
    local genome_path=$(find_s288c_genome "$REFGENOME_DIR") || exit 1
    
    if [ "$restore" = true ]; then
        restore_from_backup "$genome_path"
        exit 0
    fi
    
    validate_fasta "$genome_path" || exit 1
    
    if [ "$do_backup" = true ]; then
        create_backup "$genome_path" || exit 1
    fi
    
    local backup_file="${genome_path%_refgenome.fna}${GENOME_FILES[BACKUP_SUFFIX]}"
    
    if [ "$method" = "awk" ]; then
        reformat_headers_awk "$backup_file" "$genome_path"
    else
        reformat_headers_while "$backup_file" "$genome_path"
    fi
    
    if [ "${FILE_OPERATIONS[VERIFY]}" = true ]; then
        verify_conversion "$backup_file" "$genome_path" || exit 1
    fi
    
    log_info "Header reformatting completed successfully"
}

main "$@"
