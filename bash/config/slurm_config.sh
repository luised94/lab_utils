#!/bin/bash
# bash/config/slurm_config.sh

#' SLURM Configuration Settings
#' @description Centralized SLURM configuration for the project

# Load project configuration
source "$HOME/lab_utils/bash/config/project_config.sh"

#' Core SLURM Settings
declare -A SLURM_CONFIG=(
    # Resource Management
    [NODES]="1"
    [TASKS]="1"
    [CPUS_PER_TASK]="4"
    [MEM_PER_CPU]="50G"
    [MAX_ARRAY_SIZE]="16"
    [ARRAY_FORMAT]="1-%d%%16"
    
    # Job Control
    [NICE]="10000"
    [EXCLUDE_NODES]="c[5-22]"
    
    # Notification
    [MAIL_TYPE]="ALL"
    [MAIL_USER]="luised94@mit.edu"
    
    # Paths
    [LOG_DIR]="slurm_logs"
    [MAX_LOG_AGE]="30"
    
)

#' Module Management
declare -A SLURM_MODULES=(
    # Core Modules
    [GNU]="gnu/5.4.0"
    [SAMTOOLS]="samtools/1.10"
    [BOWTIE2]="bowtie2/2.3.5.1"
    [FASTQC]="fastqc/0.11.5"
    
    # Module Loading Order
    [LOAD_ORDER]="GNU BOWTIE2 SAMTOOLS FASTQC"
)

#' Resource Configurations
declare -A SLURM_RESOURCES=(
    # Alignment Settings
    [ALIGN_MEM]="50G"
    [ALIGN_CPUS]="4"
    [ALIGN_TIME]="24:00:00"
    
    # QC Settings
    [QC_MEM]="20G"
    [QC_CPUS]="2"
    [QC_TIME]="4:00:00"
    
    # BigWig Settings
    [BW_MEM]="30G"
    [BW_CPUS]="2"
    [BW_TIME]="8:00:00"
)

#' Reference Genome Settings
declare -A SLURM_GENOMES=(
    [BASE_DIR]="$HOME/data/REFGENS"
    [DEFAULT_GENOME]="SaccharomycescerevisiaeS288C"
    [GENOME_PATTERN]="*_refgenome.fna"
    [INDEX_SUFFIX]="_index"
)

#' Validation Functions
validate_slurm_config() {
    local log_file="$1"
    
    # Verify required directories
    if [[ ! -d "${SLURM_GENOMES[BASE_DIR]}" ]]; then
        log_error "Reference genome directory not found: ${SLURM_GENOMES[BASE_DIR]}" "$log_file"
        return 1
    fi
    
    # Verify module availability
    for module in ${SLURM_MODULES[LOAD_ORDER]}; do
        if ! module avail "${SLURM_MODULES[$module]}" 2>/dev/null; then
            log_error "Required module not available: ${SLURM_MODULES[$module]}" "$log_file"
            return 1
        fi
    done
    
    return 0
}

#' Generate SLURM Headers
#' @param job_name Character Job name
#' @param config_type Character Resource configuration type
#' @return String SLURM headers
generate_slurm_headers() {
    local job_name="$1"
    local config_type="$2"
    
    cat << EOF
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --nodes=${SLURM_CONFIG[NODES]}
#SBATCH --ntasks=${SLURM_CONFIG[TASKS]}
#SBATCH --cpus-per-task=${SLURM_RESOURCES[${config_type}_CPUS]}
#SBATCH --mem-per-cpu=${SLURM_RESOURCES[${config_type}_MEM]}
#SBATCH --time=${SLURM_RESOURCES[${config_type}_TIME]}
#SBATCH --exclude=${SLURM_CONFIG[EXCLUDE_NODES]}
#SBATCH --nice=${SLURM_CONFIG[NICE]}
#SBATCH --mail-type=${SLURM_CONFIG[MAIL_TYPE]}
#SBATCH --mail-user=${SLURM_CONFIG[MAIL_USER]}

# Load required modules
$(for module in ${SLURM_MODULES[LOAD_ORDER]}; do
    echo "module load ${SLURM_MODULES[$module]}"
done)

EOF
}
