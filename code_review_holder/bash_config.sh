
#!/bin/bash

# Add to existing or create new configuration
declare -A BAMCOMPARE_PARAMS=(
    ["BIN_SIZE"]=10
    ["NORMALIZATION"]="CPM"
    ["SCALE_METHOD"]="readCount"
    ["GENOME_SIZE"]=12157105
    ["MIN_MAPPING_QUALITY"]=20
    ["OPERATION"]="ratio"
    ["IGNORE_REGIONS"]=("chrXII")
    ["PSEUDOCOUNT"]=1
)

declare -A OUTPUT_DIRS=(
    ["BIGWIG"]="bigwig"
    ["LOGS"]="logs"
)

declare -A R_CONFIG=(
    ["SCRIPT_PATH"]="$HOME/lab_utils/next_generation_sequencing/004_bamProcessing/determineInputForArrayId.R"
)

# Add to existing module configurations
declare -A REQUIRED_MODULES=(
    ["GNU"]="gnu/5.4.0"
    ["PYTHON"]="python/2.7.13"
    ["DEEPTOOLS"]="deeptools"
    ["R"]="r/4.2.0"
)

#!/bin/bash

declare -A SRA_CONFIG=(
    ["BASE_URL"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
    ["DATA_DIR"]="$HOME/data"
)

# Study-specific configuration
declare -A EATON_2010=(
    ["BIOPROJECT"]="PRJNA117641"
    ["DESCRIPTION"]="ORC precisely positions nucleosomes at origins of replication"
    ["CONDITION"]="WT-G2-ORC"
    ["OUTPUT_FILE"]="nnNnH.fastq"
)

# Sample configuration
declare -A SAMPLES=(
    ["WT-G2-ORC-rep1.fastq.gz"]="SRR034475"
    ["WT-G2-ORC-rep2.fastq.gz"]="SRR034476"
)

#!/bin/bash

# Add to existing QC configurations
declare -A QC_PATHS=(
    ["BASE_DIR"]="$HOME/data"
    ["QC_SUBDIR"]="qualityControl"
    ["MAX_DEPTH"]=1
)

declare -A FILE_PATTERNS=(
    ["FASTQC_ZIP"]="*.zip"
)

declare -A OPERATION_DEFAULTS=(
    ["CONFIRM_TIMEOUT"]=30
    ["UNZIP_BATCH_SIZE"]=10
    ["PRESERVE_ZIP"]=true
)

declare -A BAM_QC_OUTPUTS=(
    ["FLAGSTAT"]="_bamFlagstat.txt"
    ["QUICKCHECK"]="_bamQuickcheck.txt"
    ["STATS"]="_bamStats.txt"
)

declare -A SAMTOOLS_PARAMS=(
    ["FLAGSTAT"]="-O tsv"
    ["STATS"]=""
)

declare -A QC_DIRS=(
    ["OUTPUT"]="qualityControl"
    ["LOGS"]="logs"
)

# Add to existing module configurations
declare -A REQUIRED_MODULES=(
    ["GNU"]="gnu/5.4.0"
    ["SAMTOOLS"]="samtools/1.10"
    ["FASTQC"]="fastqc/0.11.5"
)

#!/bin/bash
# bash/config/project_config.sh

declare -A PROJECT_CONFIG=(
    [REMOTE_HOST]="luria.mit.edu"
    [REMOTE_USER]="luised94"
    [REMOTE_PATH]="$HOME/data"
    [REQUIRED_DIRS]="documentation fastq"
    # File patterns
    [FASTQ_PATTERN]="*.fastq"
    [SAMPLE_PATTERN]="[0-9]{6}Bel"

    # Safety settings
    [REQUIRE_CONFIRMATION]="true"
    [STAGING_DIR_PREFIX]="bmc_staging"
    #
    # FASTQ Processing
    [FASTP_BASE_PARAMS]="--cut_window_size 4 --cut_mean_quality 20 --n_base_limit 5 --average_qual 20 --qualified_quality_phred 20 --unqualified_percent_limit 50 --html /dev/null"
    [FASTP_STANDARD_PARAMS]="--length_required 50"
    [FASTP_EATON_PARAMS]="--length_required 20"
    
    # File Patterns
    [FASTQ_INPUT_PATTERN]="*.fastq"
    [FASTQ_EXCLUDE_PATTERNS]="*unmapped* processed_*"
    [FASTQ_OUTPUT_PREFIX]="processed_"

    # FASTQ Processing
    [FASTQ_PATTERN]="*.fastq"
    [BMC_FASTQ_ID_PATTERN]="[-_]"  # Pattern for splitting
    [BMC_ID_REGEX]="[0-9]{5,6}"    # Pattern for matching ID
    [FASTQ_EXCLUDE]="*unmapped* processed_*"
    [FASTQ_PREFIX]="processed_"
    [FASTQ_SUFFIX]=".fastq"

    # Directory Structure
    [FASTQ_DIR]="fastq"
    [PROCESSED_FASTQ_DIR]="processedFastq"
    [DOC_DIR]="documentation"
    
    # Module Requirements
    [REQUIRED_MODULES]="fastp/0.20.0"
)

#!/bin/bash

declare -A GENOME_CONFIG=(
    ["BASE_DIR"]="$HOME/data/REFGENS"
    ["LOG_DIR"]="$HOME/data/REFGENS/logs"
    ["GENOME_PATTERN"]="*_refgenome.fna"
    ["INDEX_SUFFIX"]="_index"
)

declare -A SLURM_CONFIG=(
    ["NODES"]=1
    ["TASKS"]=1
    ["MEM_PER_CPU"]="20G"
    ["EXCLUDE_NODES"]="c[5-22]"
    ["MAIL_TYPE"]="ALL"
    ["MAIL_USER"]="luised94@mit.edu"
)

declare -A MODULES=(
    ["GNU"]="gnu/5.4.0"
    ["BOWTIE2"]="bowtie2/2.3.5.1"
)

# Add NCBI-specific configurations
declare -A NCBI_CONFIG=(
    ["DOWNLOAD_INCLUDES"]="genome,rna,cds,protein,gff3,gtf"
    ["ACCESSIONS"]=(
        "GCF_000146045.2"  # S. cerevisiae S288C
        "GCF_000001405.40" # Human
        "GCF_000005845.2"  # E. coli
        "GCA_002163515.1"  # S cerevisaie W303
    )
    ["GENOME_PATTERNS"]=("GCF*.fna" "GCA*.fna")
)

declare -A GENOME_NAMING=(
    ["CDS_FILE"]="cds.fna"
    ["REFGENOME_SUFFIX"]="_refgenome.fna"
    ["ORIGINAL_CDS"]="cds_from_genomic.fna"
    ["GENOMIC_PATTERN"]="*_genomic.fna"
)

declare -A GENOME_PATHS=(
    ["NCBI_DATA"]="ncbi_dataset/data"
    ["ASSEMBLY_REPORT"]="assembly_data_report.jsonl"
)

declare -A CHROMOSOME_NAMING=(
    ["PATTERN_OLD"]="chromosome"
    ["PATTERN_NEW"]="chr"
    ["SEPARATOR"]=","
)

# Add to existing genome configurations
declare -A HEADER_PATTERNS=(
    ["OLD_PREFIX"]="chromosome"
    ["NEW_PREFIX"]="chr"
    ["DELIMITER"]=","
)

declare -A GENOME_FILES=(
    ["S288C_PATTERN"]="*S288C_refgenome.fna"
    ["BACKUP_SUFFIX"]="_backup.fna"
)

# Add to existing file patterns if they exist
declare -A FILE_OPERATIONS=(
    ["BACKUP"]=true
    ["VERIFY"]=true
    ["MAX_HEADER_LENGTH"]=100
)

#!/bin/bash

# File patterns and formats
export EXPERIMENT_FILE_PATTERN="*.sh"
export EXPERIMENT_NUMBER_FORMAT="%03d"
export DATE_FORMAT="%Y%m%d"

# File naming
export FILENAME_TEMPLATE="${DATE}_${EXPERIMENT_INDEX}_experiment.md"

# Paths
export TEMPLATE_DIR="templates"
export TEMPLATE_FILE="experiment_template.md"

#!/bin/bash

# Add to existing or create new configuration
declare -A COVERAGE_PARAMS=(
    ["BIN_SIZE"]=10
    ["NORMALIZATION"]="CPM"
    ["MIN_MAPPING_QUALITY"]=20
    ["IGNORE_DUPLICATES"]=true
)
#--normalizeUsing {RPKM,CPM,BPM,RPGC}
declare -A FILE_PATTERNS=(
    ["BAM_SUFFIX"]="S288C.bam"
    ["BIGWIG_SUFFIX"]="_indivNorm.bw"
)

declare -A OUTPUT_DIRS=(
    ["BIGWIG"]="bigwig"
    ["LOGS"]="logs"
)

# Add to existing module configurations
declare -A REQUIRED_MODULES=(
    ["GNU"]="gnu/5.4.0"
    ["PYTHON"]="python/2.7.13"
    ["DEEPTOOLS"]="deeptools"
)

#!/bin/bash
declare -A BMC_CONFIG=(
    [SOURCE_FS]=""
    [TARGET_FS]=""  # More specific
    [RSYNC_OPTIONS]="-av"
    [CLEANUP_DIRS]="*D24* infosite* *_fastqc *_stats"
    [CLEANUP_FILES]="*unmapped*.fastq *.html *.zip"
    [MIN_SPACE_GB]=50
    
)
## Cleanup patterns
## BMC-specific settings
##[BMC_BASE_PATH]="/net/%s/data/bmc/public/Bell/%s"
##[BMC_FASTQ_DIR]="fastq"
##[BMC_DEFAULT_SERVER]="bmc-pub17"

#!/bin/bash
# bash/config/core_config.sh

declare -A CORE_CONFIG=(
    # Version Control
    [VERSION]="1.0.0"
    [MIN_BASH_VERSION]="4.2.0"
    # Logging
    [LOG_LEVELS]="TRACE DEBUG INFO WARNING ERROR FATAL"
    [DEFAULT_LOG_ROOT]="$HOME/logs"
    [LOG_FORMAT]="[%s] [%s] [%s] %s\n"  # timestamp, level, context, message
    # Locking
    [LOCK_TIMEOUT]="30"
    [LOCK_RETRY]="3"
    [LOCK_BASE_DIR]="/tmp/lab_utils_locks"
    # Paths
    [PROJECT_ROOT]="$HOME/lab_utils"
    [MODULE_PATH]="$HOME/lab_utils/bash/modules"
    [VERBOSE]="false"
    [RUN_SEPARATOR]="=== New Run ==="
    [TIMESTAMP_FORMAT]="%Y-%m-%d %H:%M:%S"
    [ENTRY_FORMAT]="\n%s (#%d) === %s ===\n"
    [FIRST_RUN_FORMAT]="%s (#1) === %s ===\n"
    [BUFFER_SIZE]="4096"         # Write buffer size
    [MAX_MESSAGE_LENGTH]="1024"   # Maximum message length
)

#declare -A PATHS=(
#    ["BASE_DIR"]="$HOME/lab_utils"
#    ["FUNCTIONS_DIR"]="bash/functions"
#    ["CONFIG_DIR"]="bash/config"
#    ["SCRIPTS_DIR"]="bash/scripts"
#)
#
#declare -A FILE_PATTERNS=(
#    ["FUNCTIONS"]="*.sh"
#    ["CONFIG"]="*.sh"
#    ["EXCLUDE_PATTERNS"]=(".*" "_*" "test_*")
#)
#
#declare -A LOAD_ORDER=(
#    ["PRIORITY"]=(
#        "logging.sh"
#        "lock_utils.sh"
#        h"
#    )
#    ["OPTIONAL"]=(
#        "experimental.sh"
#        "deprecated.sh"
#    )
#)
#ÃÄ Step marker
#³  Continuation
#ÀÄ Final step
#û Success
#? Failure
