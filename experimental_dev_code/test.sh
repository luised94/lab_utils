#!/bin/bash

# A simple logging function for demonstration
log_message() {
  echo "[LOG] $1"
}

# --- 1. Capture the initial state ---
# Get a list of all currently defined variable names.
# The output of compgen is sorted, which is helpful.
log_message "Capturing initial variable snapshot..."
mapfile -t initial_vars < <(compgen -A variable | sort)
# ===================================================================
# --- 2. Your script's main logic goes here ---
#    (Simulating work by defining a mix of variables)

log_message "Running main script logic..."

# Parse arguments
# EXPERIMENT_DIR is verified during slurm submission. No error handling required.
FASTQ_DIR="./example/fastq"
THREADS=$( nproc )

# Configurations
FASTQ_FILEPATH_PATTERN="consolidated*.fastq"
SAMPLE_ID_BASH_REGEX='^consolidated_(.*)_sequence\.fastq$'

# Additional required files
GENOME_DIR="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C"
GENOME_INDEX="$GENOME_DIR/SaccharomycescerevisiaeS288C_index"
BLACKLIST_BED_FILE="$HOME/data/feature_files/20250423_merged_saccharomyces_cerevisiae_s288c_blacklist.bed"

#SAMPLE_ID_SED_REGEX= "consolidated_\(.*\)_sequence\.fastq"

# Fastp configuration variables
MAXIMUM_UNQUALIFIED_BASE_PERCENT=50
MAXIMUM_N_BASE_COUNT=5
COMPLEXITY_WINDOW_SIZE=4
OVERREPRESENTATION_SAMPLING=50
CPU_THREADS="$SLURM_CPUS_PER_TASK"
MINIMUM_BASE_QUALITY=20
MINIMUM_READ_LENGTH=25

# bowtie2 configuration parameters
# All other settings are default.
# Add <m1>, <m2> and --align-paired-reads
ALIGNMENTS_TO_REPORT=1
EXTENSIONS_TO_TRY=15
SETS_OF_SEEDS=2
MAX_PENALTY_MISMATCH=4
NON_NUCLEOTIDE_PENALTY=1

# bamCoverage configuration properites
# All other settings are default.
BIN_SIZE=10
EFFECTIVE_GENOME_SIZE=12157105
MIN_MAPPING_QUALITY=20
NORM_METHOD="CPM"

BAMCOVERAGE_COMMON_PARAMS="--bam ${BLFILTERED_BAM_PATH} \
    --outFileName ${BAMCOVERAGE_PATH} \
    --binSize ${BIN_SIZE} \
    --ignoreDuplicates \
    --minMappingQuality ${MIN_MAPPING_QUALITY} \
    --normalizeUsing ${NORM_METHOD} \
    --numberOfProcessors ${SLURM_CPUS_PER_TASK}"

# Get current fastq file
PREFILTERED_FASTQ_PATH="potato.fastq"
FASTQ_BASENAME=$(basename "$PREFILTERED_FASTQ_PATH")

# Set the output path files
# Generate output name
FASTP_FILTERED_FASTQ_PATH="${FASTQ_DIR}/fastpfiltered_${SAMPLE_ID}_sequence.fastq"

FASTP_FILTERED_FASTQ_FILENAME=$(basename --suffix=.fastq "$FASTP_FILTERED_FASTQ_PATH" )
PREBLACKLIST_FILTERED_BAM_PATH="${BAM_DIR}/${FASTP_FILTERED_FASTQ_FILENAME}_to_S288C_sorted.bam"

PREBLFILT_BAM_FILENAME=$(basename --suffix=_sorted.bam "$PREBLACKLIST_FILTERED_BAM_PATH" )
BLFILTERED_BAM_PATH="${BAM_DIR}/${PREBLFILT_BAM_FILENAME}_blFiltered.bam"

BLFILTERED_BAM_NAME=$(basename --suffix=.bam "$BLFILTERED_BAM_PATH" )
BAMCOVERAGE_PATH="${EXPERIMENT_DIR}/coverage/${BLFILTERED_BAM_NAME}_${NORM_METHOD}.bw"



# A temporary variable that might be used in a loop
COUNTER=10
# ===================================================================

# --- 3. Find and process the newly created variables ---
log_message "Detecting newly defined variables..."

# Get the final list of all variable names
mapfile -t final_vars < <(compgen -A variable | sort)

# Use comm to find lines that are unique to the final list.
# -1 suppresses lines unique to the first file (initial_vars).
# -3 suppresses lines common to both files.
# The result is only the lines unique to the second file (final_vars).
log_message "The following new variables were detected:"

while IFS= read -r new_var; do
  # Skip empty lines and the 'initial_vars' and 'final_vars' themselves
  if [[ -n "$new_var" && "$new_var" != "initial_vars" && "$new_var" != "final_vars" ]]; then
    # Safely format the output and pass it to your logger
    printf -v formatted_output "%s=%q" "$new_var" "${!new_var}"
    log_message "$formatted_output"
  fi
done < <(comm -13 <(printf "%s\n" "${initial_vars[@]}") <(printf "%s\n" "${final_vars[@]}"))
