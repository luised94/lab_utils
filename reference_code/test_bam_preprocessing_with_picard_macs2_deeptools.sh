#!/usr/bin/env bash

set -euo pipefail

# === Cluster Environment Setup ===
#if [[ "$(hostname)" != "luria" ]]; then
#    echo "Error: This script must be run on luria cluster" 1>&2
#    exit 1
#fi

# Load required modules
module load picard/2.18.26 java samtools
export PICARD_JAR="/home/software/picard/picard-2.18.26/picard.jar"

if [[ ! -f "$PICARD_JAR" ]]; then
    echo "[ERROR] Picard JAR missing: $PICARD_JAR" 1>&2
    exit 10
fi

# Initialize Conda and MACS2 environment
source ~/lab_utils/core_scripts/setup_conda_and_macs2.sh || exit 1

##############################################
# Key Adjustable Parameters
##############################################

# Cluster config
THREADS=8

# File paths
declare -A SAMPLES=(
    ['test']="$HOME/data/241010Bel/alignment/processed_245018_sequence_to_S288C_sorted.bam"
    ['input']="$HOME/data/241010Bel/alignment/processed_245003_sequence_to_S288C_sorted.bam" 
    ['reference']="$HOME/data/100303Bel/alignment/processed_034475_sequence_to_S288C_sorted.bam"
)

# Reference data
GENOME_SIZE=12000000  # 1.2e7 in integer form
GENOME_FASTA="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C/SaccharomycescerevisiaeS288C_refgenome.fna"

# Output config
OUTDIR="preprocessing_test"
OUTPUT_PREFIX="test"
SUB_DIRS=("align" "predictd" "peaks" "coverage")

# Processing parameters
PVALUE=1e-6
BIN_SIZE=25
SMOOTH_LEN=75
MIN_FRAGMENT=20
MAX_FRAGMENT=300

##############################################
# Initialization & Validation
##############################################

# Create directory structure
mkdir -p "${SUB_DIRS[@]/#/$OUTDIR/}"

# Display parameter summary
echo "=== Pipeline Configuration ==="
echo "Genome Size: $GENOME_SIZE"
echo "Output Directory: $OUTDIR"
echo "Threads: $THREADS"
echo "Samples:"
for stype in "${!SAMPLES[@]}"; do
    echo " - $stype: ${SAMPLES[$stype]}"
done
echo "=============================="

# File path generators
get_deduped_path() {
    local sample_type=$1
    echo "$OUTDIR/align/${sample_type}_deduped.bam"
}

get_shifted_path() {
    local sample_type=$1
    echo "$OUTDIR/align/${sample_type}_shifted.bam"
}

##############################################
# Preprocessing Workflow
##############################################

# Initialize processed files tracker
declare -A PROCESSED_BAMS=()

# === Step 1: Duplicate Removal ===
echo -e "\n=== Marking Duplicates ==="
for sample_type in "${!SAMPLES[@]}"; do
    input="${SAMPLES[$sample_type]}"
    output=$(get_deduped_path "$sample_type")
    
    echo "Processing $sample_type..."
    echo "Input: $input"
    echo "Output: $output"

    if [[ -f "$output" ]]; then
        echo "Skipping existing: $output"
    else
        java -jar picard.jar MarkDuplicates \
            I="$input" \
            O="$output" \
            M="$OUTDIR/align/${sample_type}_dup_metrics.txt" \
            REMOVE_DUPLICATES=true
    fi

    # Validate output
    if [[ ! -s "$output" ]]; then
        echo "ERROR: Failed to create $output" >&2
        exit 1
    fi
    PROCESSED_BAMS[$sample_type]="$output"
done

# === Step 2: Fragment Analysis ===
echo -e "\n=== Analyzing Fragment Sizes ==="
declare -A FRAGMENTS=()
for sample_type in 'test' 'reference'; do
    input=$(get_deduped_path "$sample_type")
    outdir="$OUTDIR/predictd/${sample_type}_predictd"
    log_file="$outdir/macs2.log"

    echo "Processing $sample_type..."
    mkdir -p "$outdir"

    macs2 predictd \
        -i "$input" \
        -g "$GENOME_SIZE" \
        --outdir "$outdir" 2> "$log_file"

    # Extract fragment size
    frag_size=$(grep -oP 'predicted fragment length is \K\d+' "$log_file")
    
    if [[ -z "$frag_size" ]]; then
        echo "ERROR: Fragment analysis failed for $sample_type" >&2
        exit 2
    fi
    
    FRAGMENTS[$sample_type]=$frag_size
    echo "Fragment Size ($sample_type): ${frag_size}bp"
done

# === Step 3: Read Shifting ===
conda deactivate
echo "Deactivated macs2 conda environment"
if ! command -v macs2 &> /dev/null; then
    echo "MACS2 not found in environment. Deactivation confirmed"
fi

module load python/2.7.13 deeptools/3.0.1
echo "Activated python and deeptools"
echo -e "\n=== Shifting Reads ==="
declare -A SHIFTED_BAMS=()
for sample_type in 'test' 'reference'; do
    input=$(get_deduped_path "$sample_type")
    output=$(get_shifted_path "$sample_type")
    frag_size=${FRAGMENTS[$sample_type]}
    
    echo "Shifting $sample_type by ${frag_size}bp..."
    
    if [[ -f "$output" ]]; then
        echo "Skipping existing: $output"
    else
        alignmentSieve \
            -b "$input" \
            -o "$output" \
            --shift $(($frag_size/2)) -$(($frag_size/2)) \
            --numberOfProcessors "$THREADS"
    fi

    # Post-shift validation
    read_count=$(samtools view -c "$output")
    echo "Shifted BAM contains $read_count reads"
    
    if (( read_count == 0 )); then
        echo "ERROR: Empty shifted BAM: $output" >&2
        exit 3
    fi
    
    SHIFTED_BAMS[$sample_type]="$output"
done

echo -e "\n=== Preprocessing Completed ==="
