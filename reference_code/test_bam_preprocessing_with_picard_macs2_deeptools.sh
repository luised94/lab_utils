#!/usr/bin/env bash
echo "Starting script"
#set -euo pipefail

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
CONDA_ROOT=~/miniforge3

if [[ -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
    . "$CONDA_ROOT/etc/profile.d/conda.sh"
else
    echo "ERROR: Conda not found at $CONDA_ROOT" >&2
    exit 1
fi

# Activate MACS2 environment with validation
conda activate macs2_env 2> /dev/null

# Verify MACS2 installation
if ! command -v macs2 &> /dev/null; then
    echo "ERROR: MACS2 not installed. Install with:" >&2
    echo "conda install -c bioconda macs2=2.2.7.1" >&2
    exit 2
fi

echo "MACS2 environment ready in $(which python)"
#source ~/lab_utils/core_scripts/setup_conda_and_macs2.sh || exit 1

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
if [[ ! -f $GENOME_FASTA ]]; then
    echo "$GENOME_FASTA does not exist."
    exit 1
fi

# Output config
OUTDIR="$HOME/preprocessing_test"
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

get_chrom_sizes() {
    local ref_fasta="$GENOME_FASTA"
    local chrom_sizes_file="$OUTDIR/chrom.sizes"
    
    echo "=== Chromosome Size Generation ==="
    
    # Generate .fai if missing
    if [[ ! -f "${ref_fasta}.fai" ]]; then
        echo "Creating FASTA index for reference genome..."
        samtools faidx "$ref_fasta" || {
            echo "[ERROR] Failed to index reference FASTA" >&2
            exit 10
        }
    fi

    # Generate chrom.sizes if missing
    if [[ ! -f "$chrom_sizes_file" ]]; then
        echo -e "\nBuilding chromosome size file..."
        mkdir -p "$OUTDIR"
        
        echo "Input FASTA: $(basename "$ref_fasta")"
        echo "Output sizes: $chrom_sizes_file"
        
        awk '{print $1 "\t" $2}' "${ref_fasta}.fai" > "$chrom_sizes_file" || {
            echo "[ERROR] Failed to create chrom.sizes" >&2
            exit 11
        }
        
        # Validate non-empty output
        [[ -s "$chrom_sizes_file" ]] || {
            echo "[ERROR] chrom.sizes file is empty" >&2
            exit 12
        }
    fi

    # Load into associative array with debug output
    declare -gA CHROM_SIZES=()
    echo -e "\nLoading chromosome sizes:"
    while IFS=$'\t' read -r chrom size; do
        [[ -z "$chrom" ]] && continue  # Skip empty lines
        
        # Trim additional fields after first tab
        chrom="${chrom%%[[:space:]]*}"
        
        echo " - ${chrom}: ${size}bp"
        CHROM_SIZES["$chrom"]="$size"
    done < "$chrom_sizes_file"

    echo -e "\nLoaded ${#CHROM_SIZES[@]} chromosomes"
    echo "======================================="
}

shift_reads() {
    local input=$1
    local output=$2
    local shift_size=$3

    samtools view -h "$input" \
    | awk -v shift="$shift_size" '
        $0 ~ /^@/ {print; next}  # Print header lines unchanged
        {
            if ($2 % 16 == 0)  # Forward strand
                $4 = $4 + shift;
            else  # Reverse strand
                $4 = $4 - shift;
            if ($4 > 0) print;  # Only print if position is valid
        }
    ' \
    | samtools view -bS - > "$output"
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
        java -jar $PICARD_JAR MarkDuplicates \
            I="$input" \
            O="$output" \
            M="$OUTDIR/${sample_type}_dup_metrics.txt" \
            REMOVE_DUPLICATES=true
    fi
    # Validate output
    if [[ ! -s "$output" ]]; then
        echo "ERROR: Failed to create $output" >&2
        exit 2
    fi
    echo "$sample_type Deduped Reads: "
    samtools view -c "$output"
    echo "Inspect first few reads"
    samtools view "$output" | head

    PROCESSED_BAMS[$sample_type]="$output"
done

# === Step 1b: Index Deduplicated BAMs ===
echo -e "\n=== Indexing Deduplicated Files ==="
for sample_type in "${!SAMPLES[@]}"; do
    deduped_file=$(get_deduped_path "$sample_type")
    index_file="${deduped_file}.bai"
    echo "Indexing $sample_type deduped BAM..."
    
    if [[ -f "$index_file" ]]; then
        echo "Index exists: $index_file"
    else
        echo "Indexing: $deduped_file"
        samtools index "$deduped_file"
        exit_code=$?
        
        [[ $exit_code -eq 0 ]] && [[ -f "$index_file" ]] || {
            echo "[ERROR] Indexing failed for ${deduped_file} (exit $exit_code)" >&2
            exit 3
        }
    fi
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

    # Only run MACS2 if log is missing/incomplete
    if [[ ! -f "$log_file" ]] || ! grep -q 'predicted fragment length is' "$log_file"; then
        echo "Running fragment size prediction..."
        macs2 predictd \
            -i "$input" \
            -g "$GENOME_SIZE" \
            --outdir "$outdir" 2> "$log_file"
    else
        echo "Using existing prediction results in: $log_file"
    fi

    # Validate and extract fragment size
    if ! frag_size=$(grep -oP 'predicted fragment length is \K\d+' "$log_file"); then
        echo -e "\nERROR: Fragment analysis failed for $sample_type"
        echo "Debug info:"
        echo "Log file: $log_file"
        [[ -f "$log_file" ]] && tail -n 20 "$log_file"
        exit 4
    fi
    FRAGMENTS[$sample_type]=$frag_size
    echo -e "Fragment Size Analysis Complete\n  Sample: $sample_type\n  Size: ${frag_size}bp\n  Log: ${log_file/$OUTDIR/\$OUTDIR}\n"
done

# === Step 3: Initializing chromosome lengths ===
# In main workflow, after variable initialization
echo -e "\n=== Chromosome Size Initialization ==="
get_chrom_sizes
echo -e "ChromDB status: Loaded ${#CHROM_SIZES[@]} chromosomes\n"

# Example chromosome checks
test_chr="chrI"
echo "Debug: ${test_chr} length = ${CHROM_SIZES[$test_chr]:-UNDEFINED}"

# === Step 4: Read Shifting ===
conda deactivate
echo "Deactivated macs2 conda environment"
if ! command -v macs2 &> /dev/null; then
    echo "MACS2 not found in environment. Deactivation confirmed"
fi

module load python/2.7.13 deeptools/3.0.1
echo "Activated python and deeptools"
echo -e "\n=== Shifting Reads ==="
declare -A SHIFTED_BAMS=()
#for sample_type in 'test' 'reference'; do
#    input=$(get_deduped_path "$sample_type")
#    output=$(get_shifted_path "$sample_type")
#    frag_size=${FRAGMENTS[$sample_type]}
#    shift_size=$((frag_size / 2))
#    echo "Shifting $sample_type by ${shift_size}bp..."
#
#    # Validate input before shifting
#    if [[ ! -s "$input" ]]; then
#        echo "ERROR: Missing input BAM for shifting: $input" >&2
#        exit 4
#    fi
#    
#    if [[ -f "$output" ]]; then
#        echo "Skipping existing: $output"
#    else
#        # Ensure output directory exists
#        mkdir -p "$(dirname "$output")"
#        shift_reads "$input" "$output" "$shift_size"
#        #alignmentSieve \
#        #    -b "$input" \
#        #    -o "$output" \
#        #    -v \
#        #    --shift $(($frag_size/2)) -$(($frag_size/2)) \
#        #    --numberOfProcessors "$THREADS" || {
#        #        echo "Shifting failed. Check:" >&2
#        #        echo "Input: $input (size: $(du -h "$input" | cut -f1))" >&2
#        #        echo "Fragments size used: $frag_size" >&2
#        #        exit 5
#        #    }
#    fi
#    # Post-shift validation
#    if [[ -s "$output" ]]; then
#        reads=$(samtools view -c "$output")
#        echo "Shifted BAM contains $reads reads"
#        SHIFTED_BAMS[$sample_type]="$output"
#    else
#        echo "[CRITICAL] Empty shifted BAM: $output" >&2 
#        exit 5
#    fi
#
#    # Index shifted BAM
#    samtools index "$output" || {
#        echo "ERROR: Failed to index shifted BAM: $output" >&2
#        exit 6
#    }
#done
