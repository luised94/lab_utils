#!/usr/bin/env bash
#todo: Add blacklist processing
#todo: Add quality control checks
#todo: Create the R script to visualize the results of macs2 testing.

#set -euo pipefail

# === Cluster Environment Setup ===
#if [[ "$(hostname)" != "luria" ]]; then
#    echo "Error: This script must be run on luria cluster" 1>&2
#    exit 1
#fi


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
declare -A CHROM_SIZES=()

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
# Helper functions
##############################################
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

# === File Verification ===
declare -i missing_files=0
declare -A DEDUPED_BAMS=()
declare -A SHIFTED_BAMS=()
for sample_type in "${!SAMPLES[@]}"; do
    input="${SAMPLES[$sample_type]}"
    deduped_path=$(get_deduped_path "$sample_type")
    shifted_path=$(get_shifted_path "$sample_type")
    echo "Processing $sample_type..."
    echo "Input: $input"
    echo "deduped: $deduped_path"
    echo "Shifted: $shifted_path"

    if [[ ! -e "$input" ]]; then
        echo "[ERROR] Missing $sample_type file: $input" 1>&2
        missing_files=1
    fi

    if [[ ! -e "$deduped_path" ]]; then
        echo "Deduplicated file not present. Rerun test_bam_preprocessing_with_picard_macs2_deeptools.sh script."
        missing_files=1
    fi

    if [[ "$sample_type" == "test" || "$sample_type" == "reference" ]]; then
        if [[ ! -e "$shifted_path" ]]; then
            echo "Shifted file not present for sample type '$sample_type'. Rerun test_bam_preprocessing_with_picard_macs2_deeptools.sh script."
            missing_files=1
        fi
    fi
    if [[ ! -e "$shifted_path" ]]; then
        echo "Shifted file not present. Rerun test_bam_preprocessing_with_picard_macs2_deeptools.sh script."
        missing_files=1
    fi

    DEDUPED_BAMS[$sample_type]="$deduped_path"
    SHIFTED_BAMS[$sample_type]="$shifted_path"
done

(( missing_files )) && { echo "Aborting: Missing files"; echo "Run test_bam_preprocessing_with_picard_macs2_deeptools.sh script." ; exit 1; }

echo "All sample files verified successfully."

# === Coverage Track Generation ===
#declare -A COVERAGE_PATHS=(
#    ['test']="${PROCESSED_BAMS[test]}"
#    ['input']="${SAMPLES[input]%.bam}_deduped.bam"
#    ['reference']="${PROCESSED_BAMS[reference]}"
#)
#
#
## Generate multiple normalization versions
#declare -a NORMALIZATION=("" "RPKM" "CPM")  # Raw, RPKM, and CPM
#for sample_type in "${!COVERAGE_PATHS[@]}"; do
#    input_bam="${COVERAGE_PATHS[$sample_type]}"
#    
#    for norm_method in "${NORMALIZATION[@]}"; do
#        # Handle raw (unnormalized) case
#        if [[ -z "$norm_method" ]]; then
#            output_suffix="raw"
#            norm_flag=""
#        else
#            output_suffix="${norm_method,,}"
#            norm_flag="--normalizeUsing $norm_method"
#        fi
#
#        bamCoverage \
#            -b "$input_bam" \
#            -o "$OUTDIR/${sample_type}_${output_suffix}.bw" \
#            $norm_flag \
#            --binSize 25 \
#            --effectiveGenomeSize $GENOME_SIZE \
#            --smoothLength 75 \
#            --ignoreDuplicates \
#            --numberOfProcessors "$THREADS"
#    done
#done
#
## === Peak Calling ===
#declare -A MACS_PARAMS=(
#    ['test']="-t ${PROCESSED_BAMS[test]} -c ${COVERAGE_PATHS[input]}"
#    ['reference']="-t ${PROCESSED_BAMS[reference]}"
#)
#
#for sample_type in "${!MACS_PARAMS[@]}"; do
#    macs2 callpeak \
#        ${MACS_PARAMS[$sample_type]} \
#        -n "${OUTPUT_PREFIX}_${sample_type}" \
#        -g "$GENOME_SIZE" \
#        --nomodel \
#        --extsize "${FRAGMENTS[$sample_type]}" \
#        --pvalue "$PVALUE" \
#        --bdg \
#        --SPMR \
#        --outdir "$OUTDIR/${sample_type}_peaks" \
#        --keep-dup all
#
#    # Verify peak creation
#    peak_file="$OUTDIR/${sample_type}_peaks/${OUTPUT_PREFIX}_${sample_type}_peaks.narrowPeak"
#    [[ -s "$peak_file" ]] || { echo "[ERROR] No peaks found for $sample_type"; exit 3; }
#done
#
## === FRiP Calculation ===
#for sample_type in 'test' 'reference'; do
#    peaks="$OUTDIR/${sample_type}_peaks/${OUTPUT_PREFIX}_${sample_type}_peaks.narrowPeak"
#    reads_in_peaks="$OUTDIR/${sample_type}_peaks/reads_in_peaks.txt"
#    total_reads=$(samtools view -c "${PROCESSED_BAMS[$sample_type]}")
#
#    # Use raw coverage for FRiP calculation
#    bedtools intersect \
#        -a "$peaks" \
#        -b "${PROCESSED_BAMS[$sample_type]}" \
#        -c > "$reads_in_peaks" || { echo "[ERROR] bedtools intersect failed"; exit 4; }
#
#    # Add formatted float output
#    frip=$(awk -v total="$total_reads" '{sum+=$NF} END{printf "%.4f\n", sum/total}' "$reads_in_peaks")
#    echo "${sample_type^^} FRiP: $frip (using raw read counts)" > "$OUTDIR/${sample_type}_peaks/frip_score.txt"
#done
#
## === Quality Control ===
#plotFingerprint \
#    -b "${PROCESSED_BAMS[test]}" "${COVERAGE_PATHS[input]}" \
#    --labels Test Input \
#    -o "$OUTDIR/fingerprint_test_vs_input.png" \
#    --numberOfProcessors "$THREADS"
#
#plotFingerprint \
#    -b "${PROCESSED_BAMS[reference]}" \
#    --labels Reference \
#    -o "$OUTDIR/fingerprint_reference.png" \
#    --numberOfProcessors "$THREADS"
#
## === Peak Comparison ===
#bedtools intersect \
#    -a "$OUTDIR/test_peaks/${OUTPUT_PREFIX}_test_peaks.narrowPeak" \
#    -b "$OUTDIR/reference_peaks/${OUTPUT_PREFIX}_reference_peaks.narrowPeak" \
#    -wa -wb > "$OUTDIR/peak_overlap.tsv"
#
## === Sequence Extraction ===
#for sample_type in 'test' 'reference'; do
#    bedtools getfasta \
#        -fi "$GENOME_FASTA" \
#        -bed "$OUTDIR/${sample_type}_peaks/${OUTPUT_PREFIX}_${sample_type}_peaks.narrowPeak" \
#        -fo "$OUTDIR/${sample_type}_peaks/peak_sequences.fa" \
#        -name
#done
#
#echo "Pipeline completed successfully. Results in: $OUTDIR"
