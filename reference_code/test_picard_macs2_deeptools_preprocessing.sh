#!/usr/bin/env bash
#todo: Add blacklist processing
#todo: Add quality control checks
#todo: Create the R script to visualize the results of macs2 testing.

set -euo pipefail

# === Cluster Environment Setup ===
if [[ "$(hostname)" != "luria" ]]; then
    echo "Error: This script must be run on luria cluster" 1>&2
    exit 1
fi

# Load required modules
module load picard java samtools

# Initialize Conda and MACS2 environment
source ~/lab_utils/core_scripts/setup_conda_and_macs2.sh || exit 1

##############################################
# Key Adjustable Parameters
# - PVALUE=1e-6       : Peak calling significance
# - BIN_SIZE=50       : Coverage track resolution
# - SMOOTH_LEN=150    : Signal smoothing window
# - MOTIF_LEN=8       : Expected motif size
##############################################
OUTDIR="macs2_test_results"
OUTPUT_PREFIX="macs2_test"
GENOME_SIZE="1.2e7"
GENOME_FASTA="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C/SaccharomycescerevisiaeS288C_refgenome.fna"
THREADS=8
PVALUE=1e-6

# === Configuration Parameters ===
declare -A SAMPLES=(
    ['test']="$HOME/data/241010Bel/alignment/processed_245018_sequence_to_S288C_sorted.bam"
    ['input']="$HOME/data/241010Bel/alignment/processed_245003_sequence_to_S288C_sorted.bam"
    ['reference']="$HOME/data/100303Bel/alignment/consolidated_034475_sequence_to_S288C_sorted.bam"
)

# === Pipeline Setup ===
mkdir -p "$OUTDIR"/{test_peaks,reference_peaks,predictd_test,predictd_reference}

# === File Verification ===
declare -i missing_files=0
for sample_type in "${!SAMPLES[@]}"; do
    if [[ ! -e "${SAMPLES[$sample_type]}" ]]; then
        echo "[ERROR] Missing $sample_type file: ${SAMPLES[$sample_type]}" 1>&2
        missing_files=1
    fi
done

(( missing_files )) && { echo "Aborting: Missing input files"; exit 1; }
# Continue with rest of script if all files exist
echo "All sample files verified successfully."

# After file verification step
for sample_type in "${!SAMPLES[@]}"; do
    input_file="${SAMPLES[$sample_type]}"
    output_file="${input_file%.bam}_deduped.bam"
    metrics_file="${sample_type}_dup_metrics.txt"

    echo "Processing $sample_type sample..."
    
    java -jar picard.jar MarkDuplicates \
        I="$input_file" \
        O="$output_file" \
        M="$metrics_file" \
        REMOVE_DUPLICATES=true
    
    # Check command success
    if (( $? != 0 )); then
        echo "[ERROR] Failed processing ${sample_type} sample: $input_file" 1>&2
        #exit 2
    fi
done

echo "All duplicates marked successfully."

# After duplicate marking
for sample_type in "${!SAMPLES[@]}"; do
    input_file="${SAMPLES[$sample_type]}"
    deduped_file="${input_file%.bam}_deduped.bam"
    
    # Optional sorting (uncomment if needed)
    # sorted_file="${deduped_file%.bam}_sorted.bam"
    # echo "Sorting $sample_type deduped BAM..."
    # samtools sort "$deduped_file" -o "$sorted_file"
    # deduped_file="$sorted_file"  # Update reference for indexing

    echo "Indexing $sample_type deduped BAM..."
    samtools index "$deduped_file"
    
    # Check indexing success
    if (( $? != 0 )); then
        echo "[ERROR] Failed indexing ${sample_type} sample: $deduped_file" 1>&2
        exit 3
    fi
done


#declare -a FRAGMENT_SAMPLES=('test' 'reference')  # Samples needing fragment analysis
# === Fragment Size Analysis ===
declare -A FRAGMENTS=()
for sample_type in 'test' 'reference'; do
    deduped="${SAMPLES[$sample_type]%.bam}_deduped.bam"
    outdir="$OUTDIR/predictd_${sample_type}"
    cc_file="$outdir/cross_correlation.txt"

    # Run MACS2 predictd if missing results
    [[ -s "$cc_file" ]] || macs2 predictd \
        -i "$deduped" \
        -g "$GENOME_SIZE" \
        --outdir "$outdir"

    # Capture fragment size
    if frag_size=$(awk '/predicted fragment length/{print $NF}' "$cc_file"); then
        FRAGMENTS[$sample_type]=$frag_size
        echo "Fragment size ($sample_type): ${frag_size}bp"
    else
        echo "[ERROR] Failed fragment analysis for $sample_type" 1>&2
        exit 2
    fi
done

# === Read Shifting ===
declare -A PROCESSED_BAMS=()
for sample_type in 'test' 'reference'; do
    input="${SAMPLES[$sample_type]%.bam}_deduped.bam"
    output="${input%.bam}_shifted.bam"
    frag_size="${FRAGMENTS[$sample_type]}"

    alignmentSieve \
        -b "$input" \
        -o "$output" \
        --shiftToStart \
        --offsetPlus $((frag_size/2)) \
        --offsetMinus -$((frag_size/2)) \
        --numberOfProcessors "$THREADS"

    samtools index "$output"
    PROCESSED_BAMS[$sample_type]="$output"
done

# === Coverage Track Generation ===
declare -A COVERAGE_PATHS=(
    ['test']="${PROCESSED_BAMS[test]}"
    ['input']="${SAMPLES[input]%.bam}_deduped.bam"
    ['reference']="${PROCESSED_BAMS[reference]}"
)


# Generate multiple normalization versions
declare -a NORMALIZATION=("" "RPKM" "CPM")  # Raw, RPKM, and CPM

for sample_type in "${!COVERAGE_PATHS[@]}"; do
    input_bam="${COVERAGE_PATHS[$sample_type]}"
    
    for norm_method in "${NORMALIZATION[@]}"; do
        # Handle raw (unnormalized) case
        if [[ -z "$norm_method" ]]; then
            output_suffix="raw"
            norm_flag=""
        else
            output_suffix="${norm_method,,}"
            norm_flag="--normalizeUsing $norm_method"
        fi

        bamCoverage \
            -b "$input_bam" \
            -o "$OUTDIR/${sample_type}_${output_suffix}.bw" \
            $norm_flag \
            --binSize 25 \
            --effectiveGenomeSize "$GENOME_SIZE" \
            --smoothLength 75 \
            --ignoreDuplicates \
            --numberOfProcessors "$THREADS"
    done
done

# === Peak Calling ===
declare -A MACS_PARAMS=(
    ['test']="-t ${PROCESSED_BAMS[test]} -c ${COVERAGE_PATHS[input]}"
    ['reference']="-t ${PROCESSED_BAMS[reference]}"
)

for sample_type in "${!MACS_PARAMS[@]}"; do
    macs2 callpeak \
        ${MACS_PARAMS[$sample_type]} \
        -n "${OUTPUT_PREFIX}_${sample_type}" \
        -g "$GENOME_SIZE" \
        --nomodel \
        --extsize "${FRAGMENTS[$sample_type]}" \
        --pvalue "$PVALUE" \
        --bdg \
        --SPMR \
        --outdir "$OUTDIR/${sample_type}_peaks" \
        --keep-dup all

    # Verify peak creation
    peak_file="$OUTDIR/${sample_type}_peaks/${OUTPUT_PREFIX}_${sample_type}_peaks.narrowPeak"
    [[ -s "$peak_file" ]] || { echo "[ERROR] No peaks found for $sample_type"; exit 3; }
done

# === FRiP Calculation ===
for sample_type in 'test' 'reference'; do
    peaks="$OUTDIR/${sample_type}_peaks/${OUTPUT_PREFIX}_${sample_type}_peaks.narrowPeak"
    reads_in_peaks="$OUTDIR/${sample_type}_peaks/reads_in_peaks.txt"
    total_reads=$(samtools view -c "${PROCESSED_BAMS[$sample_type]}")

    # Use raw coverage for FRiP calculation
    bedtools intersect \
        -a "$peaks" \
        -b "${PROCESSED_BAMS[$sample_type]}" \
        -c > "$reads_in_peaks" || { echo "[ERROR] bedtools intersect failed"; exit 4; }

    # Add formatted float output
    frip=$(awk -v total="$total_reads" '{sum+=$NF} END{printf "%.4f\n", sum/total}' "$reads_in_peaks")
    echo "${sample_type^^} FRiP: $frip (using raw read counts)" > "$OUTDIR/${sample_type}_peaks/frip_score.txt"
done

# === Quality Control ===
plotFingerprint \
    -b "${PROCESSED_BAMS[test]}" "${COVERAGE_PATHS[input]}" \
    --labels Test Input \
    -o "$OUTDIR/fingerprint_test_vs_input.png" \
    --numberOfProcessors "$THREADS"

plotFingerprint \
    -b "${PROCESSED_BAMS[reference]}" \
    --labels Reference \
    -o "$OUTDIR/fingerprint_reference.png" \
    --numberOfProcessors "$THREADS"

# === Peak Comparison ===
bedtools intersect \
    -a "$OUTDIR/test_peaks/${OUTPUT_PREFIX}_test_peaks.narrowPeak" \
    -b "$OUTDIR/reference_peaks/${OUTPUT_PREFIX}_reference_peaks.narrowPeak" \
    -wa -wb > "$OUTDIR/peak_overlap.tsv"

# === Sequence Extraction ===
for sample_type in 'test' 'reference'; do
    bedtools getfasta \
        -fi "$GENOME_FASTA" \
        -bed "$OUTDIR/${sample_type}_peaks/${OUTPUT_PREFIX}_${sample_type}_peaks.narrowPeak" \
        -fo "$OUTDIR/${sample_type}_peaks/peak_sequences.fa" \
        -name
done

echo "Pipeline completed successfully. Results in: $OUTDIR"
