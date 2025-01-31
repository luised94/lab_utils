if [[ "$(hostname)" != "luria" ]]; then
    echo "Error: This script must be run on luria cluster"
    exit 1
fi
module load picard
module load java
module load samtools
MACS2_ENV="macs2_env"
OUTDIR="macs2_test_results"
OUTPUT_PREFIX="macs2_test"
# S. cerevisiae genome size
GENOME_SIZE="1.2e7"
# Declare associative array for sample paths. Easy for loop access for validation and unified interface
declare -A SAMPLES=(
    # Core experimental samples
    # 241010Bel is the timecourse arrest-and-release experiment for Chromatin Immunoprecipitation based assay of ORC, its atpase mutant and a suppressor.
    # TEST_SAMPLE
    ['test']="$HOME/data/241010Bel/alignment/processed_245018_sequence_to_S288C_sorted.bam"
    # INPUT_CONTROL
    ['input']="$HOME/data/241010Bel/alignment/processed_245003_sequence_to_S288C_sorted.bam"

    # Reference sample notes:
    # - From 100303Bel (2010 Eaton reference paper)
    # - Original metadata mismatch: Fourth position in sample list contains HM1108 antibody
    # REF_SAMPLE
    ['reference']="$HOME/data/100303Bel/alignment/consolidated_034475_sequence_to_S288C_sorted.bam"
)

# File verification check
missing_files=0
for sample_type in "${!SAMPLES[@]}"; do
    file_path="${SAMPLES[$sample_type]}"
    
    if [[ ! -e "$file_path" ]]; then
        echo "[ERROR] Missing $sample_type sample: $file_path" 1>&2
        missing_files=1
    fi
done

# Exit if any files are missing
if (( missing_files )); then
    echo "Aborting: Required sample files missing. Check paths and permissions." 1>&2
    exit 1
fi

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
        exit 2
    fi
done

echo "All duplicates marked successfully."


# Index all deduped BAMs
samtools index "${TEST_SAMPLE%.bam}_deduped.bam"
samtools index "${INPUT_CONTROL%.bam}_deduped.bam"
samtools index "${REF_SAMPLE%.bam}_deduped.bam"
macs2 predictd \
  -i "${TEST_SAMPLE%.bam}_deduped.bam" \
  -g "$GENOME_SIZE" \
  --outdir "$OUTDIR/predictd_test"

f=$(grep "predicted fragment length" "$OUTDIR/predictd_test/cross_correlation.txt" | awk '{print $NF}')
echo "Fragment size (test): $f bp"

alignmentSieve \
  -b "${TEST_SAMPLE%.bam}_deduped.bam" \
  -o "${TEST_SAMPLE%.bam}_shifted.bam" \
  --shiftToStart \
  --offsetPlus $((f/2)) \
  --offsetMinus -$((f/2)) \
  --numberOfProcessors 8

samtools index "${TEST_SAMPLE%.bam}_shifted.bam"

bamCoverage \
  -b "${TEST_SAMPLE%.bam}_shifted.bam" \
  -o "$OUTDIR/test_sample.bw" \
  --binSize 50 \
  --normalizeUsing RPGC \
  --effectiveGenomeSize "$GENOME_SIZE" \
  --smoothLength 50 \
  --ignoreDuplicates

bamCoverage \
  -b "${INPUT_CONTROL%.bam}_deduped.bam" \
  -o "$OUTDIR/input_control.bw" \
  --binSize 50 \
  --normalizeUsing RPGC \
  --effectiveGenomeSize "$GENOME_SIZE" \
  --smoothLength 50 \
  --ignoreDuplicates

bamCoverage \
  -b "${REF_SAMPLE%.bam}_deduped.bam" \
  -o "$OUTDIR/reference_sample.bw" \
  --binSize 50 \
  --normalizeUsing RPGC \
  --effectiveGenomeSize "$GENOME_SIZE" \
  --smoothLength 50 \
  --ignoreDuplicates

macs2 callpeak \
  -t "${TEST_SAMPLE%.bam}_shifted.bam" \
  -c "${INPUT_CONTROL%.bam}_deduped.bam" \
  -n "${OUTPUT_PREFIX}_test" \
  -g "$GENOME_SIZE" \
  --nomodel \
  --extsize "$f" \          # Single-end: extsize = fragment size
  --pvalue 1e-6 \
  --bdg \
  --outdir "$OUTDIR/test_peaks" \
  --SPMR

macs2 callpeak \
  -t "${REF_SAMPLE%.bam}_deduped.bam" \
  -n "${OUTPUT_PREFIX}_reference" \
  -g "$GENOME_SIZE" \
  --nomodel \
  --extsize "$f" \          # Use the same f as test sample (or estimate separately)
  --pvalue 1e-3 \           # Relax threshold due to no input control
  --nolambda \              # Disables local background (use with caution)
  --bdg \
  --outdir "$OUTDIR/reference_peaks" \
  --SPMR

# Test vs. Input
plotFingerprint \
  -b "${TEST_SAMPLE%.bam}_shifted.bam" "${INPUT_CONTROL%.bam}_deduped.bam" \
  --labels Test Input \
  -o "$OUTDIR/fingerprint_test.png"

# Reference Sample (no input)
plotFingerprint \
  -b "${REF_SAMPLE%.bam}_deduped.bam" \
  --labels Reference \
  -o "$OUTDIR/fingerprint_reference.png"

bedtools intersect \
  -a "$OUTDIR/test_peaks/${OUTPUT_PREFIX}_test_peaks.narrowPeak" \
  -b "$OUTDIR/reference_peaks/${OUTPUT_PREFIX}_reference_peaks.narrowPeak" \
  -wa -wb > "$OUTDIR/peak_overlap.txt"
