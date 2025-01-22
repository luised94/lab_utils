if [[ "$(hostname)" != "luria" ]]; then
    echo "Error: This script must be run on luria cluster"
    exit 1
fi
MACS2_ENV="macs2_env"
OUTDIR="macs2_test_results"
GENOME_SIZE="1.2e7"  # S. cerevisiae genome size
# From 241010Bel, any of the scripts to find appropriate files and set them manually. Should be 18 and 3. Input files were noisy and spiky relative to expected flat input.
TEST_SAMPLE="/home/luised94/data/241010Bel/alignment/processed_245018_sequence_to_S288C_sorted.bam"
INPUT_CONTROL="/home/luised94/data/241010Bel/alignment/processed_245003_sequence_to_S288C_sorted.bam"
# From 100303Bel, file is from 2010 Eaton reference paper.
# Samples are mismatched with the metadata order.
# Correct index for HM1108 antibody is 4.
REF_SAMPLE="/home/luised94/data/100303Bel/alignment/consolidated_034475_sequence_to_S288C_sorted.bam"
OUTPUT_PREFIX="macs2_test"

# For TEST_SAMPLE and INPUT_CONTROL
java -jar picard.jar MarkDuplicates \
  I="$TEST_SAMPLE" \
  O="${TEST_SAMPLE%.bam}_deduped.bam" \
  M=test_dup_metrics.txt \
  REMOVE_DUPLICATES=true

java -jar picard.jar MarkDuplicates \
  I="$INPUT_CONTROL" \
  O="${INPUT_CONTROL%.bam}_deduped.bam" \
  M=input_dup_metrics.txt \
  REMOVE_DUPLICATES=true

# For REF_SAMPLE (no input)
java -jar picard.jar MarkDuplicates \
  I="$REF_SAMPLE" \
  O="${REF_SAMPLE%.bam}_deduped.bam" \
  M=ref_dup_metrics.txt \
  REMOVE_DUPLICATES=true

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
