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
