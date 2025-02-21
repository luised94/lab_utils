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
