
# === Step 4: Read Shifting ===
echo -e "\n=== Shifting Reads ==="
conda deactivate
echo "Deactivated macs2 conda environment"
if ! command -v macs2 &> /dev/null; then
  echo "MACS2 not found in environment. Deactivation confirmed"
fi

module load python/2.7.13 deeptools/3.0.1
echo "Activated python and deeptools"
# Apply to samples
for sample_type in "${SELECTED_SAMPLE_KEYS[@]}"; do
  deduplicated_bam="$OUTDIR/align/${sample_type}_deduped.bam"
  shifted_bam="$OUTDIR/align/${sample_type}_shifted.bam"
  frag_size=${FRAGMENTS[$sample_type]}
  shift_size=$((frag_size / 2))

  # Debug
  echo -e "\n=== Processing $sample_type ==="
  echo "Deduplicated Bam: $deduplicated_bam"
  echo "Fragment size: $frag_size"
  echo "Shift amount: $shift_size bp"
  echo "Chromosome count: ${#CHROM_SIZES[@]}"

  # Sanity checks
  [[ -n "$frag_size" ]] || { echo "Fragment size missing for $sample_type"; exit 1; }
  [[ "$deduplicated_bam" =~ \.bam$ ]] || { echo "deduplicated_bam must be .bam: $deduplicated_bam" >&2; exit 1; }
  [[ "$shifted_bam" =~ \.bam$ ]] || { echo "Output must be .bam: $shifted_bam" >&2; exit 1; }

  # Check if output already exists and is valid
  if [[ -f "$shifted_bam" ]] && samtools quickcheck -v "$shifted_bam" &> /dev/null; then
    echo "Skipping existing shifted BAM: $shifted_bam"
    continue
  fi
  samtools view -h "$deduplicated_bam" | \
    awk -v shift="$shift_size" -v chrom_file="$OUTDIR/chrom.sizes" -v log_file="$OUTDIR/shift_stats.log" -f "$HOME/lab_utils/core_scripts/shift_reads.awk" | \
    samtools sort -@ $THREADS | \
    samtools view -bS - > "$shifted_bam" || {
      echo "Shifting failed for $deduplicated_bam" >&2
      rm -f "$shifted_bam"
      return 1
    }
  # Cleanup and validation
  #rm "$chrom_sizes_file"
  samtools index "$shifted_bam"
  echo "Shift validation:"
  samtools idxstats "$shifted_bam"

  if [[ -f "$shifted_bam" ]]; then
    echo "Shift Complete: $(basename "$shifted_bam")"
  else
    echo "Shift failed for $shifted_bam."
  fi
done # end bam shift for loop

