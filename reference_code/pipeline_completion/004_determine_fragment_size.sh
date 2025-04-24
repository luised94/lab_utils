
echo -e "\n=== Analyzing Fragment Sizes ==="
declare -A FRAGMENTS=()
for sample_type in "${SELECTED_SAMPLE_KEYS[@]}"; do
  echo "Processing $sample_type..."
  deduplicated_bam="$OUTDIR/align/${sample_type}_deduped.bam"
  outdir="$OUTDIR/predictd/${sample_type}_predictd"
  log_file="$outdir/macs2.log"

  mkdir -p "$outdir"

  # Only run MACS2 if log is missing/incomplete
  if [[ ! -f "$log_file" ]] || ! grep -q 'predicted fragment length is' "$log_file"; then
    echo "Running fragment size prediction..."
    macs2 predictd \
        --ifile "$deduplicated_bam" \
        --mfold 2 200 \
        --gsize "$GENOME_SIZE" \
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
done # end fragment size analysis for loop
