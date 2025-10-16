#!/bin/bash

# ============================================
# Configuration
# ============================================
echo "-------Start $0-------"
FASTQ_DIRECTORY="$HOME/data/250930Bel/fastq"
EXPECTED_LANES_PER_SAMPLE=3
EXPECTED_READ_PAIRS_PER_LANE=2  # R1 and R2
EXPECTED_FILES_PER_SAMPLE=$((EXPECTED_LANES_PER_SAMPLE * EXPECTED_READ_PAIRS_PER_LANE))

# Filename parsing indices (split on _ and -)
# Example: 250930Bel_D25-12496-2_1_sequence.fastq
# Parts: [250930Bel, D25, 12496, 2, 1, sequence.fastq]
SAMPLE_ID_START_IDX=1   # "D25"
SAMPLE_ID_END_IDX=2     # "12496"
LANE_IDX=3              # "2"
READ_PAIR_IDX=4         # "1" or "2"

# ============================================
# Discover all fastq files
# ============================================
echo "Looking for fastq files in: $FASTQ_DIRECTORY"
echo "Expected files per sample: $EXPECTED_FILES_PER_SAMPLE"
mapfile -t all_fastq_files < <(find "$FASTQ_DIRECTORY" -type f -name "*.fastq")
echo "Total fastq files found: ${#all_fastq_files[@]}"

# ============================================
# Extract unique sample IDs
# ============================================
declare -A sample_id_map  # Use associative array for uniqueness

for fastq_file in "${all_fastq_files[@]}"; do
  filename=$(basename "$fastq_file")
  echo "Filename: $filename"
  if [[ "$filename" =~ unmapped ]]; then
    echo "Skipping unmapped fastq file"
    continue
  fi

  # Split on _ and - to get components
  IFS='_-' read -ra parts <<< "$filename"
  # Build sample ID from components (D25-12496)
  sample_id="${parts[$SAMPLE_ID_START_IDX]}-${parts[$SAMPLE_ID_END_IDX]}"
  echo "Sample id: $sample_id"
  sample_id_map["$sample_id"]=1

done

# Get sorted list of unique sample IDs
mapfile -t unique_sample_ids < <(printf '%s\n' "${!sample_id_map[@]}" | sort)
echo "Unique samples found: ${#unique_sample_ids[@]}"
echo "Sample IDs:"
printf '  %s\n' "${unique_sample_ids[@]}"

# ============================================
# Verify file counts per sample
# ============================================
echo ""
echo "Verifying file counts per sample..."

declare -A sample_file_lists  # Will store space-separated file lists

for sample_id in "${unique_sample_ids[@]}"; do
  # Find all files for this sample
  mapfile -t files_for_sample < <(
    printf '%s\n' "${all_fastq_files[@]}" | grep "$sample_id" | grep -v "unmapped"
  )

  file_count=${#files_for_sample[@]}

  if [[ $file_count -ne $EXPECTED_FILES_PER_SAMPLE ]]; then
    echo "  [WARNING] $sample_id: found $file_count files, expected $EXPECTED_FILES_PER_SAMPLE"
  else
    echo "  [OK] $sample_id: $file_count files"
  fi

  # Store for later use (joining array into string)
  sample_file_lists["$sample_id"]="${files_for_sample[*]}"
done

# Build pair lists for all samples
declare -A sample_pairs  # Format: "sample_id" -> "R1_path R2_path R1_path R2_path..."

echo "Building pair lists for all samples..."
for sample_id in "${unique_sample_ids[@]}"; do
  echo "Sample: $sample_id"
  IFS=' ' read -ra sample_files <<< "${sample_file_lists[$sample_id]}"
  pairs_for_sample=""

  for lane_num in {1..3}; do
    echo "  Lane $lane_num:"
    r1="" r2=""
    for file in "${sample_files[@]}"; do
      IFS='_-' read -ra parts <<< "$(basename "$file")"
      [[ "${parts[$LANE_IDX]}" == "$lane_num" && "${parts[$READ_PAIR_IDX]}" == "1" ]] && r1="$file"
      [[ "${parts[$LANE_IDX]}" == "$lane_num" && "${parts[$READ_PAIR_IDX]}" == "2" ]] && r2="$file"
    done
    pairs_for_sample="$pairs_for_sample$r1 $r2 "
  done
  
  sample_pairs["$sample_id"]="$pairs_for_sample"
  echo "    Pairs for sample: $(printf "%s\n" $pairs_for_sample) "
done
