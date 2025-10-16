#!/bin/bash

# ============================================
# Configuration
# ============================================
echo "-------Start $0-------"
echo "Setting up configuration..."
FASTQ_LINES_PER_READ=4
FASTQ_DIRECTORY="$HOME/data/250930Bel/fastq"
# Output to documentation directory (parent of fastq folder)
documentation_dir="$(dirname "$FASTQ_DIRECTORY")/documentation"
manifest_file="$documentation_dir/paired_reads_manifest.tsv"

EXPECTED_LANES_PER_SAMPLE=3
EXPECTED_READ_PAIRS_PER_LANE=2  # R1 and R2
EXPECTED_FILES_PER_SAMPLE=$((EXPECTED_LANES_PER_SAMPLE * EXPECTED_READ_PAIRS_PER_LANE))
# Lanes in sequencer.
MIN_LANE_NUMBER=1
MAX_LANE_NUMBER=4

# Filename parsing indices (split on _ and -)
# Example: 250930Bel_D25-12496-2_1_sequence.fastq
# Parts: [250930Bel, D25, 12496, 2, 1, sequence.fastq]
SAMPLE_ID_START_IDX=1   # "D25"
SAMPLE_ID_END_IDX=2     # "12496"
LANE_IDX=3              # "2"
READ_PAIR_IDX=4         # "1" or "2"

# Ensure documentation_dir exists
mkdir -p "$documentation_dir"
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

  for lane_num in $(seq $MIN_LANE_NUMBER $MAX_LANE_NUMBER); do
    echo "  Processing lane $lane_num"
    r1="" r2=""
    for file in "${sample_files[@]}"; do
      IFS='_-' read -ra parts <<< "$(basename "$file")"
      [[ "${parts[$LANE_IDX]}" == "$lane_num" && "${parts[$READ_PAIR_IDX]}" == "1" ]] && r1="$file"
      [[ "${parts[$LANE_IDX]}" == "$lane_num" && "${parts[$READ_PAIR_IDX]}" == "2" ]] && r2="$file"
    done
    pairs_for_sample="$pairs_for_sample$r1 $r2 "
  done

  sample_pairs["$sample_id"]="$pairs_for_sample"
  echo "    Pairs for sample: "
  printf "    %s\n" $pairs_for_sample
done


echo ""
#echo "Verifying read counts in pairs..."
#
#for sample_id in "${unique_sample_ids[@]}"; do
#  echo "Sample: $sample_id"
#  IFS=' ' read -ra pair_files <<< "${sample_pairs[$sample_id]}"
#  
#  # Process pairs (every 2 files = one pair)
#  for ((i=0; i<${#pair_files[@]}; i+=2)); do
#    r1_file="${pair_files[$i]}"
#    r2_file="${pair_files[$((i+1))]}"
#
#    # Skip if files are empty strings
#    [[ -z "$r1_file" || -z "$r2_file" ]] && continue
#
#    r1_lines=$(wc -l < "$r1_file")
#    r2_lines=$(wc -l < "$r2_file")
#    r1_reads=$((r1_lines / FASTQ_LINES_PER_READ))
#    r2_reads=$((r2_lines / FASTQ_LINES_PER_READ))
#
#    if [[ $r1_reads -eq $r2_reads ]]; then
#      echo "  [OK] $(basename "$r1_file") <-> $(basename "$r2_file"): $r1_reads reads"
#    else
#      echo "  [ERROR] Read count mismatch: R1=$r1_reads, R2=$r2_reads"
#    fi
#  done
#done


echo "Writing paired reads manifest to: $manifest_file"

# Write header and data rows (sample_id, lane, R1_path, R2_path)
echo -e "sample_id\tlane\tread1_path\tread2_path" > "$manifest_file"

for sample_id in "${unique_sample_ids[@]}"; do
  IFS=' ' read -ra pair_files <<< "${sample_pairs[$sample_id]}"
  
  lane_num=1
  for ((i=0; i<${#pair_files[@]}; i+=2)); do
    r1="${pair_files[$i]}"
    r2="${pair_files[$((i+1))]}"
    [[ -z "$r1" || -z "$r2" ]] && continue  # Skip empty pairs
    
    echo -e "$sample_id\t$lane_num\t$r1\t$r2" >> "$manifest_file"
    ((lane_num++))
  done
done

echo "Manifest complete: $(wc -l < "$manifest_file") pairs written"
echo "Other scripts can detect paired-end mode by checking file existence"
