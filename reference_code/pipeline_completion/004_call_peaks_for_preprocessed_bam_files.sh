#!/usr/bin/env bash
################################################################################
# Call peaks for raw, deduplicated and the shifted bam
################################################################################
# Purpose: Generate peak calling files for preprocessed bam files
# Usage: Run as script.
# $./003_deduplicate_and_shift_bam_files.sh
# DEPENDENCIES: bash, data from 250207Bel BMC experiment and scripts 001 and 002 from the pipeline script
# OUTPUT: Duplicate files with renaming for easier identification and isolation from data directories
# NOTES: Updated the files being used in analysis, went from the data in the 241010Bel to the data in 250207Bel
# AUTHOR: LEMR
# DATE: 2025-02-07
# UPDATE: 2025-04-22
################################################################################
#set -euo pipefail

# === Cluster Environment Setup ===
#if [[ "$(hostname)" != "luria" ]]; then
#    echo "Error: This script must be run on luria cluster" 1>&2
#    exit 1
#fi

# Genome data
GENOME_SIZE=12000000  # 1.2e7 in integer form
GENOME_FASTA="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C/SaccharomycescerevisiaeS288C_refgenome.fna"
if [[ ! -f $GENOME_FASTA ]]; then
  echo "$GENOME_FASTA does not exist."
  exit 1
fi

# Core parameters used in all calls
CORE_PARAMS=(
  "-g $GENOME_SIZE"
  "--SPMR"
)

# Analysis variants with significance thresholds
declare -A PEAK_MODES=(
  ['narrow_p']="--nomodel --call-summits --pvalue 1e-6"
  ['narrow_q']="--nomodel --call-summits --qvalue 0.01"
  ['broad_p']="--broad --broad-cutoff 0.1 --nomodel --pvalue 1e-6"
  ['broad_q']="--broad --broad-cutoff 0.1 --nomodel --qvalue 0.01"
  ['auto_p']="--fix-bimodal --call-summits --pvalue 1e-6"
  ['auto_q']="--fix-bimodal --call-summits --qvalue 0.01"
)

# Output config
OUTDIR="$HOME/data/preprocessing_test"
SUB_DIRS=("align" "predictd" "peaks" "coverage")
# Create directory structure
mkdir -p "${SUB_DIRS[@]/#/$OUTDIR/}"

# Initialization of array with keys and file paths
# Create an array of files
mapfile -t FILES < <(find ~/data/preprocessing_test/align -maxdepth 1 -type f -name "*.bam" | sort)
# Count unique sample types (excluding inputs)
number_of_sample_types=$(basename -a "${FILES[@]}" | cut -d_ -f1-2 | grep -v "input" | sort -u | wc -l)

# Check if any files were found
if [ ${#FILES[@]} -eq 0 ]; then
  echo "No *.bam files found in ~/data/preprocessing_test/align" >&2
  exit 1
fi

# Initialize Conda
CONDA_ROOT=~/miniforge3
if [[ -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
  . "$CONDA_ROOT/etc/profile.d/conda.sh"
else
  echo "ERROR: Conda not found at $CONDA_ROOT" >&2
  exit 1
fi

# Activate MACS2 environment
echo "Activating MACS2 environment..."
conda activate macs2_env 2> /dev/null

# Verify MACS2 installation
if ! command -v macs2 &> /dev/null; then
  echo "ERROR: MACS2 not installed. Install with:" >&2
  echo "Run: conda install -c bioconda macs2=2.2.7.1" >&2
  exit 2
fi

echo "MACS2 environment ready in $(which python)"
# === Load fragment samples for sample types ===
# Read log file paths into array
mapfile -t FRAG_LOGS < <(find "$OUTDIR/predictd/" -type f -name "macs2.log")

# Check if any logs were found
if [[ ${#FRAG_LOGS[@]} -eq 0 ]]; then
  echo "[ERROR] Failed to find fragment size logs in $OUTDIR/predictd/" >&2
  exit 1
fi


# Validate number of log files
if [[ ${#FRAG_LOGS[@]} -ne $number_of_sample_types ]]; then
  echo "[ERROR] Expected $number_of_sample_types fragment logs, found ${#FRAG_LOGS[@]}" >&2
  echo "Found logs:" >&2
  printf '%s\n' "${FRAG_LOGS[@]}" >&2
  exit 2
fi

# Initialize fragments associative array
declare -A FRAGMENT_SIZES=()
for log in "${FRAG_LOGS[@]}"; do
  # Extract sample type from path
  sample_name=$(basename "$(dirname "$log")" | sed 's/_predictd//')

  # Extract fragment size
  frag_size=$(grep -oP 'predicted fragment length is \K\d+' "$log")

  # Check if extraction succeeded
  if [[ -z "$frag_size" ]]; then
    echo "[ERROR] Failed to extract fragment size from $log" >&2
    exit 3
  fi

  echo "Sample type: $sample_name"
  echo "Fragment size: $frag_size"

  # Validate fragment size
  if ! [[ "$frag_size" =~ ^[0-9]+$ ]] || (( frag_size <= 0 )); then
    echo "[ERROR] Invalid fragment size ($frag_size) in $log" >&2
    exit 4
  fi

  # Store in associative array
  FRAGMENT_SIZES[$sample_name]=$frag_size
done

# Find input control
INPUT_CONTROL=""
for filepath in "${FILES[@]}"; do
  filename=$(basename "$filepath")
  if [[ "$filename" == input_*_deduped.bam ]]; then
    INPUT_CONTROL="$filepath"
    echo "Found input control file: $INPUT_CONTROL"
    break
  fi
done

# Process the files if found
for filepath in "${FILES[@]}"; do
  filename=$(basename "$filepath")
  key=${filename%.bam}
  sample_type=$(echo "$key" | cut -d'_' -f1)
  sample_name=$(echo "$key" | cut -d'_' -f1-2)
  bam_type=$(echo "$key" | cut -d'_' -f3 )
  fragment_size=${FRAGMENT_SIZES[$sample_name]]}

  echo "---Sample information---"
  echo "  Filename: $filename"
  echo "  Sample_name: $sample_name"
  echo "  Sample type: $sample_type"
  echo "  Bam type: $bam_type"
  echo "  Fragment_size: $fragment_size"

  if [[ "$sample_type" == "input" && "$bam_type" == "shifted" ]]; then
    echo "Skipping input sample for shifted analysis"
    continue
  fi

done
exit 1 # Breakpoint

get_peak_name() {
  local sample=$1
  local bam_type=$2
  local mode=$3
  local has_input=$4
  echo "${OUTDIR}/peaks/${sample}_${bam_type}_${mode}_${has_input}"
}

process_bam_set() {
  local array_name=$1   # Name of the array to process
  local bam_type=$2     # raw/deduped/shifted
  local allow_input=$3  # true/false
  
  echo "=== Processing $bam_type BAMs ==="
  
  # Get array contents using indirect reference
  case $array_name in
      "SAMPLES")
          local -a keys=("${!SAMPLES[@]}")
          local -a values=("${SAMPLES[@]}")
          ;;
      "DEDUPED_BAMS")
          local -a keys=("${!DEDUPED_BAMS[@]}")
          local -a values=("${DEDUPED_BAMS[@]}")
          ;;
      "SHIFTED_BAMS")
          local -a keys=("${!SHIFTED_BAMS[@]}")
          local -a values=("${SHIFTED_BAMS[@]}")
          ;;
      *)
          echo "Unknown array: $array_name" >&2
          return 1
          ;;
  esac
  
  echo "Number of BAMs to process: ${#keys[@]}"
  echo "Available samples: ${keys[*]}"
  
  for sample in "${keys[@]}"; do
      echo -e "\nProcessing sample: $sample"
      
      # Get BAM path from appropriate array
      local bam_path
      case $array_name in
          "SAMPLES") bam_path="${SAMPLES[$sample]}" ;;
          "DEDUPED_BAMS") bam_path="${DEDUPED_BAMS[$sample]}" ;;
          "SHIFTED_BAMS") bam_path="${SHIFTED_BAMS[$sample]}" ;;
      esac
      echo "Input BAM: $bam_path"
      
      # Skip input sample for certain analyses
      if [[ "$sample" == "input" && "$bam_type" == "shifted" ]]; then
          echo "Skipping input sample for shifted analysis"
          continue
      fi
      
      # Get fragment size if available
      local ext_size=""
      if [[ -v "FRAGMENTS[$sample]" ]]; then
          ext_size="--extsize ${FRAGMENTS[$sample]}"
          echo "Using fragment size: ${FRAGMENTS[$sample]}"
      else
          echo "No fragment size available for $sample"
      fi
      
      # Process each peak calling mode
      for mode in "${!PEAK_MODES[@]}"; do
          echo -e "\n  Peak calling mode: $mode"
          
          # Without input control
          out_prefix=$(get_peak_name "$sample" "$bam_type" "$mode" "noInput")
          echo "  Output prefix: $(basename "$out_prefix")"
          
          # Skip if output exists
          if [[ -f "${out_prefix}.narrowPeak" || -f "${out_prefix}.broadPeak" ]]; then
              echo "  Skipping existing output: $(basename "$out_prefix")*Peak"
              continue
          fi
          
          echo "  Running MACS2 without input control..."
          echo "  Parameters:"
          echo "    - Treatment: $bam_path"
          echo "    - Mode: ${PEAK_MODES[$mode]}"
          echo "    - Fragment size param: $ext_size"
          
          # todo: Add blacklist processing
          macs2 callpeak \
              -t "$bam_path" \
              -n "$out_prefix" \
              "${CORE_PARAMS[@]}" \
              "${PEAK_MODES[$mode]}" \
              "$ext_size"
          
          # With input control (if allowed and available)
          if [[ "$allow_input" == true && "$sample" == "test" ]]; then
              out_prefix=$(get_peak_name "$sample" "$bam_type" "$mode" "withInput")
              echo -e "\n  Processing with input control"
              echo "  Output prefix: $(basename "$out_prefix")"
              
              if [[ -f "${out_prefix}.narrowPeak" || -f "${out_prefix}.broadPeak" ]]; then
                  echo "  Skipping existing output: $(basename "$out_prefix")*Peak"
                  continue
              fi
              
              echo "  Running MACS2 with input control..."
              echo "  Parameters:"
              echo "    - Treatment: $bam_path"
              echo "    - Control: ${DEDUPED_BAMS[input]}"
              echo "    - Mode: ${PEAK_MODES[$mode]}"
              echo "    - Fragment size param: $ext_size"
              
              macs2 callpeak \
                  -t "$bam_path" \
                  -c "${DEDUPED_BAMS[input]}" \
                  -n "$out_prefix" \
                  "${CORE_PARAMS[@]}" \
                  "${PEAK_MODES[$mode]}" \
                  "$ext_size"
          fi
      done
  done
  echo -e "\nCompleted processing $bam_type BAMs"
}


echo -e "\n=== Read in fragment sizes for peak calls ==="
# Read log file paths into array
mapfile -t FRAG_LOGS < <(
  find "$OUTDIR/predictd/" -type f -name "macs2.log"
) || {
  echo "[ERROR] Failed to find fragment size logs in $OUTDIR/predictd/" >&2
  exit 1
}

# Validate number of log files
[[ ${#FRAG_LOGS[@]} -eq 2 ]] || {
  echo "[ERROR] Expected 2 fragment logs, found ${#FRAG_LOGS[@]}" >&2
  echo "Found logs:" >&2
  printf '%s\n' "${FRAG_LOGS[@]}" >&2
  exit 2
}

# Initialize fragments associative array
declare -A FRAGMENTS=()

# Process each log file
for log in "${FRAG_LOGS[@]}"; do
  # Extract sample type from path
  sample_type=$(basename "$(dirname "$log")" | sed 's/_predictd//')
  #echo "Sample type: $sample_type"
  # Extract fragment size
  frag_size=$(grep -oP 'predicted fragment length is \K\d+' "$log") || {
      echo "[ERROR] Failed to extract fragment size from $log" >&2
      exit 3
  }
  #echo "Read fragment size: $frag_size"
  
  # Validate fragment size
  [[ "$frag_size" =~ ^[0-9]+$ ]] && (( frag_size > 0 )) || {
      echo "[ERROR] Invalid fragment size ($frag_size) in $log" >&2
      exit 4
  }
  
  # Store in associative array
  FRAGMENTS[$sample_type]=$frag_size
done

# Verify expected samples
for required in 'test' 'reference'; do
  #echo "Verifying fragment size in FRAGMENTS array: ${FRAGMENTS[$required]}"
  [[ -n "${FRAGMENTS[$required]}" ]] || {
      echo "[ERROR] Missing fragment size for $required sample" >&2
      exit 5
  }
done

# Debug output
echo "=== Fragment Sizes ==="
for sample in "${!FRAGMENTS[@]}"; do
  echo "$sample: ${FRAGMENTS[$sample]}bp"
done


echo -e "\n=== Starting Peak Calling Analysis ==="
echo "Available BAM sets:"
echo "Original BAMs: ${!SAMPLES[@]}"
echo "Deduped BAMs: ${!DEDUPED_BAMS[@]}"
echo "Shifted BAMs: ${!SHIFTED_BAMS[@]}"

# Process original BAMs
echo -e "\nProcessing original BAMs..."
process_bam_set SAMPLES "raw" true

# Process deduplicated BAMs
echo -e "\nProcessing deduplicated BAMs..."
process_bam_set DEDUPED_BAMS "deduped" true

# Process shifted BAMs
echo -e "\nProcessing shifted BAMs..."
process_bam_set SHIFTED_BAMS "shifted" true

echo -e "\n=== Peak Calling Analysis Complete ==="
