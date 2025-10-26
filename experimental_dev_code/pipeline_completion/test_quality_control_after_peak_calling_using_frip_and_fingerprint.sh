#!/bin/bash
#set -euo pipefail

# === Configuration ===
OUTDIR="$HOME/data/preprocessing_test"
THREADS=8

# Load required modules
module load samtools bedtools

# === Input Validation ===
# Check required directories and tools
[[ -d "$OUTDIR/peaks" ]] || { echo "ERROR: Missing peaks directory" >&2; exit 1; }
command -v bedtools &>/dev/null || { echo "ERROR: bedtools not found" >&2; exit 1; }
command -v samtools &>/dev/null || { echo "ERROR: samtools not found" >&2; exit 1; }

# === Define BAM Files ===
declare -A SAMPLES=(
    ['test']="$HOME/data/241010Bel/alignment/processed_245018_sequence_to_S288C_sorted.bam"
    ['input']="$HOME/data/241010Bel/alignment/processed_245003_sequence_to_S288C_sorted.bam"
    ['reference']="$HOME/data/100303Bel/alignment/consolidated_034475_sequence_to_S288C_sorted.bam"
)

##############################################
# Function definitions
##############################################
# Function to get paths (reuse from previous script)
get_deduped_path() {
    local sample_type=$1
    echo "$OUTDIR/align/${sample_type}_deduped.bam"
}

get_shifted_path() {
    local sample_type=$1
    echo "$OUTDIR/align/${sample_type}_shifted.bam"
}

##############################################
# Frip and fingerprint workflow
##############################################
# Initialize BAM arrays
declare -i missing_files=0
declare -A DEDUPED_BAMS=()
declare -A SHIFTED_BAMS=()
for sample_type in "${!SAMPLES[@]}"; do
    input="${SAMPLES[$sample_type]}"
    deduped_path=$(get_deduped_path "$sample_type")
    shifted_path=$(get_shifted_path "$sample_type")
    echo "Processing $sample_type..."
    echo "Input: $input"
    echo "deduped: $deduped_path"
    echo "Shifted: $shifted_path"

    if [[ ! -e "$input" ]]; then
        echo "[ERROR] Missing $sample_type file: $input" 1>&2
        missing_files=1
    fi

    if [[ ! -e "$deduped_path" ]]; then
        echo "Deduplicated file not present. Rerun test_bam_preprocessing_with_picard_macs2_deeptools.sh script."
        missing_files=1
    fi

    if [[ "$sample_type" == "test" || "$sample_type" == "reference" ]]; then
        if [[ ! -e "$shifted_path" ]]; then
            echo "Shifted file not present for sample type '$sample_type'. Rerun test_bam_preprocessing_with_picard_macs2_deeptools.sh script."
            missing_files=1
        fi
    fi

    DEDUPED_BAMS[$sample_type]="$deduped_path"
    SHIFTED_BAMS[$sample_type]="$shifted_path"
done

(( missing_files )) && { echo "Aborting: Missing files"; echo "Run test_bam_preprocessing_with_picard_macs2_deeptools.sh script." ; exit 1; }
echo "=== All bam files found ==="

# === Process Each Analysis Type ===
echo "=== Starting FRiP Analysis ==="

# Find all peak files
mapfile -t PEAK_FILES < <(find "$OUTDIR/peaks" -name "*_peaks.narrowPeak" -o -name "*_peaks.broadPeak")

[[ ${#PEAK_FILES[@]} -eq 0 ]] && { 
    echo "ERROR: No peak files found" >&2
    exit 2
}

echo "Found ${#PEAK_FILES[@]} peak files"
# Create output directory
FRIP_OUTPUT_DIR="$OUTDIR/frip"
mkdir -p "$FRIP_OUTPUT_DIR"

declare -A TOTAL_READS_CACHE

for peak_file in "${PEAK_FILES[@]}"; do
    # Extract information from peak filename
    peak_filename=$(basename "$peak_file")
    sample_type=${peak_filename%%_*}             # Extract sample type
    bam_type=$(echo "$peak_filename" | awk -F'_' '{print $2}')  # raw/deduped/shifted
    file_output_prefix=${peak_filename%_peaks*}
    
    echo "Processing ${peak_filename}..."
    echo "  Sample type: $sample_type"
    echo "  Bam type: $bam_type"
    echo "  Filename prefix: $file_output_prefix"
    
    # Select corresponding BAM file
    case $bam_type in
        "raw") bam_file="${SAMPLES[$sample_type]}" ;;
        "deduped") bam_file="${DEDUPED_BAMS[$sample_type]}" ;;
        "shifted") bam_file="${SHIFTED_BAMS[$sample_type]}" ;;
        *)
            echo "Unknown BAM type: $bam_type" >&2
            continue
            ;;
    esac
    
    echo "Debug message before calculating frip:"
    echo "  BAM file: $bam_file"
    
    # Verify BAM file exists
    if [[ ! -f "$bam_file" ]]; then
        echo "ERROR: Missing BAM file: $bam_file" >&2
        continue
    fi

    # Use a combined key to cache total_reads by sample and bam type
    cache_key="${sample_type}_${bam_type}"
    echo "Current cache key: $cache_key"
    if [[ -z "${TOTAL_READS_CACHE[$cache_key]}" ]]; then
        total_reads=$(samtools view -c "$bam_file")
        TOTAL_READS_CACHE[$cache_key]=$total_reads
        echo "Total reads: $total_reads"
    else
        total_reads=${TOTAL_READS_CACHE[$cache_key]}
        echo "Total reads (cached): $total_reads"
    fi

    # Define unique filenames for bedtools and FRiP results
    reads_in_peaks="$FRIP_OUTPUT_DIR/${file_output_prefix}_reads_in_peaks.txt"
    frip_output="$FRIP_OUTPUT_DIR/${file_output_prefix}_frip_score.txt"

    # Run bedtools intersect only if result is not already stored
    if [[ ! -f "$reads_in_peaks" ]]; then
        if ! bedtools intersect -a "$peak_file" -b "$bam_file" -c > "$reads_in_peaks"; then
            echo "[ERROR] bedtools intersect failed for $(basename "$bam_file")" >&2
            continue
        fi
    else
        echo "Using cached bedtools intersection result for $peak_filename"
    fi

    echo "reads_in_peaks output:"
    head -n 2 "$reads_in_peaks"

    # Skip FRiP calculation if output already exists
    if [[ -f "$frip_output" ]]; then
        echo "FRiP score already computed:"
        grep 'FRiP score:' "$frip_output"
        echo "Results saved in: $frip_output"
        continue
    fi

    # Calculate FRiP score
    frip=$(awk -v total="$total_reads" '
        { sum += $NF }
        END { printf "%.4f\n", sum/total }
    ' "$reads_in_peaks")

    # Write results to the output file
    {
        echo "BAM: $(basename "$bam_file")"
        echo "Peaks: $(basename "$peak_file")"
        echo "Total reads: $total_reads"
        echo "FRiP score: $frip"
        echo "Date: $(date)"
    } > "$frip_output"

    echo "FRiP score: $frip"
    echo "Results saved in: $frip_output"
done

echo "=== FRiP Analysis Complete ==="
echo -e "\n=== FRIP Score Summary ==="

# Define thresholds for quality assessment
declare -A FRIP_THRESHOLDS=(
    [excellent]=0.10  # >10%
    [good]=0.05       # >5%
    [moderate]=0.01   # >1%
)

# Collect all FRiP score entries once.
# Each line will contain: <filename> <FRiP_score>
frip_entries=$(find "$FRIP_OUTPUT_DIR" -name "*_frip_score.txt" -exec awk '/FRiP score:/ { print FILENAME, $NF }' {} \;)

if [[ -z "$frip_entries" ]]; then
    echo "No FRiP score files found in $FRIP_OUTPUT_DIR"
    exit 1
fi

echo "Analyzing FRiP scores..."

# --- Threshold-based Listing ---
for threshold_name in "${!FRIP_THRESHOLDS[@]}"; do
    threshold="${FRIP_THRESHOLDS[$threshold_name]}"
    echo -e "\nSamples with $threshold_name FRiP scores (>${threshold}):"
    echo "$frip_entries" | awk -v thr="$threshold" '$2 > thr { print $0 }' | sort -k2,2nr
done

# --- Overall Distribution Statistics ---
echo -e "\nFRiP Score Distribution:"
# Extract just the scores, sort them, and store in an array.
readarray -t scores < <(echo "$frip_entries" | awk '{print $2}' | sort -n)
n=${#scores[@]}

if (( n > 0 )); then
    sum=0
    for score in "${scores[@]}"; do
        sum=$(echo "$sum + $score" | bc -l)
    done
    min=${scores[0]}
    max=${scores[$((n-1))]}
    mean=$(echo "$sum / $n" | bc -l)

    # Calculate median
    if (( n % 2 == 1 )); then
        median=${scores[$((n/2))]}
    else
        mid1=${scores[$((n/2 - 1))]}
        mid2=${scores[$((n/2))]}
        median=$(echo "($mid1 + $mid2) / 2" | bc -l)
    fi

    printf "Min: %.4f\n" "$min"
    printf "Max: %.4f\n" "$max"
    printf "Mean: %.4f\n" "$mean"
    printf "Median: %.4f\n" "$median"
fi

# --- Analysis by BAM (Processing) Type ---
echo -e "\nFRiP Scores by Processing Type:"
for bam_type in "raw" "deduped" "shifted"; do
    echo -e "\n$bam_type:"
    # Filter the entries for the current bam type based on the filename.
    type_entries=$(echo "$frip_entries" | grep -i "$bam_type")
    count=$(echo "$type_entries" | wc -l)
    if (( count > 0 )); then
        sum_type=$(echo "$type_entries" | awk '{sum += $2} END {print sum}')
        avg=$(echo "$sum_type / $count" | bc -l)
        printf "Average: %.4f (%d samples)\n" "$avg" "$count"
    else
        echo "No samples found."
    fi
done

echo "=== Starting Fingerprint Analysis ==="
module purge
module load python/2.7.13 deeptools
command -v plotFingerprint &>/dev/null || {
    echo "ERROR: plotFingerprint command not found. Make sure deepTools is installed." >&2
    exit 1
}
echo "Verified plotFingerprint function."

QC_DIR="$OUTDIR/fingerprint"
mkdir -p "$QC_DIR"
echo "Created $QC_DIR directory"

declare -a samples_to_process=("test" "reference")

# Loop through the BAM processing types to create the plotFingerprint output.
# Do not output png and pdf since that requires certain libraries that may not be available on the linux cluster.
for bam_type in "original" "deduped" "shifted"; do
    echo -e "\nProcessing ${bam_type} BAMs..."

    # Declare the BAM file associative array based on type
    case "$bam_type" in
        "original") 
            declare -A bam_files=( 
                ["test"]="${SAMPLES[test]}" 
                ["reference"]="${SAMPLES[reference]}" 
                ["input"]="${SAMPLES[input]}" 
            )
            ;;
        "deduped") 
            declare -A bam_files=( 
                ["test"]="${DEDUPED_BAMS[test]}" 
                ["reference"]="${DEDUPED_BAMS[reference]}" 
                ["input"]="${DEDUPED_BAMS[input]}" 
            )
            ;;
        "shifted") 
            declare -A bam_files=( 
                ["test"]="${SHIFTED_BAMS[test]}" 
                ["reference"]="${SHIFTED_BAMS[reference]}" 
                # For shifted, we use the deduped input since input is not shifted
                ["input"]="${DEDUPED_BAMS[input]}" 
            )
            ;;
        *)
            echo "Unknown BAM type: $bam_type" >&2
            continue
            ;;
    esac

    # Process each sample for this BAM type
    for sample in "${samples_to_process[@]}"; do
        echo "Processing ${sample} ${bam_type}..."

        # Set common output filename base
        if [[ "$sample" == "test" ]]; then
            base_name="${QC_DIR}/fingerprint_${sample}_vs_input_${bam_type}"
            label1="${sample^}_${bam_type^}"
            label2="Input_${bam_type^}"
            bam_args=( "${bam_files[$sample]}" "${bam_files[input]}" )
        else
            base_name="${QC_DIR}/fingerprint_${sample}_${bam_type}"
            label1="${sample^}_${bam_type^}"
            bam_args=( "${bam_files[$sample]}" )
        fi

        # Check if the output already exists (for PNG output as an example)
        if [[ -f "${base_name}.png" ]]; then
            echo "Fingerprint plot ${base_name}.png already exists, skipping..."
        else
            echo "Generating fingerprint plot for ${sample} ${bam_type}..."
            plotFingerprint \
                -b "${bam_args[@]}" \
                --labels "$label1" ${label2:+ "$label2"} \
                --outRawCounts "${base_name}.tab" \
                --numberOfProcessors "$THREADS"
        fi
    done
done

# === Cross-comparison Fingerprints ===
# Compare processing methods for each sample (test and reference)
for sample in "${samples_to_process[@]}"; do

    echo -e "\nGenerating ${sample} processing comparison..."
    base_name="${QC_DIR}/fingerprint_${sample}_processing_comparison"
    
    if [[ -f "${base_name}.png" ]]; then
        echo "Fingerprint comparison ${base_name}.png already exists, skipping..."
    else
        plotFingerprint \
            -b "${SAMPLES[$sample]}" "${DEDUPED_BAMS[$sample]}" "${SHIFTED_BAMS[$sample]}" \
            --labels "${sample^}_Original" "${sample^}_Deduped" "${sample^}_Shifted" \
            --outRawCounts "${base_name}.tab" \
            --numberOfProcessors "$THREADS"
    fi
done

echo -e "\n=== Fingerprint Analysis Complete ==="
echo "Results saved in $QC_DIR"
