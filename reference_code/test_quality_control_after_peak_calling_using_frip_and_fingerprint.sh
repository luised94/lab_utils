#!/bin/bash
#set -euo pipefail

# === Configuration ===
OUTDIR="$HOME/preprocessing_test"
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

# Process each peak file
for peak_file in "${PEAK_FILES[@]}"; do
    # Extract information from peak filename
    peak_filename=$(basename "$peak_file")
    sample_type=${peak_filename%%_*}  # Extract sample type
    #bam_type=$(echo "$filename" | grep -oP '(?<=_)[^_]+(?=_[^_]+_[^_]+_peaks)')  # raw/deduped/shifted
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
        *) echo "Unknown BAM type: $bam_type" >&2; continue ;;
    esac

    echo "Debug message before calculating frip:"
    echo "  BAM file: $bam_file"

    # Verify BAM file exists
    [[ -f "$bam_file" ]] || {
        echo "ERROR: Missing BAM file: $bam_file" >&2
        continue
    }

    total_reads=$(samtools view -c "$bam_file")
    echo "Total reads: $total_reads"
    reads_in_peaks="$FRIP_OUTPUT_DIR/reads_in_peaks.txt"
    bedtools intersect \
        -a "$peak_file" \
        -b "$bam_file" \
        -c > "$reads_in_peaks" || {
            echo "[ERROR] bedtools intersect failed for $(basename "$bam_file")" >&2
            continue
        }
    echo "reads_in_peaks output:"
    cat $reads_in_peaks | head -n 2

    # Calculate FRiP score
    frip=$(awk -v total="$total_reads" '
        {sum+=$NF} 
        END{printf "%.4f\n", sum/total}' "$reads_in_peaks")

    # Define the output file
    output_file="$FRIP_OUTPUT_DIR/${file_output_prefix}_frip_score.txt"

    # Save results if the file does not exist
    [[ ! -f "$output_file" ]] && {
            echo "BAM: $(basename "$bam_file")"
            echo "Peaks: $(basename "$peak_file")"
            echo "Total reads: $total_reads"
            echo "FRiP score: $frip"
            echo "Date: $(date)"
    } > "$output_file"

    echo "FRiP score: $frip"
    echo "Results saved in: $output_file"
done

echo "=== FRiP Analysis Complete ==="

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
