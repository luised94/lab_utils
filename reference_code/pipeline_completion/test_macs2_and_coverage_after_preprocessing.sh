#!/usr/bin/env bash
#todo: Add blacklist processing
#todo: Add quality control checks
#todo: Create the R script to visualize the results of macs2 testing.

#set -euo pipefail

# === Cluster Environment Setup ===
#if [[ "$(hostname)" != "luria" ]]; then
#    echo "Error: This script must be run on luria cluster" 1>&2
#    exit 1
#fi

##############################################
# Key Adjustable Parameters
##############################################

# Cluster config
THREADS=8

# File paths
declare -A SAMPLES=(
    ['test']="$HOME/data/241010Bel/alignment/processed_245018_sequence_to_S288C_sorted.bam"
    ['input']="$HOME/data/241010Bel/alignment/processed_245003_sequence_to_S288C_sorted.bam" 
    ['reference']="$HOME/data/100303Bel/alignment/processed_034475_sequence_to_S288C_sorted.bam"
)
declare -A CHROM_SIZES=()

# Reference data
GENOME_SIZE=12000000  # 1.2e7 in integer form
GENOME_FASTA="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C/SaccharomycescerevisiaeS288C_refgenome.fna"
if [[ ! -f $GENOME_FASTA ]]; then
    echo "$GENOME_FASTA does not exist."
    exit 1
fi

# Define normalization methods
declare -a NORMALIZATION=("" "RPKM" "CPM")  # Raw, RPKM, and CPM

# Output config
OUTDIR="$HOME/preprocessing_test"
OUTPUT_PREFIX="test"
SUB_DIRS=("align" "predictd" "peaks" "coverage")

# Processing parameters
PVALUE=1e-6
BIN_SIZE=25
SMOOTH_LEN=75
MIN_FRAGMENT=20
MAX_FRAGMENT=300

##############################################
# Helper functions
##############################################
# File path generators
get_deduped_path() {
    local sample_type=$1
    echo "$OUTDIR/align/${sample_type}_deduped.bam"
}

get_shifted_path() {
    local sample_type=$1
    echo "$OUTDIR/align/${sample_type}_shifted.bam"
}

# Function to generate BigWig name
get_bigwig_name() {
    local sample=$1
    # 'raw', 'deduped', or 'shifted'
    local bam_type=$2
    local norm_method=$3
    local suffix="raw"
    [[ -n "$norm_method" ]] && suffix="${norm_method,,}"
    echo "$OUTDIR/coverage/${sample}_${bam_type}_${suffix}.bw"
}

get_chrom_sizes() {
    local ref_fasta="$GENOME_FASTA"
    local chrom_sizes_file="$OUTDIR/chrom.sizes"

    echo "=== Chromosome Size Generation ==="

    # Generate .fai if missing
    if [[ ! -f "${ref_fasta}.fai" ]]; then
        echo "Creating FASTA index for reference genome..."
        samtools faidx "$ref_fasta" || {
            echo "[ERROR] Failed to index reference FASTA" >&2
            exit 10
        }
    fi

    # Generate chrom.sizes if missing
    if [[ ! -f "$chrom_sizes_file" ]]; then
        echo -e "\nBuilding chromosome size file..."
        mkdir -p "$OUTDIR"

        echo "Input FASTA: $(basename "$ref_fasta")"
        echo "Output sizes: $chrom_sizes_file"

        awk '{print $1 "\t" $2}' "${ref_fasta}.fai" > "$chrom_sizes_file" || {
            echo "[ERROR] Failed to create chrom.sizes" >&2
            exit 11
        }

        # Validate non-empty output
        [[ -s "$chrom_sizes_file" ]] || {
            echo "[ERROR] chrom.sizes file is empty" >&2
            exit 12
        }
    fi

    # Load into associative array with debug output
    echo -e "\nLoading chromosome sizes:"
    while IFS=$'\t' read -r chrom size; do
        [[ -z "$chrom" ]] && continue  # Skip empty lines

        # Trim additional fields after first tab
        chrom="${chrom%%[[:space:]]*}"

        echo " - ${chrom}: ${size}bp"
        CHROM_SIZES["$chrom"]="$size"
    done < "$chrom_sizes_file"

    echo -e "\nLoaded ${#CHROM_SIZES[@]} chromosomes"
    echo "======================================="
}

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
            
            macs2 callpeak \
                -t "$bam_path" \
                -n "$out_prefix" \
                ${CORE_PARAMS[@]} \
                ${PEAK_MODES[$mode]} \
                $ext_size
            
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
                    ${CORE_PARAMS[@]} \
                    ${PEAK_MODES[$mode]} \
                    $ext_size
            fi
        done
    done
    echo -e "\nCompleted processing $bam_type BAMs"
}

##############################################
# Initialization & Validation
##############################################

# Create directory structure
mkdir -p "${SUB_DIRS[@]/#/$OUTDIR/}"

# Display parameter summary
echo "=== Pipeline Configuration ==="
echo "Genome Size: $GENOME_SIZE"
echo "Output Directory: $OUTDIR"
echo "Threads: $THREADS"
echo "Samples:"
for stype in "${!SAMPLES[@]}"; do
    echo " - $stype: ${SAMPLES[$stype]}"
done
echo "=============================="

# === File Verification ===
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

echo "All sample files verified successfully."

module load python/2.7.13 deeptools/3.0.1

# Process original BAMs
echo "=== Generating coverage for original BAMs ==="
for sample_type in "${!SAMPLES[@]}"; do
    input="${SAMPLES[$sample_type]}"
    echo "Processing: $input"
    for norm_method in "${NORMALIZATION[@]}"; do
        output=$(get_bigwig_name "$sample_type" "raw" "$norm_method")
        echo "File will be output to: $output"
        [[ -f "$output" ]] && {
            echo "Skipping existing: $output"
            continue
        }

        echo "Processing $sample_type (raw) with ${norm_method:-no} normalization"
        norm_flag=""
        [[ -n "$norm_method" ]] && norm_flag="--normalizeUsing $norm_method"
        echo "Using norm flag: $norm_flag"
        bamCoverage \
            -b "$input" \
            -o "$output" \
            $norm_flag \
            --binSize 25 \
            --effectiveGenomeSize "$GENOME_SIZE" \
            --smoothLength 75 \
            --ignoreDuplicates \
            --numberOfProcessors "$THREADS"
    done
done

# Process deduped BAMs
echo -e "\n=== Generating coverage for deduplicated BAMs ==="
for sample_type in "${!DEDUPED_BAMS[@]}"; do
    input="${DEDUPED_BAMS[$sample_type]}"
    echo "Processing: $input"
    for norm_method in "${NORMALIZATION[@]}"; do
        output=$(get_bigwig_name "$sample_type" "deduped" "$norm_method")
        echo "File will be output to: $output"
        [[ -f "$output" ]] && {
            echo "Skipping existing: $output"
            continue
        }

        echo "Processing $sample_type (deduped) with ${norm_method:-no} normalization"
        norm_flag=""
        [[ -n "$norm_method" ]] && norm_flag="--normalizeUsing $norm_method"
        echo "Using norm flag: $norm_flag"
        bamCoverage \
            -b "$input" \
            -o "$output" \
            $norm_flag \
            --binSize 25 \
            --effectiveGenomeSize "$GENOME_SIZE" \
            --smoothLength 75 \
            --ignoreDuplicates \
            --numberOfProcessors "$THREADS"
    done
done

# Process shifted BAMs (skip input)
echo -e "\n=== Generating coverage for shifted BAMs ==="
for sample_type in 'test' 'reference'; do
    input="${SHIFTED_BAMS[$sample_type]}"
    echo "Processing: $input"
    [[ -z "$input" ]] && continue  # Skip if no shifted BAM exists
    
    for norm_method in "${NORMALIZATION[@]}"; do
        output=$(get_bigwig_name "$sample_type" "shifted" "$norm_method")
        echo "File will be output to: $output"
        [[ -f "$output" ]] && {
            echo "Skipping existing: $output"
            continue
        }

        echo "Processing $sample_type (shifted) with ${norm_method:-no} normalization"
        norm_flag=""
        [[ -n "$norm_method" ]] && norm_flag="--normalizeUsing $norm_method"
        echo "Using norm flag: $norm_flag"
        bamCoverage \
            -b "$input" \
            -o "$output" \
            $norm_flag \
            --binSize 25 \
            --effectiveGenomeSize "$GENOME_SIZE" \
            --smoothLength 75 \
            --ignoreDuplicates \
            --numberOfProcessors "$THREADS"
    done
done

# Summary of generated files
#echo -e "\n=== Coverage Generation Summary ==="
#find "$OUTDIR/coverage" -name "*.bw" -exec basename {} \; | sort

# Create array of bigwig files
declare -a BIGWIG_FILES
mapfile -t BIGWIG_FILES < <(find "$OUTDIR/coverage" -name "*.bw" | sort)

# Validate array
if [[ ${#BIGWIG_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No bigwig files found in $OUTDIR/coverage" >&2
    exit 3
fi

# Display found files
#echo "Found ${#BIGWIG_FILES[@]} bigwig files:"
#for bw in "${BIGWIG_FILES[@]}"; do
#    echo "  $bw"
#done

# Purge existing modules
echo "Purging modules..."
module purge

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
    echo "conda install -c bioconda macs2=2.2.7.1" >&2
    exit 2
fi

echo "MACS2 environment ready in $(which python)"

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
