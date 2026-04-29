#!/bin/bash
# Description: Download datasets from uniprot as proteome or  specific proteins
# The set of proteins will be compared to a dataset of IDRs from shark.
# Use curl, wget and rest api to obtain datasets
# Usage: Run on the command line or as script.

####################
# Configuration
####################

#!/bin/bash

# Script to download Saccharomyces cerevisiae protein data from UniProt
# Author: Data Analysis Script
# Date: $(date +%Y-%m-%d)
# Purpose: Download yeast proteome or specific proteins for IDR string matching analysis

set -e  # Exit on any error

# Configuration
DATA_DIR="data"
ORGANISM_ID="559292"  # Saccharomyces cerevisiae taxon ID
BASE_URL="https://rest.uniprot.org/uniprotkb"

# Color output for better readability
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Create data directory
mkdir -p "$DATA_DIR"

# Usage function
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -p, --proteome          Download complete S. cerevisiae proteome"
    echo "  -s, --specific IDS      Download specific proteins (comma-separated UniProt IDs)"
    echo "  -r, --reviewed-only     Download only reviewed entries (default: all)"
    echo "  -f, --format FORMAT     Output format (fasta, tsv, json) [default: fasta]"
    echo "  -h, --help             Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 -p                                    # Download entire yeast proteome"
    echo "  $0 -s P12345,Q67890                     # Download specific proteins"
    echo "  $0 -p -r                                # Download only reviewed yeast proteins"
    echo "  $0 -s P12345 -f tsv                     # Download specific protein in TSV format"
}

# Default values
DOWNLOAD_PROTEOME=false
SPECIFIC_IDS=""
REVIEWED_ONLY=false
FORMAT="fasta"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--proteome)
            DOWNLOAD_PROTEOME=true
            shift
            ;;
        -s|--specific)
            SPECIFIC_IDS="$2"
            shift 2
            ;;
        -r|--reviewed-only)
            REVIEWED_ONLY=true
            shift
            ;;
        -f|--format)
            FORMAT="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validate format
case $FORMAT in
    fasta|tsv|json|xml)
        ;;
    *)
        echo -e "${RED}Error: Invalid format '$FORMAT'. Supported: fasta, tsv, json, xml${NC}"
        exit 1
        ;;
esac

# Function to download with curl and error handling
download_data() {
    local url=$1
    local output_file=$2
    local description=$3
    
    echo -e "${BLUE}Downloading $description...${NC}"
    echo "URL: $url"
    echo "Output: $output_file"
    
    # Use curl with proper headers and error handling
    if curl -H "Accept: text/plain; format=$FORMAT" \
            -H "User-Agent: BashScript/1.0" \
            --fail \
            --show-error \
            --location \
            --output "$output_file" \
            "$url"; then
        
        echo -e "${GREEN}û Successfully downloaded to $output_file${NC}"
        
        # Show file info
        if [[ -f "$output_file" ]]; then
            file_size=$(du -h "$output_file" | cut -f1)
            echo "  File size: $file_size"
            
            # Count entries based on format
            case $FORMAT in
                fasta)
                    if command -v grep >/dev/null 2>&1; then
                        seq_count=$(grep -c "^>" "$output_file" 2>/dev/null || echo "Unable to count")
                        echo "  Number of sequences: $seq_count"
                    fi
                    ;;
                tsv)
                    if command -v wc >/dev/null 2>&1; then
                        line_count=$(($(wc -l < "$output_file") - 1))  # Subtract header
                        echo "  Number of entries: $line_count"
                    fi
                    ;;
            esac
        fi
        echo
        return 0
    else
        echo -e "${RED}? Failed to download $description${NC}"
        return 1
    fi
}

# Function to download complete proteome
download_proteome() {
    echo -e "${YELLOW}=== Downloading S. cerevisiae Proteome ===${NC}"
    
    # Build query
    local query="organism_id:$ORGANISM_ID"
    if [[ "$REVIEWED_ONLY" == true ]]; then
        query="$query AND reviewed:true"
        local file_suffix="_reviewed"
        local desc_suffix=" (reviewed only)"
    else
        local file_suffix=""
        local desc_suffix=""
    fi
    
    # Construct URL
    local url="${BASE_URL}/search?query=${query}&format=${FORMAT}"
    local output_file="$DATA_DIR/saccharomyces_cerevisiae_proteome${file_suffix}.${FORMAT}"
    
    download_data "$url" "$output_file" "S. cerevisiae proteome${desc_suffix}"
}

# Function to download specific proteins
download_specific() {
    echo -e "${YELLOW}=== Downloading Specific Proteins ===${NC}"
    
    # Convert comma-separated IDs to array
    IFS=',' read -ra ID_ARRAY <<< "$SPECIFIC_IDS"
    
    local success_count=0
    local total_count=${#ID_ARRAY[@]}
    
    for protein_id in "${ID_ARRAY[@]}"; do
        # Trim whitespace
        protein_id=$(echo "$protein_id" | xargs)
        
        local url="${BASE_URL}/${protein_id}?format=${FORMAT}"
        local output_file="$DATA_DIR/${protein_id}.${FORMAT}"
        
        if download_data "$url" "$output_file" "protein $protein_id"; then
            ((success_count++))
        fi
    done
    
    echo -e "${BLUE}Downloaded $success_count out of $total_count proteins${NC}"
}

# Main execution logic
echo -e "${GREEN}=== UniProt S. cerevisiae Data Downloader ===${NC}"
echo "Data directory: $DATA_DIR"
echo "Format: $FORMAT"
echo "Reviewed only: $REVIEWED_ONLY"
echo

# Check what to download
if [[ "$DOWNLOAD_PROTEOME" == true && -n "$SPECIFIC_IDS" ]]; then
    echo -e "${YELLOW}Both proteome and specific IDs specified. Downloading both...${NC}"
    echo
    download_proteome
    download_specific
elif [[ "$DOWNLOAD_PROTEOME" == true ]]; then
    download_proteome
elif [[ -n "$SPECIFIC_IDS" ]]; then
    download_specific
else
    echo -e "${RED}Error: Please specify either -p (proteome) or -s (specific IDs)${NC}"
    usage
    exit 1
fi

echo -e "${GREEN}=== Download Complete ===${NC}"
echo "Files saved in: $DATA_DIR/"
ls -lh "$DATA_DIR"/*."$FORMAT" 2>/dev/null || echo "No files with .$FORMAT extension found"

echo
echo -e "${BLUE}=== Usage Notes ===${NC}"
echo " For IDR string matching, the complete proteome gives you comprehensive coverage"
echo " Specific protein downloads are faster but require knowing target proteins beforehand"
echo " FASTA format is best for sequence analysis and string matching"
echo " TSV format provides tabular data with customizable fields"

# Suggest next steps based on what was downloaded
if [[ "$DOWNLOAD_PROTEOME" == true ]]; then
    echo
    echo -e "${YELLOW}=== Suggested Next Steps ===${NC}"
    echo "1. Use grep or other tools to search for your IDR sequences in the proteome"
    echo "2. Example: grep -B1 'YOURSEQUENCE' $DATA_DIR/saccharomyces_cerevisiae_proteome*.fasta"
    echo "3. Consider using BLAST for fuzzy matching if exact matches aren't found"
fi
