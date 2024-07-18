#!/bin/bash
#DESCRIPTION: Consolidate fastq files that are spread into multiple files into one.
#USAGE: $000_consolidateFastq.sh EatonBel
#This is only required to be run once if there are more than one fastq file per sample
#This alters the number of - in the string which affects downstream processing.
#TODO: Rerun to ensure that the script works generally. Need to remember that it uses cut. 
#TODO: Incorporate the read checking inside the for loop, not tested currently
#TODO: Test the logic for concatenating the files. 
# Define the directory to process
ABSOLUTE_PATH_OF_DIR="$HOME/data/$1"
OUTPUT_DIR="${ABSOLUTE_PATH_OF_DIR}fastq/"
mkdir -p "$OUTPUT_DIR"

# INITIALIZE_ARRAY
mapfile -t FASTQ_PATHS < <(find "${ABSOLUTE_PATH_OF_DIR}" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \) | sort )

#Could substitute cut with awk -F '-' '{print $2}'
UNIQUE_IDS=($(printf '%s\n' "${FASTQ_PATHS[@]}" | cut -d- -f2 | uniq ))

# Concatenate FASTQ files for each unique ID, filtering paths containing the ID
# and writing to output files named D24-{UNIQUE_ID}_NA_sequence.fastq
for UNIQUE_ID in ${UNIQUE_IDS[@]}; do
    OUTPUT_FILE="${OUTPUT_DIR}D24-${UNIQUE_ID}_NA_sequence.fastq"
    echo "Processing ID: ${UNIQUE_ID}, Output: ${OUTPUT_FILE}"
    for FASTQ_PATH in "${FASTQ_PATHS[@]}"; do
        if [[ $FASTQ_PATH =~ $UNIQUE_ID ]]; then 
            cat "$FASTQ_PATH" >> $OUTPUT_FILE
            if [ $? -eq 0 ]; then
                rm "$FASTQ_PATH"
                echo "Appended and deleted $FASTQ_PATH to $OUTPUT_FILE"
            else
                echo "Error appending $FASTQ_PATH. Not deleted"
            fi
        fi
    done
done

echo "Files processed: $(echo ${#FASTQ_PATHS[@]})"
echo "Number of Unique IDS: $(echo ${#UNIQUE_IDS[@]})"
echo "Files in directory: $(find "${ABSOLUTE_PATH_OF_DIR}" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \) | sort | wc -l)"
