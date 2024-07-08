#!/bin/bash
#This is only required to be run once if there are more than one fastq file per sample 
#This alters the number of - in the string which affects downstream processing.
# Define the directory to process
ABSOLUTE_PATH_OF_DIR="$HOME/data/$1"
OUTPUT_DIR="${ABSOLUTE_PATH_OF_DIR}fastq/"
mkdir -p "$OUTPUT_DIR"

# INITIALIZE_ARRAY
mapfile -t FASTQ_PATHS < <(find "${ABSOLUTE_PATH_OF_DIR}" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \) | sort )

#Could substitute cut with awk -F '-' '{print $2}'
UNIQUE_IDS=($(printf '%s\n' "${FASTQ_PATHS[@]}" | cut -d- -f2 | uniq ))

#TODO Incorporate the read checking inside the for loop, not tested currently
#files_passed_read_check=0
for UNIQUE_ID in ${UNIQUE_IDS[@]}; do 
	OUTPUT_FILE="${OUTPUT_DIR}D24-${UNIQUE_ID}_NA_sequence.fastq"
	echo "Processing ID: ${UNIQUE_ID}, Output: ${OUTPUT_FILE}"
	#FILTERED_PATHS=($(printf '%s\n' "${FASTQ_PATHS[@]}" | grep "${UNIQUE_ID}")) 
	#echo "cat "${FILTERED_PATHS[@]}" >> ${OUTPUT_FILE}"
	echo "Loop starting"
	#total_reads=0
	for FASTQ_PATH in "${FASTQ_PATHS[@]}"; do
		if [[ $FASTQ_PATH =~ $UNIQUE_ID ]]; then 
			cat "$FASTQ_PATH" >> "$OUTPUT_FILE"
			#reads_in_file=$(grep -c "^@" "$FASTQ_PATH")
			#total_reads=$((total_reads + reads_in_file))		
		fi
	done




done 
echo "Files processed: $(echo ${#FASTQ_PATHS[@]})"
echo "Number of Unique IDS: $(echo ${#UNIQUE_IDS[@]})"


#This also works to do a quick check of the concatenated files since I didnt test them but had to move forward.
#uncatenated_files=$(find ~/data/240304Bel/ -type d -name "*D24*" -exec find {} -type f -name "*.fastq" \;)
#grep -c "^@" $(find ~/data/240304Bel/fastq -type f )
#grep -c "^@" ${uncatenated_files[@]}
#
