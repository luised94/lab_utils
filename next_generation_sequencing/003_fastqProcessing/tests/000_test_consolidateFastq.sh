#STATUS: REMOVE.
# Define the directory to process
ABSOLUTE_PATH_OF_DIR="$HOME/data/$1"
OUTPUT_DIR="${ABSOLUTE_PATH_OF_DIR}fastq/"
mkdir -p "$OUTPUT_DIR"

# INITIALIZE_ARRAY
mapfile -t FASTQ_PATHS < <(find "${ABSOLUTE_PATH_OF_DIR}" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \) | sort )

#Could substitute cut with awk -F '-' '{print $2}'
#TODO may have to replace this by adding basename, ##*/ or awk with $NF 
UNIQUE_IDS=($(printf '%s\n' "${FASTQ_PATHS[@]}" | cut -d- -f2 | uniq ))

#${1}_D24-${UNIQUE_ID}_NA_sequence.fastq 
for UNIQUE_ID in ${UNIQUE_IDS[@]}; do 
	OUTPUT_FILE="${OUTPUT_DIR}D24-${UNIQUE_ID}_NA_sequence.fastq"
	echo "Processing ID: ${UNIQUE_ID}, Output: ${OUTPUT_FILE}"
	FILTERED_PATHS=($(printf '%s\n' "${FASTQ_PATHS[@]}" | grep "${UNIQUE_ID}")) 
	echo "cat "${FILTERED_PATHS[@]}" >> ${OUTPUT_FILE}"
	echo "Loop output"
	for FASTQ_PATH in "${FASTQ_PATHS[@]}"; do
		if [[ $FASTQ_PATH =~ $UNIQUE_ID ]]; then 
			echo "cat "$FASTQ_PATH" >> "$OUTPUT_FILE""
		fi
	done
	
done 
echo "Files processed: $(echo ${#FASTQ_PATHS[@]})"
echo "Number of Unique IDS: $(echo ${#UNIQUE_IDS[@]})"
