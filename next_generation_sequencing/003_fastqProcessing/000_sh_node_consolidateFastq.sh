DIR_TO_PROCESS="$1"

# Define the log directory
ABSOLUTE_PATH_OF_DIR="$HOME/data/$DIR_TO_PROCESS/"

# INITIALIZE_ARRAY
mapfile -t FASTQ_PATHS < <(find "${ABSOLUTE_PATH_TO_DIR}" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \) | sort )
mapfile -t UNIQUE_IDS < <(printf '%s\n' "${FASTQ_PATHS[@]}" | cut -d- -f2 | uniq )

#${1}_D24-${UNIQUE_ID}_NA_sequence.fastq 
for UNIQUE_ID in ${UNIQUE_IDS[@]}; do 
#	echo $UNIQUE_ID
	printf '%s\n' "${FASTQ_PATHS[@]}" | grep "${UNIQUE_ID}" | xargs -I {} echo "cat {} >> ${DIR_TO_PROCESS}fastq/D24-${UNIQUE_ID}_NA_sequence.fastq"
done
