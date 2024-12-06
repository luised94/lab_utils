#STATUS:
#!/bin/bash
if [ $# -ne 1 ]; then 
	echo -e "No directory provided. Provide the basename of a directory with no / in the ~/data directory"
fi

DIRECTORY_TO_PROCESS="$(find -H $HOME/data -maxdepth 1 -type d -name "$1")"
LOG_DIR="${DIRECTORY_TO_PROCESS}/logs/"

if [ ! -d "$DIRECTORY_TO_PROCESS" ]; then
  echo "Error: DIRECTORY_TO_PROCESS not found!"
  exit 1
fi

echo "Will process $DIRECTORY_TO_PROCESS"

EXTENSIONS=(".out" ".err")

for EXT in "${EXTENSIONS[@]}"; do
	
	# INITIALIZE_ARRAY
	mapfile -t LOG_PATHS < <(find "${LOG_DIR}" -type f -name "*${EXT}" ! -name "*consolidated*" )
	
	#Could substitute cut with awk -F '-' '{print $2}'
	if [ ${#LOG_PATHS[@]} -eq 0 ]; then
		echo "No Log files to consolidate. Exiting\n"
		continue
	fi

	#Alternative to awk: echo "2024-04-15-19-25_filtering_8980861_8980862_1.out" | cut -d'_' -f1-3 --output-delimiter='_'
	UNIQUE_JOBIDS=($(printf '%s\n' "${LOG_PATHS[@]}" | awk -F'/' '{print $NF}' | awk -F'_' '{print $3}' | uniq ))
	echo "Total number of unique JOB_IDS is: ${#UNIQUE_JOBIDS[@]}"
	
	for JOB_ID in ${UNIQUE_JOBIDS[@]}; do 
		echo "Processing ${JOB_ID}"
		# Find all .out files for the current job ID
		JOB_LOG_PATHS=($(printf '%s\n' "${LOG_PATHS[@]}" | grep "${JOB_ID}" ))
		echo "Total number of LOG Files with ${JOB_ID} is ${#JOB_LOG_PATHS[@]}"
		echo "First LOG PATH is ${JOB_LOG_PATHS[0]}"
		# Extract the oldest date ID for the current job ID
		OLDEST_DATE_ID=$(printf '%s\n' "${JOB_LOG_PATHS[@]}" | awk -F'/' '{print $NF}' | awk -F'_' '{print $1}' | sort -r | head -n 1)
	
		# Construct the output file name
		OUTPUT_FILE="${LOG_DIR}${OLDEST_DATE_ID}_${JOB_ID}_consolidated$EXT"
		echo -e "Output file is ${OUTPUT_FILE}\n"
		
		# Consolidate the log files for the current job ID
		for LOG_FILE in "${JOB_LOG_PATHS[@]}"; do
#			echo "Outputting ${LOG_FILE}."
	  		echo -e "\n\n----- $(basename "${LOG_FILE}") -----\n" >> "${OUTPUT_FILE}"
	  		cat "${LOG_FILE}" >> "${OUTPUT_FILE}"
		done
		rm "${JOB_LOG_PATHS[@]}"
	done
done
echo "Number of files in ${LOG_DIR} is: $(ls ${LOG_DIR} | wc -l)"
