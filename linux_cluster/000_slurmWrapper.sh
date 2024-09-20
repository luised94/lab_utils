#!/bin/bash
#USAGE: 
#~/lab_utils/next_generation_sequencing/linux_cluster/000_sh_node_slurmWrapper.sh 1-N%16 000_scriptToRun.sh 240304Bel
echo -e "Executing from $(pwd) \n"

if [ $# -ne 3 ]; then
    echo -e "Usage: \n $0 <array_number> <script_name> <DIRECTORY_TO_PROCESS_to_process> \n" 
    echo -e 'Array number is 1, an integer (for a specific task) or a range (1-N%16) depending on the number of array tasks to create.(--array= option for SBATCH) \n' 
    echo -e 'script_name is basename of lab_utils script ( eg 003_sh_slurm_alignFastq.sh ) \n' 
    echo -e 'Directory is name of directory without /. Must be in ~/data ( eg 240304Bel ) \n' 
    exit 1
fi

array_range="$1"
SCRIPT_TO_RUN=$(find $HOME/lab_utils -type f -name "$2")
DIRECTORY_TO_PROCESS="$(find -H $HOME/data -maxdepth 1 -type d -name "$3")"

echo -e "Will process ${DIRECTORY_TO_PROCESS}\n"
echo -e "Will run ${SCRIPT_TO_RUN} \n"

# Check if script and DIRECTORY_TO_PROCESS exist
if [ ! -f "$SCRIPT_TO_RUN" -o ! -d "$DIRECTORY_TO_PROCESS" ]; then
  echo -e "Error: Script or DIRECTORY_TO_PROCESS not found!\n"
  echo -e "Echo statement above without a noun shows which one is missing.\n"
  exit 1
fi


timeid=$(date +%Y%m%d%M%S)
#awk '{print substr($0, index($0, last"/")) "/"}' <<< "$TEST_DIR"
#echo $TEST_DIR | rev | cut -d/ -f1 | rev | xargs -I {} echo {}/
JOB_ID=$(sbatch --parsable --array="$array_range" "$SCRIPT_TO_RUN" "${DIRECTORY_TO_PROCESS##*/}/" "${timeid}")
echo "Time is ${timeid}"
echo "JOB is ${JOB_ID}"
echo "View logs using vim ${DIRECTORY_TO_PROCESS}/logs/*_${JOB_ID}.out."
echo "View standard error using vim ${DIRECTORY_TO_PROCESS}/logs/*_${JOB_ID}.err."


#EXTRACT_TO_SCRIPT: cleanupscript
#sbatch --dependency=afterany:$SLURM_JOB_ID cleanup_script.sh
#slurm_files=$(find . -maxdepth 1 -type f -name "slurm*.out")
#SBATCH --job-name=main_job
#SBATCH --output=/dev/null
#TARGET_DIR="${1:-.}"
# Find and remove SLURM output files in the specified directory
#find "$TARGET_DIR" -name 'slurm-*.out' -exec rm {} +


