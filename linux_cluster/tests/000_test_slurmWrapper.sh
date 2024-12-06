#STATUS: REMOVE.
#!/bin/bash
#USAGE: 
#~/lab_utils/next_generation_sequencing/linux_cluster/000_sh_node_slurmWrapper.sh 1-N%16 000_scriptToRun.sh 240304Bel
if [ $# -ne 3 ]; then
	echo "Usage: $0 <array_number> <<script_name> <DIRECTORY_TO_PROCESS_to_process>" 
	echo 'Array number is an integer or range (1-N%16) depending on the number of array tasks to create (--array= option for SBATCH)' 
	echo 'script_name is basename of lab_utils script ( eg 000_sh_node_test_slurmWrapper.sh )' 
	echo 'Directory is DIRECTORY_TO_PROCESS to process with / ( eg 240304Bel/ )' 
	exit 1
fi

array_range="$1"
SCRIPT_TO_RUN=$(find $HOME/lab_utils -type f -name "$2")
DIRECTORY_TO_PROCESS="$(find -H $HOME/data -maxdepth 1 -type d -name "$3")"

echo "Will process $DIRECTORY_TO_PROCESS"
echo "Will run $SCRIPT_TO_RUN"

# Check if script and DIRECTORY_TO_PROCESS exist
if [ ! -f "$SCRIPT_TO_RUN" -o ! -d "$DIRECTORY_TO_PROCESS" ]; then
  echo "Error: Script or DIRECTORY_TO_PROCESS not found!"
  exit 1
fi


#awk '{print substr($0, index($0, last"/")) "/"}' <<< "$TEST_DIR"
#echo $TEST_DIR | rev | cut -d/ -f1 | rev | xargs -I {} echo {}/
JOB_ID=$(sbatch --parsable --array="$array_range" "$SCRIPT_TO_RUN" "${DIRECTORY_TO_PROCESS##*/}/" )

echo "JOB is ${JOB_ID}"
echo "View logs using vim ${DIRECTORY_TO_PROCESS}/logs/*_${JOB_ID}_*_1.out."


#slurm_files=$(find . -maxdepth 1 -type f -name "slurm*.out")



