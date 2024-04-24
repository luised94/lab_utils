#!/bin/bash
#SBATCH -N 1 # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1 # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE.
#SBATCH --exclude=c[5-22]
#SBATCH --mem-per-cpu=20G # amount of RAM per node#
#SBATCH --array=1-49%16
#USAGE: First, determine this by running the INITIALIZE_ARRAY, modify the array number. Then, from anywhere, run 'sbatch ~/lab_utils/next_generation_sequencing/test_001_filterFastq.sh <dir>'
#SETUP
DIR_TO_PROCESS="$1"

# Define the log directory
LOG_DIR="$HOME/data/$DIR_TO_PROCESS/logs"

# Ensure the log directory exists
mkdir -p "$LOG_DIR"
timeid=$(date "+%Y-%m-%d-%M-%S")
# Construct the file names
OUT_FILE="${LOG_DIR}/${timeid}_filtering_${SLURM_ARRAY_JOB_ID}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
ERR_FILE="${LOG_DIR}/${timeid}_filtering_${SLURM_ARRAY_JOB_ID}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

# Redirect stdout and stderr to the respective files
exec >"$OUT_FILE" 2>"$ERR_FILE"

#LOG
echo "SLURM_JOB_ID=${SLURM_JOB_ID}, SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}, SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

echo "Started from $(pwd)"
echo "START TIME: $(date "+%Y-%m-%d-%M-%S")"
DIR_TO_PROCESS="$HOME/data/$DIR_TO_PROCESS"
echo "Executing for $DIR_TO_PROCESS"

#MODULE_LOAD
module purge
module load fastp/0.20.0

#INITIALIZE_ARRAY
mapfile -t fastq_paths < <(find "${DIR_TO_PROCESS}" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \))

#INPUT_OUTPUT
fastq_path=${fastq_paths[$SLURM_ARRAY_TASK_ID-1]}
output_path=$(echo "$fastq_path" | cut -d/ -f7 | xargs -I {} echo "${DIR_TO_PROCESS}processed-fastq/processed_{}")
json_path=$(echo "$fastq_path" | cut -d/ -f7 | xargs -I {} echo "${DIR_TO_PROCESS}logs/processed_{}")
#LOG
echo "Starting fitering"
echo "FASTQ_FILE: $fastq_path"
echo "OUTPUT_FILE: $output_path"

#Confirmed fastp works from node using example file
#COMMAND_TO_EXECUTE 
echo "COMMAND_OUTPUT_START"
fastp -i "${fastq_path}"  -o "${output_path}" --json "${json_path%.fastq}.json" --html /dev/null --cut_window_size 4 --cut_mean_quality 20 --n_base_limit 5 --average_qual 20 --length_required 50 --qualified_quality_phred 20 --unqualified_percent_limit 50

#LOG
echo "COMMAND_OUTPUT_END"
echo "Filtering completed"
echo "END TIME: $(date "+%Y-%m-%d-%M-%S")"
