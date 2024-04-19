#!/bin/bash
#SBATCH -N 1 # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1 # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE.
#SBATCH --exclude=c[5-22] #Required by MIT
#SBATCH --mem-per-cpu=50G # amount of RAM per node
#SBATCH --cpus-per-task=4
#SBATCH --array=1-392%16 # Required %16 to limit number of tasks created to 16 for resource management purposes. Needs to be calculated before hand. TODO: Write wrapper script to dynamically calculate taks number.:w
#SBATCH --nice=10000 #Required by MIT
#USAGE: First, determine this by running the INITIALIZE_ARRAY and multiplying by number of genomes, modify the array number. For test, leave at 1-2 to test array creation. Then, from anywhere, run 'sbatch ~/data/lab_utils/next_generation_sequencing/slurm_002_alignFastq.sh <dir>'
#SETUP
DIR_TO_PROCESS="$1"

# Define the log directory
LOG_DIR="$HOME/data/$DIR_TO_PROCESS/logs"

# Ensure the log directory exists
mkdir -p "$LOG_DIR"
timeid=$(date "+%Y-%m-%d-%M-%S")
# Construct the file names
OUT_FILE="${LOG_DIR}/${timeid}_qualityControl_${SLURM_ARRAY_JOB_ID}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
ERR_FILE="${LOG_DIR}/${timeid}_qualityControl_${SLURM_ARRAY_JOB_ID}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

# Redirect stdout and stderr to the respective files
exec >"$OUT_FILE" 2>"$ERR_FILE"

#LOG
echo "SLURM_JOB_ID=${SLURM_JOB_ID}, SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}, SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

echo "Started from $(pwd)"
echo "START TIME: $(date "+%Y-%m-%d-%M-%S")"
DIR_TO_PROCESS="$HOME/data/$DIR_TO_PROCESS"
echo "Executing for $DIR_TO_PROCESS"
REFGENOME_DIR="$HOME/data/REFGENS"

#MODULE_LOAD
module purge
module load gnu/5.4.0
module load python/2.7.13
module load deeptools 

#INITIALIZE_ARRAY
mapfile -t BAM_PATHS < <(find "${DIR_TO_PROCESS}" -type f -name "*.bam" )

echo "NUMBEROFFILESPROCESSED: $(find "${DIR_TO_PROCESS}" -type f -name "*.bam" | wc -l ) "
#INPUT_OUTPUT
#fastq_path=${FASTQ_PATHS[$SLURM_ARRAY_TASK_ID-1]}
#output_path=$(echo "$fastq_path" | cut -d/ -f7 | xargs -I {} echo "${DIR_TO_PROCESS}processed-fastq/processed_{}")

#LOG

echo "Starting coverage output"

#COMMAND_TO_EXECUTE 
echo "COMMAND_OUTPUT_START"
OUTPUT_FILE=${DIR_TO_PROCESS}bigwig/"$(echo ${BAM_PATHS[$SLURM_ARRAY_TASK_ID]%.bam}.bw | cut -d/ -f7)"

bamCoverage -b ${BAM_PATHS[$SLURM_ARRAY_TASK_ID]} -o ${OUTPUT_FILE} --binSize 10 --normalizeUsing RPKM

#LOG
echo "COMMAND_OUTPUT_END"
echo "Quality control check completed"
echo "END TIME: $(date "+%Y-%m-%d-%M-%S")"
