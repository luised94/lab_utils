#!/bin/bash
#SBATCH -N 1 # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1 # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE.
#SBATCH --exclude=c[5-22] #Required by MIT
#SBATCH --mem-per-cpu=50G # amount of RAM per node
#SBATCH --cpus-per-task=4
#SBATCH --nice=10000 #Required by MIT
#DESCRIPTION: Use fastqc program on all fastq files, pre and post processing with fastp.
#USAGE: Use with slurm wrapper specifying the number of tasks. 
#SETUP
DIR_TO_PROCESS="$1"
timeid=$2

# Define the log directory
LOG_DIR="$HOME/data/$DIR_TO_PROCESS/logs"

# Ensure the log directory exists
mkdir -p "$LOG_DIR"
timeid=$(date "+%Y-%m-%d-%M-%S")
# Construct the file names
OUT_FILE="${LOG_DIR}/qualityControl_${SLURM_ARRAY_JOB_ID}.out"
ERR_FILE="${LOG_DIR}/qualityControl_${SLURM_ARRAY_JOB_ID}.err"

# Redirect stdout and stderr to the respective files
exec >> "$OUT_FILE" 2>> "$ERR_FILE"
echo "TASK_START"

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
module load samtools/1.10
module load fastqc/0.11.5

#INITIALIZE_ARRAY
mapfile -t PROCESSEDFQ_PATHS < <(find "${DIR_TO_PROCESS}" -type f -name "processed_*.fastq" )
mapfile -t FASTQ_PATHS < <(find "${DIR_TO_PROCESS}" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \))
#LOG

echo "Starting Quality Control"
#COMMAND_TO_EXECUTE 
echo "COMMAND_OUTPUT_START"
fastqc --outdir=${DIR_TO_PROCESS}qualityControl/ ${FASTQ_PATHS[$SLURM_ARRAY_TASK_ID-1]} 
fastqc --outdir=${DIR_TO_PROCESS}qualityControl/ ${PROCESSEDFQ_PATHS[$SLURM_ARRAY_TASK_ID-1]} 

#LOG
echo "COMMAND_OUTPUT_END"
echo "Aligning completed"
echo "END TIME: $(date "+%Y-%m-%d-%M-%S")"
echo "TASK_END"
