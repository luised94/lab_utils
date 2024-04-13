#!/bin/bash
#SBATCH -N 1 # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1 # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE.
#SBATCH --exclude=c[5-22]
#SBATCH --mem-per-cpu=20G # amount of RAM per node#
#SBATCH --array=1-4%16
#USAGE: From anywhere, run 'sbatch ~/data/lab_utils/next_generation_sequencing/test_001_filterFastq.sh <dir>'

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

# Your script's commands follow...
echo "SLURM_JOB_ID=${SLURM_JOB_ID}, SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}, SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

echo "Started from $(pwd)"
echo "START TIME: $(date "+%Y-%m-%d-%M-%S")"
refgenome_dir="$HOME/data/REFGENS"
echo "Executing from $DIR_TO_PROCESS"

module purge
module load 
module load 

mapfile -t fastq_paths < <(find "${refgenome_dir}" -type f \( -name "*unmapped*" \) -prune -o -name  "*.fastq")

fastq_path=${genome_paths[$SLURM_ARRAY_TASK_ID-1]}

echo "Starting indexing for $genome_path"
#Confirmed bowtie2-build works
#bowtie2-build $fastq_path "${fastq_path%.fastq}_processed.fastq"
fastp -i in1.fastq  -o out1.fastq --qualified_quality_phred 20 --unqualified_percent_limit 50

echo "Indexing completed"
echo "END TIME: $(date "+%Y-%m-%d-%M-%S")"

