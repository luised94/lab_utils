#!/bin/bash
#SBATCH -N 1 # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1 # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGEÂ
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE.
#SBATCH --exclude=c[5-22]
#SBATCH --mem-per-cpu=20G # amount of RAM per node#
#USAGE: From anywhere, run 'sbatch ~/data/lab_utils/next_generation_sequencing/test_bt2build_refgenomes.sh'

# Define the log directory
LOG_DIR="$HOME/data/REFGENS/logs"

# Ensure the log directory exists
mkdir -p "$LOG_DIR"

# Construct the file names
OUT_FILE="${LOG_DIR}/testIndexing_${SLURM_ARRAY_JOB_ID}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
ERR_FILE="${LOG_DIR}/testIndexing_${SLURM_ARRAY_JOB_ID}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

# Redirect stdout and stderr to the respective files
exec >> "$OUT_FILE" 2>> "$ERR_FILE"
echo "TASK_START"

# Your script's commands follow...
echo "Processing with SLURM_JOB_ID=${SLURM_JOB_ID}, SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}, SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

echo "Started from $(pwd)"
echo "START TIME: $(date "+%Y-%m-%d-%M-%S")" 
refgenome_dir="$HOME/data/REFGENS"
echo "Executing from $refgenome_dir"

module purge
module load gnu/5.4.0
module load bowtie2/2.3.5.1

mapfile -t genome_paths < <(find "${refgenome_dir}" -type f -name "*_refgenome.fna")

genome_path=${genome_paths[$SLURM_ARRAY_TASK_ID-1]}

echo "Starting indexing for $genome_path"
#Confirmed bowtie2-build works 
echo $genome_path "${genome_path%_refgenome.fna}_index"

echo "Indexing completed"
echo "END TIME: $(date "+%Y-%m-%d-%M-%S")"
echo -e "TASK_END"
