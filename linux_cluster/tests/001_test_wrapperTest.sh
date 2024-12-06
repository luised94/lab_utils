#STATUS:
#!/bin/bash
#SBATCH -N 1 # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1 # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE.
#SBATCH --exclude=c[5-22]
#SBATCH --mem-per-cpu=20G # amount of RAM per node#

#SETUP
DIR_TO_PROCESS="$1"

# Define the log directory
LOG_DIR="$HOME/data/$DIR_TO_PROCESS/logs"

# Ensure the log directory exists
mkdir -p "$LOG_DIR"

timeid=$(date "+%Y-%m-%d-%M-%S")
# Construct the file names
OUT_FILE="${LOG_DIR}/${timeid}_aligning_${SLURM_ARRAY_JOB_ID}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
ERR_FILE="${LOG_DIR}/${timeid}_aligning_${SLURM_ARRAY_JOB_ID}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

# Redirect stdout and stderr to the respective files
exec >"$OUT_FILE" 2>"$ERR_FILE"

echo "Started from $(pwd)"
echo "START TIME: $(date "+%Y-%m-%d-%M-%S")"
DIR_TO_PROCESS="$HOME/data/$DIR_TO_PROCESS"
echo "Executing for $DIR_TO_PROCESS"
REFGENOME_DIR="$HOME/data/REFGENS"

#LOG
echo "SLURM_JOB_ID=${SLURM_JOB_ID}, SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}, SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

#INITIALIZE_ARRAY
mapfile -t FASTQ_PATHS < <(find "${DIR_TO_PROCESS}" -type f -name "processed_*.fastq" )
mapfile -t GENOME_PATHS < <(find "${REFGENOME_DIR}" -type f -name "*_refgenome.fna")

#INPUT_OUTPUT
#fastq_path=${FASTQ_PATHS[$SLURM_ARRAY_TASK_ID-1]}
#output_path=$(echo "$fastq_path" | cut -d/ -f7 | xargs -I {} echo "${DIR_TO_PROCESS}processed-fastq/processed_{}")

#LOG
# Perform indexing by getting quotient and remainder
GENOME_INDEX=$(( ($SLURM_ARRAY_TASK_ID - 1) / ${#FASTQ_PATHS[@]}))
FASTQ_INDEX=$(( ($SLURM_ARRAY_TASK_ID - 1)  % ${#FASTQ_PATHS[@]}))
echo "Processing ${GENOME_INDEX} and ${FASTQ_INDEX}"

#Get names for output
FASTQ_ID=$(echo "${FASTQ_PATHS[$FASTQ_INDEX]}" | cut -d_ -f2 )
GENOME_NAME=$( echo "${GENOME_PATHS[$GENOME_INDEX]}" | cut -d_ -f1 | rev | cut -d/ -f1 | rev )
echo "$FASTQ_ID | $GENOME_NAME"

echo "${FASTQ_PATHS[$FASTQ_INDEX]}"
echo "${GENOME_PATHS[$GENOME_INDEX]}"


