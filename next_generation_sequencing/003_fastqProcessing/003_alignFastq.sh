#!/bin/bash
#SBATCH -N 1 # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1 # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE.
#SBATCH --exclude=c[5-22] #Required by MIT
#SBATCH --mem-per-cpu=50G # amount of RAM per node
#SBATCH --cpus-per-task=4
#SBATCH --nice=10000 #Required by MIT
#DESCRIPTION: Align fastq files to multiple genomes using Bowtie2 and Samtools via Slurm Sbatch command.
#USAGE: sbatch --array=1-N 003_alignFastq.sh "240808Bel"
# where N is (number of genomes) * (number of fastq files)
#OR 
#USE VIA 000_slurmWrapper.sh 
#000_slurmWrapper.sh 1-N%16 003_alignFastq.sh "240808Bel"

set -e 
#set -u

#Global variables
#REFGENOME_DIR="$HOME/data/REFGENS"
#EXPERIMENT_DIR=""
#LOG_DIR=""
#LOG_FILE=""

#validate_input() {
#    if [[ $# -ne 1 || "$1" == */ ]]; then
#    echo "Usage: sbatch --array=<array-range> $0 <experiment_name>" >&2
#        echo "Error: Invalid input. Please provide a single argument without a trailing slash. >&2
#        echo "Example: sbatch --array=1-20%16 $0 240808Bel"
#        exit 1
#    fi
#    
#    EXPERIMENT_DIR="$HOME/data/$1"
#    LOG_DIR="$HOME/data/$EXPERIMENT_DIR/logs"
#
#}
#setup_logging() {
#    mkdir -p "$LOG_DIR"
#    local timeid=$(date "+%Y%m%d_%H%M%S")
#
#    LOG_FILE="${LOG_DIR}/aligning_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${timeid}.log"
#
#    # Redirect stdout and stderr to a temporary file
#    exec 3>&1 4>&2
#    exec 1>"$LOG_FILE" 2>&1
#}
#
#log_message() {
#    local timestamp=$(date "+%Y-%m-%d_%H:%M:%S")
#    echo "[$timestamp] $1"
#}
#
#load_modules() {
#    log_message "Loading required modules"
#    module purge
#    module load gnu/5.4.0 bowtie2/2.3.5.1 samtools/1.10
#}
#
#initialize_arrays() {
#    mapfile -t FASTQ_PATHS < <(find "$EXPERIMENT_DIR
#SETUP
DIR_TO_PROCESS="$1"

REFGENOME_DIR="$HOME/data/REFGENS"
LOG_DIR="$HOME/data/${DIR_TO_PROCESS}/logs"
# Define the log directory

# Ensure the log directory exists
mkdir -p "$LOG_DIR"
timeid=$(date "+%Y-%m-%d-%M-%S")
# Construct the file names
OUT_FILE="${LOG_DIR}/aligning_${SLURM_ARRAY_JOB_ID}.out"
ERR_FILE="${LOG_DIR}/aligning_${SLURM_ARRAY_JOB_ID}.err"

# Redirect stdout and stderr to the respective files
exec >> "$OUT_FILE" 2>> "$ERR_FILE"
echo "TASK_START"

#LOG
echo "SLURM_JOB_ID=${SLURM_JOB_ID}, SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}, SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

echo "Started from $(pwd)"
echo "START TIME: $(date "+%Y-%m-%d-%M-%S")"
DIR_TO_PROCESS="$HOME/data/$DIR_TO_PROCESS"
echo "Executing for $DIR_TO_PROCESS"


#INITIALIZE_ARRAY
mapfile -t FASTQ_PATHS < <(find "${DIR_TO_PROCESS}" -type f -name "*.fastq" )
mapfile -t GENOME_PATHS < <(find "${REFGENOME_DIR}" -type f -name "*_refgenome.fna")

#INPUT_OUTPUT
#fastq_path=${FASTQ_PATHS[$SLURM_ARRAY_TASK_ID-1]}
#output_path=$(echo "$fastq_path" | cut -d/ -f7 | xargs -I {} echo "${DIR_TO_PROCESS}processed-fastq/processed_{}")

echo "Number of files to process: $(find "${DIR_TO_PROCESS}" -type f -name "*.fastq" | wc -l )"
echo "Number of files to process: $(find "${REFGENOME_DIR}" -type f -name "*_refgenome.fna" | wc -l )"
#LOG
# Perform indexing by getting quotient and remainder
GENOME_INDEX=$(( ($SLURM_ARRAY_TASK_ID - 1) / ${#FASTQ_PATHS[@]}))
FASTQ_INDEX=$(( ($SLURM_ARRAY_TASK_ID - 1)  % ${#FASTQ_PATHS[@]}))
echo "Processing ${GENOME_INDEX} and ${FASTQ_INDEX}"

#Get names for output
FASTQ_ID=$(echo "${FASTQ_PATHS[$FASTQ_INDEX]%.fastq}" | awk -F'/' '{print $NF}'  )
GENOME_NAME=$( echo "${GENOME_PATHS[$GENOME_INDEX]}" | awk -F'/' '{print $NF}' | cut -d_ -f1 )
#GENOME_NAME=$( echo "${GENOME_PATHS[$GENOME_INDEX]}" | cut -d_ -f1 | rev | cut -d/ -f1 | rev )
echo "FASTQ_ID and GENOME_NAME: "${FASTQ_ID}" | "$GENOME_NAME""

echo "Starting alignment"
#COMMAND_TO_EXECUTE 
echo "COMMAND_OUTPUT_START"
#bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]
#For testing grab a subset of the data using seqtk:seqtk sample [-2] [-s seed=11] <in.fa> <frac>|<number>
#seqtk sample -s100 your_fastq_file.fastq 0.1 > subsampled_fastq_file.fastq

module purge
module load gnu/5.4.0 bowtie2/2.3.5.1 samtools/1.10

# Align FASTQ reads to genome, convert to sorted BAM, and index
# Uses Bowtie2 for alignment, Samtools for conversion, sorting, and indexing
# Outputs: {FASTQ_ID}_{GENOME_NAME}.bam and its index file
bowtie2 -x ${GENOME_PATHS[$GENOME_INDEX]%_refgenome.fna}_index -U ${FASTQ_PATHS[$FASTQ_INDEX]} -p $SLURM_CPUS_PER_TASK -q --mp 4 --met-stderr |
samtools view -@ ${SLURM_CPUS_PER_TASK} -b - |
samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${DIR_TO_PROCESS}alignment/${FASTQ_ID}_${GENOME_NAME}.bam -

samtools index ${DIR_TO_PROCESS}alignment/${FASTQ_ID}_${GENOME_NAME}.bam

#LOG
echo "COMMAND_OUTPUT_END"
echo "Aligning completed"
echo "END TIME: $(date "+%Y-%m-%d-%M-%S")"
echo "TASK_END"
