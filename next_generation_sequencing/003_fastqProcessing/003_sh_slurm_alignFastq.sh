#!/bin/bash
#SBATCH -N 1 # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1 # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE.
#SBATCH --exclude=c[5-22] #Required by MIT
#SBATCH --mem-per-cpu=50G # amount of RAM per node
#SBATCH --cpus-per-task=4
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
OUT_FILE="${LOG_DIR}/${timeid}_aligning_${SLURM_ARRAY_JOB_ID}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
ERR_FILE="${LOG_DIR}/${timeid}_aligning_${SLURM_ARRAY_JOB_ID}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

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
module load bowtie2/2.3.5.1
module load samtools/1.10

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
FASTQ_ID=$(echo "${FASTQ_PATHS[$FASTQ_INDEX]%.fastq}" | awk -F'/' '{print $NF}'  )
GENOME_NAME=$( echo "${GENOME_PATHS[$GENOME_INDEX]}" | awk -F'/' '{print $NF}' | cut -d_ -f1 )
#GENOME_NAME=$( echo "${GENOME_PATHS[$GENOME_INDEX]}" | cut -d_ -f1 | rev | cut -d/ -f1 | rev )
echo "${FASTQ_ID} | $GENOME_NAME"

echo "Starting alignment"
#COMMAND_TO_EXECUTE 
echo "COMMAND_OUTPUT_START"
#bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]
#For testing grab a subset of the data using seqtk:seqtk sample [-2] [-s seed=11] <in.fa> <frac>|<number>
#seqtk sample -s100 your_fastq_file.fastq 0.1 > subsampled_fastq_file.fastq
bowtie2 -x ${GENOME_PATHS[$GENOME_INDEX]%_refgenome.fna}_index -U ${FASTQ_PATHS[$FASTQ_INDEX]} -p $SLURM_CPUS_PER_TASK -q --mp 4 --met-stderr |
samtools view -@ ${SLURM_CPUS_PER_TASK} -b - |
samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${DIR_TO_PROCESS}alignment/${FASTQ_ID}_${GENOME_NAME}.bam -

samtools index ${DIR_TO_PROCESS}alignment/${FASTQ_ID}_${GENOME_NAME}.bam

#LOG
echo "COMMAND_OUTPUT_END"
echo "Aligning completed"
echo "END TIME: $(date "+%Y-%m-%d-%M-%S")"
