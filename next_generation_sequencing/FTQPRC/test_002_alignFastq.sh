#!/bin/bash
#SBATCH -N 1 # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1 # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE.
#SBATCH --exclude=c[5-22]
#SBATCH --mem-per-cpu=50G # amount of RAM per node#
#SBATCH --cpus-per-task=4
#SBATCH --array=1-10%16
#SBATCH --nice=10000
#USAGE: First, determine this by running the INITIALIZE_ARRAY and multiplying by number of genomes, modify the array number. For test, leave at 1-2 to test array creation. Then, from anywhere, run 'sbatch ~/data/lab_utils/next_generation_sequencing/test_002_alignFastq.sh <dir>'
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
module load bowtie2/2.3.5.1
module load samtools/1.10

#INITIALIZE_ARRAY
mapfile -t fastq_paths < <(find "${DIR_TO_PROCESS}" -type f -name "processed_*.fastq" )
mapfile -t genomes_paths < <(find "${REFGENOME_DIR}" -type f -name "*_refgenome.fna")

#INPUT_OUTPUT
#fastq_path=${fastq_paths[$SLURM_ARRAY_TASK_ID-1]}
#output_path=$(echo "$fastq_path" | cut -d/ -f7 | xargs -I {} echo "${DIR_TO_PROCESS}processed-fastq/processed_{}")

#LOG
# Perform indexing by getting quotient and remainder
GENOME_INDEX=$(( ($SLURM_ARRAY_TASK_ID - 1) / ${#fastq_paths[@]}))
FASTQ_INDEX=$(( ($SLURM_ARRAY_TASK_ID - 1)  % ${#fastq_paths[@]}))
echo "Processing ${GENOME_INDEX} and ${FASTQ_INDEX}"

#Get names for output
fastq_ID=$(echo "${fastq_paths[$fastq_index]}" | cut -d_ -f3 )
genome_name=$( echo "${genomes_paths[$genome_index]}" | cut -d_ -f1 | rev | cut -d/ -f1 | rev )
echo "$(basename $fastq_ID) | $(basename $genome_name)"

echo "Starting fitering"
#COMMAND_TO_EXECUTE 
echo "COMMAND_OUTPUT_START"
#bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]
echo "bowtie2 -x ${genomes_paths[$GENOME_INDEX]%_refgenome.fna}_index -U ${fastq_paths[$FASTQ_INDEX]} -p $SLURM_CPUS_PER_TASK -q --mp 4 --met-stderr - |"
echo "samtools view -@ ${SLURM_CPUS_PER_TASK} -bS |"
echo "samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${DIR_TO_PROCESS}alignment/${fastq_ID}_${genome_name}.bam "

echo "samtools index ${DIR_TO_PROCESS}alignment/${fastq_ID}_${genome_name}.bam"

#LOG
echo "COMMAND_OUTPUT_END"
echo "Aligning completed"
echo "END TIME: $(date "+%Y-%m-%d-%M-%S")"
