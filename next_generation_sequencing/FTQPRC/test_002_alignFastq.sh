#!/bin/bash
#SBATCH -N 1 # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1 # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE.
#SBATCH --exclude=c[5-22]
#SBATCH --mem-per-cpu=20G # amount of RAM per node#
#SBATCH --array=1-8%16
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
refgenome_dir="$HOME/data/REFGENS"

#MODULE_LOAD
module purge
module load bowtie2/2.3.5.1

#INITIALIZE_ARRAY
mapfile -t fastq_paths < <(find "${DIR_TO_PROCESS}" -type f -name "processed_*.fastq" )
mapfile -t genomes_paths < <(find "${refgenome_dir}" -type f -name "*_refgenome.fna")

#INPUT_OUTPUT
#fastq_path=${fastq_paths[$SLURM_ARRAY_TASK_ID-1]}
#output_path=$(echo "$fastq_path" | cut -d/ -f7 | xargs -I {} echo "${DIR_TO_PROCESS}processed-fastq/processed_{}")

#LOG
echo "Starting fitering"
echo "FASTQ_FILE: $fastq_path"
echo "OUTPUT_FILE: $output_path"

GENOME_INDEX=$(($SLURM_ARRAY_TASK_ID / ${#fastq_paths[@]}))
FASTQ_INDEX=$(($SLURM_ARRAY_TASK_ID % ${#fastq_paths[@]}))
echo "Processing ${GENOME_INDEX} and ${FASTQ_INDEX}"
#TODO: Crete and echo the output file. Figure out how to output to bam directly,
#COMMAND_TO_EXECUTE 
echo "COMMAND_OUTPUT_START"
#bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]
echo "bowtie2 -x  ${genomes_paths[$GENOME_INDEX]%_refgenome.fna}_index -U ${FASTQ_FILES[$FASTQ_INDEX]} -S ${DIR_TO_PROCESS}alignment/${FASTQ_FILES[$FASTQ_INDEX]##*/}_$(basename ${genomes_paths[$GENOME_INDEX]})_aligned.sam -q --mp 4 --met-stderr"

#LOG
echo "COMMAND_OUTPUT_END"
echo "Filtering completed"
echo "END TIME: $(date "+%Y-%m-%d-%M-%S")"
