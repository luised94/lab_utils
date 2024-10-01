#!/bin/bash
#SBATCH -N 1 # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1 # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE.
#SBATCH --exclude=c[5-22] #Required by MIT
#SBATCH --mem-per-cpu=50G # amount of RAM per node
#SBATCH --cpus-per-task=4
#SBATCH --nice=10000 #Required by MIT
#DESCRIPTION: Generate bigwig files that be used for genome track plots.
#USAGE: Use with slurm wrapper script.
#SETUP
DIR_TO_PROCESS="$1"
timeid=$2

# Define the log directory
LOG_DIR="$HOME/data/$DIR_TO_PROCESS/logs"

# Ensure the log directory exists
mkdir -p "$LOG_DIR"
# Construct the file names
OUT_FILE="${LOG_DIR}/${timeid}_qualityControl_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
ERR_FILE="${LOG_DIR}/${timeid}_qualityControl_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

# Redirect stdout and stderr to the respective files
exec > "$OUT_FILE" 2> "$ERR_FILE"
echo "TASK_START"

#LOG
echo "SLURM_JOB_ID=${SLURM_JOB_ID}, SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}, SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

echo "Started from $(pwd)"
echo "START TIME: $(date "+%Y-%m-%d-%M-%S")"
DIR_TO_PROCESS="$HOME/data/$DIR_TO_PROCESS"
echo "Executing for $DIR_TO_PROCESS"

#MODULE_LOAD
module purge
module load gnu/5.4.0
module load python/2.7.13
module load deeptools 
module load r/4.2.0

#INPUT_OUTPUT
echo "Starting coverage output"
echo "COMMAND_OUTPUT_START"
# R cat command used to return the input and sample paths also returns a NULL. Grep filters it out.
# Returns two paths, first path is the path to sample bam file. Second path is the path to its respective or backup input file.
BASE_DIR=$(basename ${DIR_TO_PROCESS})
readarray -t sample_paths < <(Rscript ~/lab_utils/next_generation_sequencing/004_bamProcessing/determineInputForArrayId.R ${BASE_DIR} ${SLURM_ARRAY_TASK_ID} | grep -v '^NULL$')
# Check if the sample_paths array is non-empty
if [ ${#sample_paths[@]} -eq 0 ]; then
    echo "Error: No output was generated by the R script."
    exit 1
fi

# Check each element in the samples array for emptiness
for sample in "${sample_paths[@]}"; do
    if [ -z "$sample" ]; then
        echo "Error: One of the output lines is empty."
        exit 1
    fi
done
echo "All samples are valid."

SAMPLE=${sample_paths[0]}
INPUT=${sample_paths[1]}
echo "Sample path: ${SAMPLE}"
echo "Input path: ${INPUT}"
# Add awk statement to process sample and input names. 
OUTPUT_FILE=${DIR_TO_PROCESS}bigwig/"${timeid}_$(echo ${SAMPLE%.bam} | awk -F'/' '{print $NF}' | awk -F'_' '{print $1}' )_$(echo ${INPUT%.bam} | awk -F'/' '{print $NF}' | awk -F'_' '{print $1}')_bamcomp.bw"
echo "Output file: ${OUTPUT_FILE}"
HALF_CPU=$((SLURM_CPUS_PER_TASK/2))
echo "Using ${HALF_CPU} CPUS"
bamCompare -b1 ${SAMPLE} -b2 ${INPUT} \
    -o ${OUTPUT_FILE} \
    --binSize 10 \
    --normalizeUsing CPM \
    --scaleFactorsMethod readCount \
    --effectiveGenomeSize 12157105 \
    --ignoreDuplicates \
    --minMappingQuality 30 \
    --operation ratio \
    --ignoreForNormalization chrXII \
    --numberOfProcessors ${HALF_CPU}
#bamCompare -b1 ${SAMPLE} -b2 ${INPUT} \
#    -o ${OUTPUT_FILE} \
#    --binSize 10 \
#    --normalizeUsing CPM \
#    --scaleFactorsMethod SES \
#    --effectiveGenomeSize 12157105 \
#    --ignoreDuplicates \
#    --minMappingQuality 30 \
#    --operation log2 \
#    --pseudocount 1 \
#    --ignoreForNormalization chrXII
#    --numberOfProcessors ${SLURM_CPUS_PER_TASK}/2 
#${SLURM_CPUS_PER_TASK}
#--blackListFileName yeast_blacklist.bed \
#--extendReads 150 \

echo "COMMAND_OUTPUT_END"
echo "Quality control check completed"
echo "END TIME: $(date "+%Y-%m-%d-%M-%S")"
echo "TASK_END"
