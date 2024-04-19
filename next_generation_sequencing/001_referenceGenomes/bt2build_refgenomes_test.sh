#!/bin/bash
#SBATCH -N 1 # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1 # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGEÂ
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE.
#SBATCH --exclude=c[5-22]
#SBATCH --mem-per-cpu=20G # amount of RAM per node#
#SBATCH --job-name=index_genome_test
#SBATCH --array=1-4%16
#SBATCH --error=~/data/REFGENS/logs/indexing_%A_%a.err
#SBATCH --output=~/data/REFGENS/logs/indexing_%A_%a.out
##############################################
#USAGE: From ~/data/REFGENS/, run 'sbatch bt2build_refgenomes_test.sh'

echo "START TIME: $(date "+%Y-%m-%d-%M-%S")" 
echo "From $(pwd)"

module purge
module load gnu/5.4.0
module load bowtie2/2.3.5.1

mapfile -t genome_paths < <(find . -type f -name "*_refgenome.fna")

genome_path=${genome_paths[$SLURM_ARRAY_TASK_ID-1]}

echo "Starting indexing for $genome_path"
#Confirmed bowtie2-build works 
echo $genome_path "${genome_path%_refgenome.fna}_index"

echo "Indexing completed"
echo "START TIME: $(date "+%Y-%m-%d-%M-%S")"
