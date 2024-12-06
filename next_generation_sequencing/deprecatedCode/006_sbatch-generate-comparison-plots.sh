#STATUS:
#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1                      # Number of tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --mail-type=END           # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE 
#SBATCH --mail-user=luised94@mit.edu  # Email to which notifications will be sent. Equivalent to the -M option in SGE. You must replace [] with your email address.
#SBATCH --exclude=c[5-22]
#SBATCH --mem-per-cpu=90G         # amount of RAM per node
#############################################

module purge 
module add gnu/5.4.0
module add r/4.2.0

cd /home/luised94/data/221024Bel_CHIP

export R_LIBS=/home/luised94/R/x86_64-pc-linux-gnu-library/4.2

Rscript --vanilla /home/luised94/data/rscripts/generate-comparison-coverage-plots.R $1 > /home/luised94/data/221024Bel_CHIP/logs/slurm-${SLURM_JOBID}.Rout 2>&1

cat /home/luised94/data/221024Bel_CHIP/logs/slurm-${SLURM_JOBID}.Rout >> /home/luised94/data/221024Bel_CHIP/logs/slurm-${SLURM_JOBID}.out

rm /home/luised94/data/221024Bel_CHIP/logs/slurm-${SLURM_JOBID}.Rout 

cd /home/luised94/data/rscripts

find . -maxdepth 1 -name 'slurm*' -delete
