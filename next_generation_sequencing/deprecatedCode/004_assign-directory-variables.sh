#STATUS:
#!/bin/bash
echo Working directory is $(pwd)

module purge
module add gnu/5.4.0
module add r/4.2.0

date

cd /home/luised94/data/221024Bel_CHIP

echo Now, working directory is $(pwd)

Rscript -e 'source("/home/luised94/data/rscripts/3-assign-directory-variables.r")'

date

echo Bash script complete
