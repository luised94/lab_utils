#STATUS: REMOVE.
#!/bin/bash

module purge 
module add gnu/5.4.0
module add r/4.2.0

date

cd /home/luised94/data/221024Bel_CHIP

Rscript --vanilla -e 'source("../rscripts/package-installation.R")'
Rscript --vanilla -e 'source("../rscripts/project-startup.R")'


date

echo Bash script complete
