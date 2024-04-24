
#Running R from command line, use default interpreter 
Rscript -e 'renv::run("/home/luised94/data/R-scripts/1-package-installation.r", project = "/home/luised94/data/221024Bel_CHIP/")'


Rscript -e 'renv::run("/home/luised94/data/R-scripts/2-ngs-project-startup.r", project = "/home/luised94/data/221024Bel_CHIP/")'

Rscript -e 'renv::run("../R-scripts/3-assign-directory-variables.r", project = "/home/luised94/data/221024Bel_CHIP/")'