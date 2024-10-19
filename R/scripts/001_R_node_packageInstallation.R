# May need to install run in command line because of systemfonts, textshaping, and ragg:
# sudo apt-get install libharfbuzz-dev libfribidi-dev libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libboost-all-dev build-essential
renv::install(c("tidyverse", "R.utils", "ggplot2", "BiocManager", "remotes","devtools", "svglite"))
repository <- "github::"
user <- "RGLab/"
packages <- c("RProtoBufLib", "cytolib", "flowCore")
packages_to_install <- paste0(repository, user, packages)
renv::install(c(packages_to_install, "ggcyto"))
renv::install(c("tidyverse", "R.utils", "ggplot2", "BiocManager", "remotes","devtools"))

options(repos = BiocManager::repositories())

library(flowCore)

renv::snapshot()

