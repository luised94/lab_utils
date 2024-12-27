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
# Installs R libraries to library_location.
# Inside R source the file. 
# This should be setup in 001_setupR/000_installingR4.2.0.sh
#home_directory <- Sys.getenv("HOME")
#library_location <- paste(home_directory, "R/x86_64-pc-linux-gnu-library/4.2", sep = "/")

install.packeges(c("tidyverse", "R.utils", "ggplot2", "BiocManager", "renv"),
        lib = library_location)
# If there is an error because of the location of the file, you can use .libPaths() to adjust the location that will be used to install. See one line below.
#.libPaths( library_location, .libPaths())

library(BiocManager)
bioconductor_packages_to_install <- c("QuasR", "GenomicAlignments", "Gviz", "rtracklayer", "ShortRead")
BiocManager::install(bioconductor_packages_to_install, lib = library_location)
#Renv has to be installed overall before it can be used to install 
#install.packeges("renv")
renv::init(bioconductor = "3.16")
# May need to install run in command line because of systemfonts, textshaping, and ragg:
# sudo apt-get install libharfbuzz-dev libfribidi-dev libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libboost-all-dev build-essential
renv::install(c("tidyverse", "R.utils", "ggplot2", "BiocManager", "remotes","devtools"))

options(repos = BiocManager::repositories())
#Wasted a bunch of time trying to figure out how to install flowCore dependencies from source. Solved it by installing through github.
devtools::install_github("RGLab/RProtoBufLib")
devtools::install_github("RGLab/cytolib")
devtools::install_github("RGLab/flowCore")
renv::install("ggcyto")

library(flowCore)

renv::snapshot()

#library(BiocManager)
#bioconductor_packages_to_install <- c("QuasR", "GenomicAlignments", "Gviz", "rtracklayer", "ShortRead")
#BiocManager::install(bioconductor_packages_to_install
#		lib = library_location)
