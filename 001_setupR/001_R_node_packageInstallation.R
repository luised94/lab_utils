
# Installs R libraries to library_location.
# Inside R source the file. 
# This should be setup in 000_installingR4.2.0.sh
#home_directory <- Sys.getenv("HOME")
#library_location <- paste(home_directory, "R/x86_64-pc-linux-gnu-library/4.2", sep = "/")

#Renv has to be installed overall before it can be used to install 
#install.packeges("renv")
renv::init(bioconductor = "3.16")
# May need to install run in command line because of systemfonts, textshaping, and ragg:
# sudo apt-get install libharfbuzz-dev libfribidi-dev libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libboost-all-dev build-essential
renv::install(c("tidyverse", "R.utils", "ggplot2", "BiocManager"))

options(repos = BiocManager::repositories())

#
#renv::install("BH@1.75.0-0")
renv::install(c("flowCore"))

library(flowCore)

renv::snapshot()

#library(BiocManager)
#bioconductor_packages_to_install <- c("QuasR", "GenomicAlignments", "Gviz", "rtracklayer", "ShortRead")
#BiocManager::install(bioconductor_packages_to_install
#		lib = library_location)
