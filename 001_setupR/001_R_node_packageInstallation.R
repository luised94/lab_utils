
# Installs R libraries to library_location.
# Inside R source the file. 
# This should be setup in 000_installingR4.2.0.sh
#home_directory <- Sys.getenv("HOME")
#library_location <- paste(home_directory, "R/x86_64-pc-linux-gnu-library/4.2", sep = "/")

install.packeges(c("tidyverse", "R.utils", "ggplot2", "BiocManager", "renv"),
		lib = library_location)
# If there is an error because of the location of the file, you can use .libPaths() to adjust the location that will be used to install. See one line below.
#.libPaths( library_location, .libPaths())

library(BiocManager)
bioconductor_packages_to_install <- c("QuasR", "GenomicAlignments", "Gviz", "rtracklayer", "ShortRead")
BiocManager::install(bioconductor_packages_to_install
		lib = library_location)
