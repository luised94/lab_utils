# 
#
#home_directory <- Sys.getenv("HOME")
#library_location <- paste(home_directory, "R/x86_64-pc-linux-gnu-library/4.2", sep = "/")
#Check that library to installs is correct.
#.libPaths()

#TODO: Need to add renv for this module.
library(BiocManager)
bioconductor_packages_for_fastq <- c()
bioconductor_packages_to_install_tracks <- c("QuasR", "GenomicAlignments", "Gviz", "rtracklayer", "ShortRead")
BiocManager::install(bioconductor_packages_to_install,
        lib = library_location)
bioconductor_packages_to_install_peaks <- c()
bioconductor_packages_to_install_motifs <- c()
bioconductor_packages_to_install_visualization <- c()
bioconductor_packages_to_install_statistics <- c()
bioconductor_packages_to_install_reproducibility <- c()
