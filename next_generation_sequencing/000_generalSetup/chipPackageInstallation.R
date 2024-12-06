#STATUS:
# Initialize renv for reproducibility
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::init()

# Load BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(BiocManager)

# Define package lists for different analysis tasks
bioconductor_packages_for_fastq <- c("ShortRead", "Rsubread", "Biostrings", "dada2")
bioconductor_packages_to_install_tracks <- c("QuasR", "GenomicAlignments", "Gviz", "rtracklayer")
bioconductor_packages_to_install_peaks <- c("ChIPseeker", "ChIPpeakAnno", "DiffBind", "normR", "mosaics", "csaw")
bioconductor_packages_to_install_motifs <- c("motifStack", "TFBSTools", "JASPAR2020", "universalmotif", "memes")
bioconductor_packages_to_install_visualization <- c("ggbio", "ComplexHeatmap", "EnhancedVolcano", "ggplot2", "ggcoverage", "gggenome")
bioconductor_packages_to_install_statistics <- c("DESeq2", "edgeR", "limma")

# Additional CRAN packages
cran_packages <- c("ggplot2", "rmarkdown", "knitr", "tidyverse", "furrr")

# Combine all packages into a single vector
all_packages_to_install <- c(
  bioconductor_packages_for_fastq,
  bioconductor_packages_to_install_tracks,
  bioconductor_packages_to_install_peaks,
  bioconductor_packages_to_install_motifs,
  bioconductor_packages_to_install_visualization,
  bioconductor_packages_to_install_statistics,
  cran_packages
)

# Remove duplicates
all_packages_to_install <- unique(all_packages_to_install)

# Install all packages
BiocManager::install(all_packages_to_install)

# Load all installed packages
invisible(lapply(all_packages_to_install, library, character.only = TRUE))

# Print confirmation message
cat("All packages have been installed and loaded successfully.\n")

# Print session info for reproducibility
sessionInfo()

# Snapshot the current state of the project
renv::snapshot()

# TODO: Add your analysis workflow below
# Consider using R Markdown for a reproducible report
# Example:
# rmarkdown::render("analysis_workflow.Rmd")

# TODO: Consider implementing a Snakemake or Nextflow pipeline for full reproducibility
