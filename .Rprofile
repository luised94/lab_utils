message("Sourcing .Rprofile...")
# Source renv when working on my local computer.
# Prevent renv activation while on cluster.
lscpu_output_chr <- system("lscpu", intern = TRUE)
cpu_info <- lscpu_output_chr[grep("Model name:", lscpu_output_chr)]
if(!grepl("Xeon", cpu_info))
{
  source("renv/activate.R")
}

# Set to use nvim or vim as the viewer of the help page
options(help_type = "text")
options(pager = file.path(Sys.getenv("HOME"), "lab_utils/bash/helper/nvim_R_pager.sh"))

# Remove previously defined objects to avoid polluting environment.
rm(list = ls())
