library(flowCore)
library(tidyverse)
library(svglite)
<<<<<<< HEAD
library(ggplot2)

path_to_fcs_files <- "/mnt/c/Users/Luis/Dropbox (MIT)/Lab/Experiments/Yeast Genetics/2023_10_04 Cdc6 Overexpression/CL4_Cdc6Oe_CarbonSourceShift_01"
fcs_file_paths <- list.files(path_to_fcs_files, pattern = "\\.fcs$", full.names = TRUE)
#fcs_data <- read.FCS(fcs_file_paths[1])

#TODO: Initialize variables. for each one, use length to assign module class. then maye use switch or some other tidyverse to replace values with variables.
#REMINDER: Calculated the types of samples that would account for my FACS samples. Seems like I had 4 genotypes, 2 carbon sources and 7 timepoints.
genotypes <- c("ORC4", "orc4-R267A", "ORC1", "orc1-K485T")
carbon_sources <- c("GLU", "GAL")
time_points <- c(seq(0, 150, by = 30))
df_sample_info <- as.data.frame(expand.grid(genotype_Orc4 = genotypes, 
					    carbonSources = carbon_sources, 
					    timePoints = time_points)) %>% arrange(genotype_Orc4,carbonSources,timePoints)
df_sample_info$filePath <- gtools::mixedsort(fcs_file_paths)#
is_WTandGLU <- df_sample_info$genotype_Orc4 == "ORC4" & df_sample_info$carbonSources == "GLU"
subset_df <- df_sample_info %>% filter(is_WTandGLU)
for (i in 1:nrow(subset_df)) {
#	flow_data <- read.FCS(df_sample_info$filePath[i])
#	exprs(flow_data[,'FITC-Width')
	formatted_output <- sprintf("%s | %s ", i, basename(subset_df$filePath[i]))
	print(formatted_output)
}
