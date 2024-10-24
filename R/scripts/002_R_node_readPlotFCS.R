# functions moved
# revisit
#library(flowCore)
#library(tidyverse)
#library(svglite)
#library(ggplot2)
#library(ggridges)
#
#path_to_fcs_files <- "/mnt/c/Users/Luis/Dropbox (MIT)/Lab/Experiments/Yeast Genetics/2023_10_04 Cdc6 Overexpression/CL4_Cdc6Oe_CarbonSourceShift_01"
##=======
#
##' Extract FITC Width Data from an FCS File
##'
##' This function reads an FCS file and extracts the fluorescence data
##' corresponding to the FITC width of individual cells.
##'
##' @param file_path The path to the FCS file.
##' @return A data frame with one column, FITC_Width, containing the extracted data.
##' @examples
##' extract_fitc_width_from_fcs("path/to/your/file.fcs")
##' @export
#extract_fitc_width_from_fcs <- function(file_path) {
#	    if (!file.exists(file_path)) {
#		            stop("File does not exist: ", file_path)
#    }
#    
#    fcs_data <- tryCatch({
#	            read.FCS(file_path, transformation = FALSE)
#		        }, error = function(e) {
#				        stop("Failed to read FCS file: ", e$message)
#		        })
#        
#        if (!'FITC-Width' %in% colnames(exprs(fcs_data))) {
#		        stop("FITC-Width data not found in the file.")
#	    }
#        
#        data_frame <- data.frame(FITC_Width = exprs(fcs_data[,'FITC-Width']))
#	    return(data_frame)
#}
#
##' Generate a Current Datetime String for File Naming
##'
##' This function returns a string representing the current date and time,
##' formatted for use in file names. The format used is "YYYY-MM-DD-HH-MM-SS",
##' and the system's local time zone is assumed.
##'
##' @return A character string of the current date and time.
##' @examples
##' get_current_datetime_string()
##' @export
#get_current_datetime_string <- function() {
#	          return(format(Sys.time(), "%Y-%m-%d-%H-%M-%S"))
#}
######
#
#windows_user <- list.files("/mnt/c/Users")[grepl(Sys.info()[["user"]], list.files("/mnt/c/Users"), ignore.case = TRUE)]
#path_to_fcs_files <- sprintf("/mnt/c/Users/%s/Dropbox (MIT)/Lab/Experiments/Yeast Genetics/2023_10_04 Cdc6 Overexpression/CL4_Cdc6Oe_CarbonSourceShift_01", windows_user)
#fcs_file_paths <- list.files(path_to_fcs_files, pattern = "\\.fcs$", full.names = TRUE)
##fcs_data <- read.FCS(fcs_file_paths[1])
#
##TODO: Initialize variables. for each one, use length to assign module class. then maye use switch or some other tidyverse to replace values with variables.
##REMINDER: Calculated the types of samples that would account for my FACS samples. Seems like I had 4 genotypes, 2 carbon sources and 7 timepoints.
#genotypes <- c("ORC4", "orc4-R267A", "ORC1", "orc1-K485T")
#carbon_sources <- c("GLU", "GAL")
#time_points <- c(seq(0, 150, by = 30))
#df_sample_info <- as.data.frame(expand.grid(genotype_Orc4 = genotypes, 
#					    carbonSources = carbon_sources, 
#					    timePoints = time_points)) %>% arrange(genotype_Orc4,carbonSources,timePoints)
#df_sample_info$filePath <- gtools::mixedsort(fcs_file_paths)#
#is_WTandGLU <- df_sample_info$genotype_Orc4 == "ORC4" & df_sample_info$carbonSources == "GLU"
#subset_df <- df_sample_info %>% filter(is_WTandGLU)
#for (i in 1:nrow(subset_df)) {
#
##Check the first six to ensure creation went well. 
#for (i in 1:6) {
##	flow_data <- read.FCS(df_sample_info$filePath[i])
##	exprs(flow_data[,'FITC-Width'])
#	formatted_output <- sprintf("%s | %s ", i, basename(subset_df$filePath[i]))
#	print(formatted_output)
#}
#
#subset_df$fitcData <- map(subset_df$filePath, extract_fitc_width_from_fcs)
#combined_data <- bind_rows(subset_df$fitcData, .id = "sample_id")
#
#head(combined_data)
#dim(combined_data)
#str(combined_data)
#summary(combined_data)
##Plot the data on the same plot and use color gradient to see. 
#output_directory <- sprintf("/mnt/c/Users/%s/Dropbox (MIT)/Lab/Experiments/Yeast Genetics/2023_10_04 Cdc6 Overexpression/", windows_user)
#plot_output <- paste0(output_directory, get_current_datetime_string(), "_densityPlot.svg")
##svglite(plot_output, width = 6, height = 4)
#ggplot(combined_data, aes(x = `FITC.Width`, y = after_stat(density), fill = sample_id)) +
#	geom_density(alpha = .2) +
#	xlim(0, median(combined_data$FITC.Width) + 1500) + 
#	labs(title = "Fluorescence Intensity Distribution", x = "Fluorescence Intensity", y = "Density")
##dev.off()
#
#
## Plot a staggered version using an offset of the samples on the same plot. 
#combined_data$offset <- as.numeric(as.factor(combined_data$sample_id)) * 100 # Create the offset by multiplying sample id with 10
#plot_output <- paste0(output_directory, get_current_datetime_string(), "_plot_count.svg")
