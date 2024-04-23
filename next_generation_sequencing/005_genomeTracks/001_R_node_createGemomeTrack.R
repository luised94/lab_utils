#Load packages using through a list of strings and suppress the messages, return a TRUE if loading was succesful
package_list <- c("QuasR", "GenomicAlignments", "Gviz", "rtracklayer", "ShortRead", "tidyverse")
package_was_loaded <- unlist(
    suppressPackageStartupMessages(
        lapply(
            package_list, 
	    library, 
    	    character.only = TRUE, 
            logical.return=TRUE, 
	    quietly = TRUE
    )
  )
)

#Determine which packages were not loaded, print all packages loaded or which packages were not loaded. 
packages_not_loaded <- package_list[!package_was_loaded]
if (length(packages_not_loaded) == 0) {
    print("All packages loaded.")
} else {
  lapply(
    packages_not_loaded,
    function(x) { 
      message <- paste(x, "Package did not install") 
      print(message)
    }
  )
}

# 
# Create directory variables and dataframe with sample information.
working_directory <- paste(Sys.getenv("HOME"), "data", "240304Bel", sep ="/")
documentation_dir <- paste(working_directory, "documentation", sep ="/")
path_to_sample_info <- list.files(documentation_dir, pattern = "table", full.names = TRUE)
df_sample_info <- as.data.frame(read.csv(path_to_sample_info, 
				header = TRUE)
		  )

directory_of_refgenomes <- paste(Sys.getenv("HOME"), "data", "REFGENS", sep = "/")

# Read in S. cerevisiae S288C reference genome and process into dataframe.
genome_file_path <- list.files(directory_of_refgenomes, pattern = "S288C_refgenome.fna", full.names = TRUE, recursive = TRUE)
df_sacCer_refGenome <- readFasta(genome_file_path)
df_sacCer_refGenome <- data.frame(chrom = names(as(df_sacCer_refGenome, "DNAStringSet")), 
			       basePairSize = width(df_sacCer_refGenome)) %>% filter(chrom != "chrM")

# Process chromosome names to turn into chr<num> format and create the chromosome identifier value.
parts_by_comma <- unlist(strsplit(df_sacCer_refGenome$chrom, ","))
chromosome_names <- parts_by_comma[!grepl("complete", parts_by_comma)]
chromosome_number <- sub(".*chromosome ", "", chromosome_names)
df_sacCer_refGenome$chrom <- paste("chr", chromosome_number, sep = "") 

chromosome_ID <- unlist(lapply(strsplit(chromosome_names, " "), '[[', 1))
df_sacCer_refGenome$chrom_ID <- chromosome_ID

 
#Create GRanges object to read in a particular chromosome
chromosome_to_plot <- 14
genomeRange_to_get <- GRanges(seqnames=c(df_sacCer_refGenome$chrom_ID[chromosome_to_plot]), 
		ranges = IRanges(start = 1, 
		end = df_sacCer_refGenome$basePairSize[chromosome_to_plot]), 
		strand = "*")


options(ucscChromosomeNames=FALSE) # Has to be run every time if you are using chromosome ID to get tracks from the bigwig file.
bigwig_directory <- paste(working_directory, "bigwig", sep = "/")

# Create list with samples to plot. 
experiments_to_plot <- list(c("nnAYI", "nnNYV", "WnNYV"), c("nnAYI", "nnNNA", "nnNYA"), c("WnANI", "WnAYM", "WnNYM"), c("WnANI", "nnNYV", "WnNYV", "4nNYV", "44NYV"),  c("WnANI", "nnAYM", "WnAYM", "4nAYM", "44AYM"))
descriptive_names_for_plots <- c("controlV5", "controlAuxinTreatment", "controlMcmCellCycle", "comparingORCWT4R4PS", "comparinsMcmWT4R4PS")

for (experiment_index in 1:length(experiments_to_plot)) {
	all_tracks_to_plot <- list(GenomeAxisTrack(name = paste("Chr ", chromosome_to_plot, " Axis", sep = "") ) )
	df_sample_info_subset <- df_sample_info %>% filter(short_name %in% unlist(experiments_to_plot[experiment_index])) %>% arrange(antibody, ORC4_Rescue, Suppressor, Cell_Cycle, Auxin_treatment)
	#Create all of the tracks and append them to the all_tracks_to_plot variable
	for (sample_index in 1:nrow(df_sample_info_subset)) {
		initial_matches <- list.files(bigwig_directory, pattern = as.character(df_sample_info_subset$sample_ID[sample_index]), full.names = TRUE, recursive = TRUE) 	
		path_to_bigwig <- initial_matches[grepl("S288C", initial_matches)][1]
		bigwig_to_plot <- import(con = path_to_bigwig, which = genomeRange_to_get)
		track_to_plot <- DataTrack(bigwig_to_plot, type = "l", name = df_sample_info_subset$short_name[sample_index], chromosome = df_sacCer_refGenome$chrom_ID[chromosome_to_plot])
		all_tracks_to_plot <- append(all_tracks_to_plot, track_to_plot)                                                                          
	}
	# Determine the scale of the plotTracks plot by getting the max # TODO figure out if there is a way to normalize the samples 
	MAX <- -Inf
	for (track in 1:length(all_tracks_to_plot)) {
	  if(class(all_tracks_to_plot[[track]]) != "GenomeAxisTrack"){
	    if(max(all_tracks_to_plot[track][[1]]@data) > MAX) MAX <- max(all_tracks_to_plot[track][[1]]@data)
	  }
	}
	plot_output_dir <- paste(working_directory, "plots", sep = "/")
	date_file_created <- stringr::str_replace_all(Sys.time(), pattern = ":| |-", replacement="")  
	svg(paste(plot_output_dir, "/", date_file_created, "_", descriptive_names_for_plots[experiment_index], ".svg", sep = ""))
	plotTracks(all_tracks_to_plot, main = "Complete View of Chromosome 14", chromosome = df_sacCer_refGenome$chrom_ID[chromosome_to_plot], ylim = c(0, MAX))
	dev.off()
}

# Define path_to_bigwig for all of the files









# IGNORE FOR NOW
# TODO Figure out how to download hawkins timing data and then how to process it.  
#track_start <- 1
#track_end <- df_sacCer_refgenome$size[chromosome_to_plot]
#GRanges(seqnames = genome_df$chrom[index_], ranges = IRanges(start = , end = end), strand = "*"),
#AnnotationTrack(get(genome_df$origin_gr_var[index_]), name = "Origins"),
#
#assign(df_names[i], get(df_names[i]) %>% as.data.frame() %>% dplyr::select(1:7) %>% filter(!is.na(Chromosome)))
#feature_df %>% filter(Chromosome == index_) %>% data.frame(),














#Condition that determines if sample is
#is_input <- fastq_ids$antibody == 'input'
#
##Conditions that determine if sample is negative control
#is_negative<- (fastq_ids$antibody == 'V5' & fastq_ids$complement == 'none') |
#  (fastq_ids$antibody == 'Myc' & fastq_ids$auxin == 'yes') |
#  (fastq_ids$antibody == 'UM174' & fastq_ids$cellcycle == 'Noc') |
#  (fastq_ids$antibody == 'UM174' & fastq_ids$complement == 'none' & fastq_ids$auxin == 'yes')

#Creating a factor vector depending on the conditions.
#fastq_ids$Sample_type <-  factor(case_when(is_input ~ 'input',
#                 is_negative ~ 'negative',
#                 TRUE ~ 'experiment'))
#factor(case_when(as.numeric(rownames(fastq_ids)) <= 24 ~ 'A',
#                 as.numeric(rownames(fastq_ids)) > 24 ~ 'B'))
#
#fastq_ids <- fastq_ids %>% mutate(Pool = factor(case_when(as.numeric(rownames(fastq_ids)) <= 24 ~ 'A',
#                                          as.numeric(rownames(fastq_ids)) > 24 ~ 'B')))


# URL to download hawkins timing data
#hawkins_timing_url <- "https://ars.els-cdn.com/content/image/1-s2.0-S2211124713005834-mmc2.xlsx"

##### EXTRACT_TO_SNIPPET
#Define files that I want to load
#annofile_extension <- c("*.bed", "*.gff3", "*.xls", "*.tab")
#
##Create list, find files in feature_folder that contain annotations of interest
#feature_files <- c()
#for (extension in 1:length(annofile_extension)){
#  feature_files<- c(feature_files,
#                    list.files(feature_folder, pattern = annofile_extension[extension]))
#}
#  if (grepl("xls|xlsx", tools::file_ext(filepath))){
#    assign(df_names[i], readxl::read_excel(filepath))
#    if(df_names[i] == "timing"){
#      assign(df_names[i], get(df_names[i]) %>% as.data.frame() %>% dplyr::select(1:7) %>% filter(!is.na(Chromosome)))
#    }
#
#  } else if(tools::file_ext(filepath) == "bed"){
#    assign(df_names[i], readBed(filepath))
#    print(seqlevelsStyle(eval(parse(text = df_names[i]))))
#
#  } else if(tools::file_ext(filepath) == "gff3"){
#    assign(df_names[i], toGRanges(makeTxDbFromGFF(filepath), format = 'gene'))
#    print(seqlevelsStyle(eval(parse(text = df_names[i]))))
#
#  } else if(tools::file_ext(filepath) == "tab"){
#    assign(df_names[i], read.delim(filepath, header = FALSE))
#    if(df_names[i] == "SGD_features"){
#      assign(df_names[i], get(df_names[i]) %>% as.data.frame() %>% dplyr::select(c(9,10,11,12,2,4,5)))
#    }
#  }
#}
##### EXTRACT_TO_SNIPPET

##### EXTRACT_TO_SNIPPET
##Create the variables that will be assign the directory locations
##Some processing has to be done based on how they were created. Remove ./data/, -files and replace and - with _
#dir_variables <- paste(paste(stringr::str_replace(stringr::str_replace(list.dirs(data_folder, recursive = FALSE), pattern = "./data/", replace = ""), pattern = "-files", replace=""), "_directory", sep = ""))
#dir_variables <- stringr::str_replace(dir_variables, pattern = "-", "_")
##Create the data frame to iterate through
#dir_df <- data.frame(names = dir_names, variables = dir_variables)
#
##Initialize file extensions in same order as dir_variables based on initial script
#file_extension <- c(".bam", ".bw", ".fastq", ".html", ".fasta", ".bed", ".fastq")
#
##Create the variables depending on whether they exist or not
#if (!all(unlist(lapply(dir_variables, exists)))) {
#
#  dir_list <- c()
#  #For all of the directories, create the string variable with the path, assign it and add it to dir list
#  for (i in 1:length(dir_df$names)){
#
#    folder_variable <- paste(data_folder, dir_df$names[i], sep = "")
#    assign(dir_df$variables[i], folder_variable)
#    dir_list <- append(dir_list, folder_variable)
#  }
#
#  print("Directory variables created.")
#
#} else {
#
#  print("All variables assigned.")
#}
##### EXTRACT_TO_SNIPPET
##### EXTRACT_TO_SNIPPET
##Output bigwig files by chunking through chromosomes ----
##Uses mapply to go through bam files and their bw connections. Uses function bamReadPosAndQwidthByChromToRle mostly taken from the introduction to Rsamtools manual
##https://www.bioconductor.org/packages/release/bioc/vignettes/Rsamtools/inst/doc/Rsamtools-Overview.pdf
#
#mapply(function(bam_files, path_to_bigwignection) {
#  print(paste("Working on:", bam_files, sep = " "))
#  cvg <-
#    RleList(
#      mapply(
#        FUN = bamReadPosAndQwidthByChromToRle,
#        chrom_name = df_sacCer_refGenome_df$chrom,
#        bam_file = bam_files,
#        chromosome_size = df_sacCer_refGenome_df$size,
#        SIMPLIFY = FALSE
#      )
#    )
#  print(paste("Exporting:", path_to_bigwignection, sep = " "))
#  export.bw(con = path_to_bigwignection,
#            object = cvg,
#            format = "bw")
#}
#bamReadPosAndQwidthByChromToRle <- function(chrom_name, bam_file, chromosome_size, ...)
#{
#  #Commented message used for diagnostic purposes
#  #print(paste("Scanning", chrom_name, "Length is", chromosome_size, sep = " "))
#  #Define parameters to read. See SAM manual for details. POS is first position. QWidth is length of read essentially
#  param <- ScanBamParam(
#    what = c('pos', 'qwidth'),
#    which = GRanges(chrom_name, IRanges(1, chromosome_size)),
#    flag = scanBamFlag(isUnmappedQuery = FALSE)
#  )
#  #Read in the bam file according to param
#  x <- scanBam(bam_file, ..., param = param)[[1]]
#  #Get coverage of an IRanges object constructed using  pos and width
#  coverage(IRanges(x[["pos"]], width = x[["qwidth"]]))
#}
#
###### EXTRACT_TO_SNIPPET
# BASH_SECTION
 #Rcript -e 'source("/home/luised94/data/rscripts/3-assign-directory-variables.r")'
