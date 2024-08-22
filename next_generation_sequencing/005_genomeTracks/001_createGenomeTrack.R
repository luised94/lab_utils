#DESCRIPTION: 
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
working_directory <- paste(Sys.getenv("HOME"), "data", "240808Bel", sep ="/")
documentation_dir <- paste(working_directory, "documentation", sep ="/")
feature_file_directory <- paste(Sys.getenv("HOME"), "data", "feature_files", sep = "/")
path_to_sample_info <- list.files(documentation_dir, pattern = "sample_table", full.names = TRUE)
df_sample_info <- as.data.frame(read.delim(path_to_sample_info, header = TRUE))
df_sample_info$sample_ID <- as.character(13119:13151)
#TODO: Define the comparisons and plots to be generated in my sampleConfig.R template. 
#TODO: Find the best way to have the same levels and factors when I read in the sampleGridConfig.R file. This is defined by the categories list variable. Can grab that during 002_loadSampleGrid and use it to equalize. 
#TODO: Use 002_loadSampleGrid to open the sample_table given a directory. Use conditional statement to determine the ID. Could potentially use system function to call find and print with awk statement. R would require list.files(), strsplit, and grabbing regular expression for 5 digits.
#Process G1 and G2 values into Alpha and Nocodazole
#df_sample_info$Cell_Cycle <- ifelse(df_sample_info$Cell_Cycle == "G1", "Alpha", 
#                             ifelse(df_sample_info$Cell_Cycle == "G2", "Nocodazole", 
#                             df_sample_info$Cell_Cycle))
#
#df_sample_info$antibody <- ifelse(df_sample_info$antibody == "174", "74", 
#                             ifelse(df_sample_info$antibody == "185", "85", 
#                             df_sample_info$antibody))
#df_sample_info$short_name <- do.call(paste0, lapply(df_sample_info[setdiff(names(df_sample_info), c("short_name", "sample_ID"))], substr, 1, 1))

directory_of_refgenomes <- paste(Sys.getenv("HOME"), "data", "REFGENS", sep = "/")

# Read in S. cerevisiae S288C reference genome and process into dataframe.
genome_file_path <- list.files(directory_of_refgenomes, pattern = "S288C_refgenome.fna", full.names = TRUE, recursive = TRUE)
df_sacCer_refGenome <- readFasta(genome_file_path)
df_sacCer_refGenome <- data.frame(chrom = names(as(df_sacCer_refGenome, "DNAStringSet")), 
                   basePairSize = width(df_sacCer_refGenome)) %>% filter(chrom != "chrM")

# Process chromosome names to turn into chr<num> format and create the chromosome identifier value.
#TODO: Dont have to do this since they were adjusted using a previous bash script. 
#parts_by_comma <- unlist(strsplit(df_sacCer_refGenome$chrom, ","))
#chromosome_names <- parts_by_comma[!grepl("complete", parts_by_comma)]
#chromosome_number <- sub(".*chromosome ", "", chromosome_names)
#df_sacCer_refGenome$chrom <- paste("chr", chromosome_number, sep = "") 
##TODO: Output this processed genome. So I can easily read it again.
#chromosome_ID <- unlist(lapply(strsplit(chromosome_names, " "), '[[', 1))
##Used to extract from Bigwig file
#df_sacCer_refGenome$chrom_ID <- chromosome_ID

options(ucscChromosomeNames=FALSE) # Has to be run every time if you are using chromosome ID to get tracks from the bigwig file.
bigwig_directory <- paste(working_directory, "bigwig", sep = "/")

# Create list with samples to plot. 
#TODO: need to automate the creation of the experiments to plot. 
#TODO: Need to ensure that the names I use in the columns are compatible with my short name convention for subsetting df

experiments_to_plot <- list(
    c("lnnnMI", "onnnMA", "lnnnMA", "lnnyMA"),
    c("onnnMI", "onnnG7", "on7nG7", "onnnM7", "on7nM7")
)

descriptive_names_for_plots <- c(
    "AlFA_comparison",
    "MCM7_MCMpoly"
)
#experiments_to_plot <- list(
#    c("lnnnMI", "onnnMA", "lnnnMA", "lnnyMA"),
#    c("lnnnMI", "onnnMH", "lnnnMH", "lnnyMH"), 
#    c("lnnnMI", "lwnnG7", "lwnnM7"),
#    c("onnnMI", "onnnG7", "on7nG7", "onnnM7", "on7nM7"), 
#    c("onnnMI", "onnnG1", "on7nG1", "onnnM1", "on7nM1"), 
#    c("onnnMI", "onnnGC", "on7nGC", "onnnMC", "on7nMC"), 
#    c("lnnnMI", "lwnnG1", "lw2nG1", "lw2nM1"), 
#    c("lnnnMI", "lwnnGC", "lw2nGC", "lw2nMC"), 
#    c("lnnnMI", "lwnnG7", "lw2nG7", "lw2nM7")
#) 

#descriptive_names_for_plots <- c(
#    "AlFA_comparison",
#    "ORC_comparison", 
#    "Cell_Cycle_MCMpoly",
#    "MCM7tag_MCMpoly", 
#    "MCM7tag_HA11", 
#    "MCM7tag_HAC", 
#    "MCM2tag_HA11", 
#    "MCM2tag_HAC", 
#    "MCM2tag_MCMpoly"
#)

#Create GRanges object to read in a particular chromosome
#chromosome_to_plot <- 12
chromosome_to_plot <- 10
#chromosome_to_plot <- 14
genomeRange_to_get <- GRanges(seqnames=c(df_sacCer_refGenome$chrom[chromosome_to_plot]), 
        ranges = IRanges(start = 1, 
        end = df_sacCer_refGenome$basePairSize[chromosome_to_plot]), 
        strand = "*")

#Create variables to name plot
main_title_of_plot_track <- paste("Complete View of Chrom", as.character(chromosome_to_plot, 
                  sep = " "))
date_plot_created <- stringr::str_replace_all(Sys.time(), pattern = ":| |-", replacement="")  

plot_output_dir <- paste(working_directory, "plots", sep = "/")
for (experiment_index in 1:length(experiments_to_plot)) {
    all_tracks_to_plot <- list(GenomeAxisTrack(name = paste("Chr ", chromosome_to_plot, " Axis", sep = "") ) )
    #Subset the dataframe using strings    
    #TODO: Need some way to aggregate the column names I use from the sampleConfig.R files to keep track of variables I use to keep consistent experiment to experiment.
    #Files are already ordered from the sampleGridConfig.R file. No need to arrange.
    df_sample_info_subset <- df_sample_info %>% filter(short_name %in% unlist(experiments_to_plot[experiment_index])) #%>% arrange(antibody, rescue_allele, mcm_tag, cell_cycle, auxin_treatment)
    #print("Head output of subset dataframe")
    #print(head(df_sample_info_subset))
    #Create all of the tracks and append them to the all_tracks_to_plot variable
    for (sample_index in 1:nrow(df_sample_info_subset)) {
        #cat(sprintf("Sample index: %s", sample_index), "\n")
    #    cat(sprintf("Sample ID: %s", as.character(df_sample_info_subset$sample_ID[sample_index])), "\n")
        #cat(sprintf("Sample name: %s", as.character(df_sample_info_subset$short_name[sample_index])), "\n")
        #initial_matches <- list.files(bigwig_directory, pattern = as.character(df_sample_info_subset$sample_ID[sample_index]), full.names = TRUE, recursive = TRUE)     
        #path_to_bigwig <- initial_matches[grepl("S288C", initial_matches)]
        #bigwig_to_plot <- import(con = path_to_bigwig, which = genomeRange_to_get)
        #track_to_plot <- DataTrack(bigwig_to_plot, type = "l", name = df_sample_info_subset$short_name[sample_index], chromosome = df_sacCer_refGenome$chrom[chromosome_to_plot])
        #all_tracks_to_plot <- append(all_tracks_to_plot, track_to_plot)                                                                          
        #cat("Using sample index to subset")
        #cat(sprintf("Bigwig file to access: %s", path_to_bigwig), "\n")

        sample_ID_pattern <- df_sample_info_subset$sample_ID[df_sample_info_subset$short_name == unlist(experiments_to_plot[experiment_index])[sample_index]]
        initial_matches <- list.files(bigwig_directory, pattern = as.character(sample_ID_pattern), full.names = TRUE, recursive = TRUE)     
        path_to_bigwig <- initial_matches[grepl("S288C", initial_matches)]
        bigwig_to_plot <- import(con = path_to_bigwig, which = genomeRange_to_get)
        sample_short_name <- df_sample_info_subset$short_name[df_sample_info_subset$short_name == unlist(experiments_to_plot[experiment_index])[sample_index]]
        track_to_plot <- DataTrack(bigwig_to_plot, type = "l", name = sample_short_name, chromosome = df_sacCer_refGenome$chrom[chromosome_to_plot])
        all_tracks_to_plot <- append(all_tracks_to_plot, track_to_plot)                                                                          

        cat(sprintf("Experiment_to_plot sample: %s", unlist(experiments_to_plot[experiment_index])[sample_index]), "\n")
        cat(sprintf("Sample Short_name: %s", sample_short_name), "\n")
        cat(sprintf("Sample ID_pattern: %s", sample_ID_pattern), "\n")
        cat("Using experiments_to_plot subset","\n")
        cat(sprintf("Bigwig file to access: %s", path_to_bigwig), "\n")
    }
#TODO figure out if there is a way to normalize the samples 
#    # Determine the scale of the plotTracks plot by getting the max 
#    MAX <- -Inf
#    for (track in 1:length(all_tracks_to_plot)) {
#        print("Name of track")
#        print(names(all_tracks_to_plot[[track]]))
#      if(class(all_tracks_to_plot[[track]]) != "GenomeAxisTrack"){
#        if(max(all_tracks_to_plot[track][[1]]@data) > MAX) MAX <- max(all_tracks_to_plot[track][[1]]@data)
#      }
#    }
    #Generate the plot
#    svg(paste(plot_output_dir, "/", date_plot_created, "_", chromosome_to_plot, "_", descriptive_names_for_plots[experiment_index], ".svg", sep = ""))
#    plotTracks(all_tracks_to_plot, main = main_title_of_plot_track, chromosome = df_sacCer_refGenome$chrom[chromosome_to_plot], ylim = c(0, MAX))
#    dev.off()
}
for (sample_index in 1:nrow(df_sample_info)) {
    all_tracks_to_plot <- list(GenomeAxisTrack(name = paste("Chr ", chromosome_to_plot, " Axis", sep = "") ) )
    sample_ID_pattern <- df_sample_info$sample_ID[sample_index]
    initial_matches <- list.files(bigwig_directory, pattern = as.character(sample_ID_pattern), full.names = TRUE, recursive = TRUE)     
    print(initial_matches[1])
    path_to_bigwig <- initial_matches[grepl("S288C", initial_matches)]
    print("Name of the bigwig path")
    print(path_to_bigwig)
    if (length(path_to_bigwig) > 0){
        bigwig_to_plot <- import(con = path_to_bigwig, which = genomeRange_to_get)
        sample_short_name <- df_sample_info$short_name[sample_index]
        track_to_plot <- DataTrack(bigwig_to_plot, type = "l", name = sample_short_name, chromosome = df_sacCer_refGenome$chrom[chromosome_to_plot])
        all_tracks_to_plot <- append(all_tracks_to_plot, track_to_plot)                                                                          
    }
    #Generate the plot
    print("Name of the plot to be generated")
    print(paste(plot_output_dir, "/", date_plot_created, "_", chromosome_to_plot, "_", df_sample_info$short_name[sample_index], ".svg", sep = ""))
    #svg(paste(plot_output_dir, "/", date_plot_created, "_", chromosome_to_plot, "_", df_sample_info$short_name[sample_short_name], ".svg", sep = ""))
    #plotTracks(all_tracks_to_plot, main = main_title_of_plot_track, chromosome = df_sacCer_refGenome$chrom[chromosome_to_plot], ylim = c(0, 100000))
    #dev.off()

}
