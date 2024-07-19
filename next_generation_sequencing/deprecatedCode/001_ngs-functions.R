#Functions used to create a pipeline for Next-Generation Sequencing analysis. 
#Includes filter and trim function, align, sort and index function. 
#Designed to be used with a list of files and the destination files for the output files of each function. 
#Also meant to be used with mapply to process multiple files continuously.

#Filter and trim function that requires the file to be processed (filelocation) and the file to be outputed
#Destination file (destination_file) is created if it is not missing. 
#Uses QuasR and ShortRead packages. Outputs uncompressed fastq files.

myFilterAndTrim <- function(filelocation, destination_file, folder = './'){
  # #check and install QuasR and ShortRead packages if not installed.
  # if (!require(QuasR) | !require(ShortRead)) {
  #   install.packages('QuasR', 'ShortRead')
  # }
  # 
  # #Load libraries that are used in function 
  # library(QuasR)
  # library(ShortRead)
  
  ##Create destination file, by default it is the current working directory
  #If it wasn't supplied, create the destination file based on which folder was provided 
  #Create based on whether the folder is the working directory or not
  if (missing(destination_file)){
    if (folder == './'){
      destination_file <- paste('./', sprintf("processed_%s", basename(filelocation )), sep = "")
    } else {
      destination_file <- paste(folder, sprintf("/processed_%s", basename(filelocation )), sep = "")
    }
  } 
    
  ## open input stream
  stream <- open(FastqStreamer(filelocation))
  on.exit(close(stream))
  
  repeat {
    ## input chunk
    tmp_fq <- yield(stream)
    
    #If input is 0, stop the function
    if (length(tmp_fq) == 0){
      print("Processed all")
      break
    }
    
    
    ##Trim reads that contain 'N'
    tmp_fq <- clean(tmp_fq)
    ##Filter duplicated reads 
    tmp_fq <- tmp_fq[!srduplicated(tmp_fq)]
    ##Trim reads
    tmp_fq <- trimEnds(tmp_fq, a ="5")
    ##Calculate meanq_reads to obtain reads with mean quality larger than 20
    meanq_reads <- rowSums(as(quality(tmp_fq), "matrix"), na.rm = TRUE)/width(tmp_fq)
    #meanq_reads[is.na(meanq_reads)] <- 0]
    #Had to include a filter for NA values. Not sure why some of them have it. 
    tmp_fq <- tmp_fq[meanq_reads[!is.na(meanq_reads)] >= 20]
    
    ## append to destination
    ##Do not compress, output quality as FastqQuality 
    writeFastq(tmp_fq, destination_file, "a", compress = FALSE , qualityType = "FastqQuality")
  }
  
  
}

#Includes an align sort and index function 
#Must provide most of the parameters required by bowtie2_samtools. Stops if they are not explicitly provided.
#Sorts and indexes files by default!
#Uses Rbowtie2 but also requires stringr.

alignSortandIndex <- function(index_dir, index_prefix, fastq_file, alignment_file, alignment_params, output_type = ".bam", sort_and_index=TRUE){
  #check and install QuasR and ShortRead packages if not installed.
  # if (!require(Rbowtie2)|!require(stringr)) {
  #   print("Rbowtie2 or stringr not installed. Installing with install.packages(c('Rbowtie2', 'stringr'))")
  #   install.packages(c('Rbowtie2','stringr'))
  # }
  
  #Load libraries that are used in function 
  # library(Rbowtie2)
  # library(stringr)
  
  #Check for input parameters.
  if (missing(index_dir)){
    stop("Execution stopped. Please provide diretory for index")
  } 
  if(missing(index_prefix)){
    stop("Execution stopped. Please provide index prefix")
  }
  if(missing(alignment_file)){
    stop("Execution stopped. Please provide alignment file path. If sam, specify output_type as .sam")
  }
  if(missing(alignment_params)){
    stop("Execution stopped. Please provide alignment parameters as a string.")
  }
  
  #Check that alignment file path and output_type are same
  if (!(paste(".", tools::file_ext(alignment_file),sep = "") == output_type)){
    stop("Execution stopped. Alignment file and output type is different. Please make sure alignment ends with .bam or .fastq")
  }
  
  #Align the files using bowtie2_samtools to output as BAM or SAM. Uses Rsamtools if samtools is not installed.
  bowtie2_samtools(bt2Index = file.path(tools::file_path_as_absolute(index_dir), index_prefix),
                                      output = str_replace(alignment_file, pattern = output_type, replacement = ""),
                                      outputType = str_replace(output_type, pattern = ".", replacement = ""),
                                      seq1 = fastq_file,
                                      seq2 = NULL,
                                      bamFile = NULL,
                                      alignment_params,
                                      overwrite = TRUE)
  
  #Sort and Index the created bam files for further analysis. Print if file is not to be indexed and sorted.
  if(sort_and_index){
    print("Sorting and Indexing")
    sortBam(alignment_file, destination = paste(tools::file_path_sans_ext(alignment_file), "_sorted", sep = ''))
    indexBam(paste(tools::file_path_sans_ext(alignment_file), "_sorted.bam", sep = ''))
  } else {
    print(paste(paste(tools::file_path_sans_ext(alignment_file), "_sorted.bam", sep = ''), "not sorted and indexed."))
  }
  
}

#Reads in using Rsamtools functions a bam file by chunking through the chromosomes
#Provide chromosome name and chromosome size
#Can modify the ScanBamParam to read other ways
bamReadPosAndQwidthByChromToRle <- function(chrom_name, bam_file, chromosome_size, ...)
{
  #Commented message used for diagnostic purposes
  #print(paste("Scanning", chrom_name, "Length is", chromosome_size, sep = " "))
  #Define parameters to read. See SAM manual for details. POS is first position. QWidth is length of read essentially
  param <- ScanBamParam(
    what = c('pos', 'qwidth'),
    which = GRanges(chrom_name, IRanges(1, chromosome_size)),
    flag = scanBamFlag(isUnmappedQuery = FALSE)
  )
  #Read in the bam file according to param
  x <- scanBam(bam_file, ..., param = param)[[1]]
  #Get coverage of an IRanges object constructed using  pos and width
  coverage(IRanges(x[["pos"]], width = x[["qwidth"]]))
}
#ckage_to_check <- c("QuasR", "GenomicAlignments","Gviz","rtracklayer","ShortRead")
 BASH_SECTION
#Rcript -e 'source("/home/luised94/data/rscripts/3-assign-directory-variables.r")'#Simple function that returns the time in year,month,date,hour,minute,second with underscores in between. Add to files names 
underscoreDate <- function(){
  return(stringr::str_replace_all(Sys.time(), pattern = ":| |-", replacement="_"))
}

loadPackages <- function(package_list){
  package_was_loaded <- unlist(suppressPackageStartupMessages(lapply(package_list, library, character.only = TRUE, logical.return=TRUE, quietly = TRUE)))
  if (length(lapply(package_list[!package_was_loaded], function(x){
    paste(x, "Package did not install")
  })) == 0) print("All packages loaded.")
}

assignChrTracks <- function(genome_df, index_, feature_df, range_){
  assign(genome_df$gr_var[index_], GRanges(seqnames=genome_df$chrom[index_], ranges = IRanges(start= 1, end = genome_df$size[index_]), strand = "*"), envir = .GlobalEnv)
  assign(genome_df$chr_df[index_], feature_df %>% filter(Chromosome == index_) %>% data.frame(), envir = .GlobalEnv)
  start <- get(genome_df$chr_df[index_]) %>% .$Position - range_
  end <- get(genome_df$chr_df[index_]) %>% .$Position + range_
  assign(genome_df$origin_gr_var[index_], GRanges(seqnames = genome_df$chrom[index_], ranges = IRanges(start = start, end = end), strand = "*"), envir = .GlobalEnv)
  assign(genome_df$origin_track_var[index_], AnnotationTrack(get(genome_df$origin_gr_var[index_]), name=paste(sacCer3_df$chrom[index_], "Origins", sep = " ")), envir = .GlobalEnv)
}
