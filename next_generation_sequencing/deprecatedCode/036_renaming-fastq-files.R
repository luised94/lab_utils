#STATUS:
#Rename fastq files 

# fastq_ids<- read.table('./scripts/fastq_info.txt')
library(ShortRead)
library(tidyverse)
#Renaming files ---- 

fastq_ids <- read.table("../rscripts/fastq_info.txt") #or add an extra ../ depending on current working dir. Most scripts will be set relative to 221024Bel_CHIP
knitr::kable(fastq_ids)
#Obtain files that will be renamed 
fastq_files_to_rename <- list.files(path ='.', pattern = "*.fastq$", recursive = TRUE, full.names = TRUE)[!grepl("unmapped", list.files(path ='.', pattern = "*.fastq$", recursive = TRUE))]
print(fastq_files_to_rename)
#Subset the dataframe. I think it was recycling the TRUE values 
bmc_ids <- fastq_ids[!grepl("nnNnH", fastq_ids$sname),]
prefix <- './data/fastq-files/'
grepl()


for (i in 1:length(fastq_files_to_rename)) {
  correct_name <- bmc_ids$sname[as.logical(na.omit(unlist(lapply(bmc_ids$BMC_ID2, grepl, x = fastq_files_to_rename[i]))))]
  #correct_name <- any(grepl(fastq_ids$BMC_ID2[i], fastq_files_to_rename[grepl(fastq_ids$BMC_ID2[i], fastq_files_to_rename)]))
  print(correct_name)
  
  if (file.exists(fastq_files_to_rename[i])) {
    # print(i)
    print(fastq_files_to_rename[i])
    destination_file <- paste(prefix, correct_name, ".fastq", sep = '')
    print(destination_file)
    
    #File.rename will move the file to new location
    #file.rename(fastq_files_to_rename[i], replacement_name)
    #Cat should work but decided to use ShortRead package just in case
    #cat(fastq_files_to_rename[i], file =  correct_name, append = TRUE)
    #Using ShortRead package just in case
    ## open input stream
    stream <- open(FastqStreamer(fastq_files_to_rename[i]))
    on.exit(close(stream))

    repeat {
      ## input chunk
      tmp_fq <- yield(stream)

      #If input is 0, stop the function
      if (length(tmp_fq) == 0) {
        print("Processed all")
        break
      }
      writeFastq(
        tmp_fq,
        destination_file,
        mode = "a",
        compress = FALSE ,
        qualityType = "FastqQuality"
      )
    }
  }
}

#Get the newly created files.
fastq_files_to_analyze <- list.files(path ='./data/fastq-files', pattern = "*.fastq$", recursive = FALSE, full.names = TRUE)[!grepl("unmapped", list.files(path ='./data/fastq-files', pattern = "*.fastq$", recursive = FALSE))]


sum(as.numeric(format(file.size(fastq_files_to_rename), scientific=TRUE))) - sum(as.numeric(format(file.size(fastq_files_to_analyze), scientific=TRUE)))
sum(as.numeric(format(file.size(fastq_files_to_analyze), scientific=TRUE)))
# replacement_names <- c()
# for (i in 1:length(fastq_files_to_rename)){
#   
#   print(fastq_files_to_rename[i])
#   replacement_name <- paste(prefix, paste(fastq_ids$sname[i], ".fastq", sep = ''), sep = '')
#   print(replacement_name)
#   print(fastq_ids$V1[i])
#   replacement_names <- append(replacement_names, replacement_name)
#   print(i)
# }
#After running successfully, the gunzip command was used to unzip the eaton data
# fastq_ids$file_location <- list.files(path ='.', pattern = "*.fastq.gz|*.fastq$", recursive = TRUE, full.names = TRUE)[!grepl("unmapped", list.files(path ='.', pattern = "*.fastq.gz|*.fastq$", recursive = TRUE))]

