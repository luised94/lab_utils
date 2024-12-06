#STATUS: REMOVE.
#Changing the name of the fastq files. 
#See C:\Users\Luis\Dropbox (MIT)\Lab\Experiments\ATP hydrolysis by ORC and Cdc6\ORC\ORC1_4 interface\ORC4 R267A\Assays\Yeast Genetics\2022_08_26_CHIP-seq_ORC-MCM-supps\downloading-bmc-data.md

library(tidyverse)

#Read in the file with the info for the fastq sequences
fastq_ids <- utils::read.table("C:/Users/Luis/Downloads/fastq-ids.tab")

#Files were moved to be in same folder as this script. Initially in my downloads folder. Initially had to be refered to in my downloads folder.
# fastq_ids <- bind_cols(utils::read.table("C:/Users/Luis/Downloads/fastq-ids.tab"), utils::read.table("C:/Users/Luis/Downloads/221024Bel_CHIP.tab"))
column_names <- c("BMC_ID1", "V2", "Index_Seq1", "BMC_ID2", "V3", "Index_Seq2")
names(fastq_ids) <- column_names


#Create the treatment vector and separate column two in treatment columns
treatments <- c("complement", "suppressor", "cellcycle", "auxin", "antibody", "index")
fastq_ids <- fastq_ids %>% separate("V2", treatments, sep = "_", remove = FALSE)

#Create a shorter named column that still contains treatment info. 
#Each letter will stand for the treatment used corresponding to the same order as treatment vector made above.
short_names <- str_c(str_sub(fastq_ids$complement,1,1), str_sub(fastq_ids$suppressor,1,1),str_sub(fastq_ids$cellcycle,1,1),str_sub(fastq_ids$auxin,1,1),str_sub(fastq_ids$antibody,1,1))
fastq_ids$sname <- short_names

#Condition that determines if sample is 
is_input <- fastq_ids$antibody == 'input'

#Conditions that determine if sample is negative control 
is_negative<- (fastq_ids$antibody == 'V5' & fastq_ids$complement == 'none') |
  (fastq_ids$antibody == 'Myc' & fastq_ids$auxin == 'yes') |
  (fastq_ids$antibody == 'UM174' & fastq_ids$cellcycle == 'Noc') |
  (fastq_ids$antibody == 'UM174' & fastq_ids$complement == 'none' & fastq_ids$auxin == 'yes')

#Creating a factor vector depending on the conditions.  
fastq_ids$Sample_type <-  factor(case_when(is_input ~ 'input',
                 is_negative ~ 'negative',
                 TRUE ~ 'experiment'))
factor(case_when(as.numeric(rownames(fastq_ids)) <= 24 ~ 'A',
                 as.numeric(rownames(fastq_ids)) > 24 ~ 'B'))

fastq_ids <- fastq_ids %>% mutate(Pool = factor(case_when(as.numeric(rownames(fastq_ids)) <= 24 ~ 'A',
                                          as.numeric(rownames(fastq_ids)) > 24 ~ 'B')))

#Creating the fastq rows for SRR files (combined by ) from eaton paper ----
fastq_ids <- fastq_ids %>% add_row(complement = c("none"),
                             suppressor = c("none"), 
                             cellcycle = c("Noc"),
                             auxin = c("n"),  
                             antibody = c("HM1108"),
                             index= c("49"), 
                             sname = c("nnNnH"),
                             Sample_type = c("eaton"),
                             Pool = c("C"),
                             V2 = c("none_none_Noc_no_HM1108"))

#Write dataframe to file to reference in cluster
write.table(fastq_ids, file = "./scripts/fastq_info.txt", sep = "\t")

#Read with function below
read.table('./scripts/fastq_info.txt')

# list.files(path ='.', pattern = "*.fastq$", recursive = TRUE)[!grepl("unmapped", list.files(path ='.', pattern = "*.fastq$", recursive = TRUE))]
