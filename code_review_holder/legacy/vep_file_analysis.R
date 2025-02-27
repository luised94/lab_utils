library(tidyverse)
library(plyr)
library(xlsx)
vep_file <- read.delim("C:/Users/liusm/Desktop/Exp6_vep_analysis.txt", sep = "\t")

vep_file <- vep_file[1:dim(vep_file)[1],c(2:7,10,15:21,29,33)]
vep_columns <- colnames(vep_file)
vep_file %>% mutate_at(vep_columns[c(3,4,7)], factor) -> vep_file
levels(vep_file[,3]) #Check levels to filter by consequence correctly. 
vep_file %>% filter(Consequence %in% levels(vep_file[,3])[c(3,5,6)] & Existing_variation == "-") -> vep_filt

vep_array <- array(data = NA, length(vep_filt$SIFT))
for (i in 1:length(vep_filt$SIFT)){
  parentheses_idx <- gregexpr(pattern ='(',vep_filt$SIFT[i], fixed=TRUE)[[1]][1]
  SIFT_string <- substring(vep_filt$SIFT[i], parentheses_idx, nchar(vep_filt$SIFT[i]))
  vep_array[i] <- as.numeric(substring(SIFT_string, 2,nchar(SIFT_string)-1)) 
}

vep_filt['SIFT_Score'] <- vep_array
vep_filt_sorted <- arrange(vep_filt, vep_filt$Existing_variation,vep_filt$IMPACT,vep_filt$SIFT_Score)[,c(3,4,5,6,8,10,11,13,15,16)]

write.xlsx(as.data.frame(vep_filt), "Exp6_top_variants.xlsx", row.names = FALSE)
