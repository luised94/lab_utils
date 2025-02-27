library(dplyr)
library(tidyverse)
library(ggplot2)
library(xlsx)
library(gridExtra)
library(ggpubr)
library(plyr)
library(raster)

#Paths need to be updated
path_1 = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/Query Output/ARS_info_2020-04-08.csv"
path_2 = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/Query Output/ORF_info_2020-04-06.csv"
path_3 = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/Query Output/Chr_info_2020-04-06.csv"

read.csv(path_1, header = TRUE) %>% mutate_if(is.factor, as.character) -> ARS_df
read.csv(path_2, header = TRUE) %>% mutate_if(is.factor, as.character) -> ORF_df
read.csv(path_3, header = TRUE) -> CHR_df

ARS_df$ChrNum <- as.factor(ARS_df$ChrNum)
ARS_df <- arrange(ARS_df, ChrNum, Genome_Start)

ORF_df$ChrNum <- as.factor(ORF_df$ChrNum)
ORF_df <- filter(ORF_df, ORF_df$ChrNum != "chrmt")
ORF_df <- arrange(ORF_df, ChrNum, Genome_Start)
ORF_df$ChrNum <- factor(ORF_df$ChrNum)

# str(ARS_df)
# str(ORF_df)

#Determine # of ARS per Chromosome 
# arsperChr <- ARS_df %>% group_by(ChrNum) %>% count() 
# plotting <- ARS_df %>% ggplot(aes(ChrNum, y = ..count..)) + geom_bar() + theme_classic() + labs(title = "#ARS per Chromosome", x = "Chrom Num", y = "# of ARS")
# plotting
#Initialize array with # of elements = to # of ARS and then create a column to put gene hits
ARS_array <- array(data = NA, length(ARS_df$ChrNum))
ARS_df['Gene_Upstream'] <- ARS_array
ARS_df['TXN_Direction_Upstream'] <- ARS_array
ARS_df['Gene_Downstream'] <- ARS_array
ARS_df['TXN_Direction_Downstream'] <- ARS_array

#Create variables to use in for loop
ARS_str_1000 = ARS_df$Genome_Start - 1000
ARS_end_1000 = ARS_df$Genome_end + 1000
Genome_St = ORF_df$Genome_Start
Genome_End = ORF_df$Genome_End

#For all of the ARS, determine genes that are in the same chromosome, then compare Start and End positions with ARS gene start and 
#end positions then filter ORF dataframe and select Gene Name and Strand for further quantification.
for (ars in 1:length(ARS_array)){
  same_chromosome <- ORF_df$ChrNum == ARS_df$ChrNum[ars] #determine genes in same chromosome to current ARS
  #logic for determining whether gene is upstream accounting for directionality
  is_upstream_plus = Genome_St < ARS_str_1000[ars] & Genome_End > ARS_str_1000[ars] & Genome_St < ARS_end_1000[ars] & 
    Genome_End < ARS_end_1000[ars] 
  is_upstream_minus =  Genome_St > ARS_str_1000[ars] & Genome_End < ARS_str_1000[ars] & Genome_St < ARS_end_1000[ars] & 
    Genome_End < ARS_end_1000[ars]
  #filter dataframe to gene that is in front of ARS
  gene_holder_upstream <- ORF_df %>% filter((same_chromosome) & (is_upstream_plus | is_upstream_minus)) %>% dplyr::select(GeneName,Strand)
  #if it is not empty, then add the value to the column of ARS_df
  if( !is_empty(gene_holder_upstream[[1]])){
    ARS_df$Gene_Upstream[[ars]] <- gene_holder_upstream['GeneName'][[1]]
    ARS_df$TXN_Direction_Upstream[[ars]] <- gene_holder_upstream['Strand'][[1]]
  }
  
  #logic for determining whether gene is upstream accounting for directionality
  is_downstream_plus = Genome_St > ARS_str_1000[ars] & Genome_End > ARS_str_1000[ars] & Genome_St < ARS_end_1000[ars] & 
    Genome_End > ARS_end_1000[ars]
  is_downstream_minus = Genome_St > ARS_str_1000[ars] & Genome_End > ARS_str_1000[ars] & Genome_St > ARS_end_1000[ars] &
    Genome_End < ARS_end_1000[ars]
  #filter dataframe to gene that is after the ARS 
  gene_holder_downstream <- ORF_df %>% filter((same_chromosome) & (is_downstream_plus | is_downstream_minus)) %>% dplyr::select(GeneName,Strand)
  #if it is not empty, then add the value to the column of ARS_df
    if( !is_empty(gene_holder_downstream[[1]])){
    ARS_df$Gene_Downstream[[ars]] <- gene_holder_downstream['GeneName'][[1]]
    ARS_df$TXN_Direction_Downstream[[ars]] <- gene_holder_downstream['Strand'][[1]]
    }
}

write.xlsx(ARS_df, "ARS_dataframe.xlsx")



total_ARS = length(ARS_df$DBID)
ARS_gene_none <- sum(is.na(ARS_df$TXN_Direction_Downstream) & is.na(ARS_df$TXN_Direction_Upstream))/total_ARS
ARS_gene_upstream <- sum(!is.na(ARS_df$TXN_Direction_Upstream) & is.na(ARS_df$TXN_Direction_Downstream))/total_ARS
ARS_gene_downstream <- sum(!is.na(ARS_df$TXN_Direction_Downstream) & is.na(ARS_df$TXN_Direction_Upstream))/total_ARS


ARS_df$TXN_Direction_Upstream[is.na(ARS_df$TXN_Direction_Upstream)] <- 2
ARS_df$TXN_Direction_Downstream[is.na(ARS_df$TXN_Direction_Downstream)] <- 2
ARS_df$Gene_Upstream[is.na(ARS_df$Gene_Upstream)] <- 2
ARS_df$Gene_Downstream[is.na(ARS_df$Gene_Downstream)] <- 2

ARS_df['ARS_TXN_type'] <- ARS_array

ARS_df$ARS_TXN_type[which(ARS_df$TXN_Direction_Upstream == 1 & ARS_df$TXN_Direction_Downstream == 1)] <- "Codirectional_1"
ARS_df$ARS_TXN_type[which(ARS_df$TXN_Direction_Upstream == -1 & ARS_df$TXN_Direction_Downstream == -1)] <- "Codirectional_2"
ARS_df$ARS_TXN_type[which(ARS_df$TXN_Direction_Upstream == 1 & ARS_df$TXN_Direction_Downstream == -1)] <- "Convergent"
ARS_df$ARS_TXN_type[which(ARS_df$TXN_Direction_Upstream == -1 & ARS_df$TXN_Direction_Downstream == 1)] <-  "Divergent"
ARS_df$ARS_TXN_type[which((ARS_df$TXN_Direction_Upstream == 1 | ARS_df$TXN_Direction_Upstream == -1) & ARS_df$TXN_Direction_Downstream == 2)] <- "Upstream"
ARS_df$ARS_TXN_type[which(ARS_df$TXN_Direction_Upstream == 2 & (ARS_df$TXN_Direction_Downstream == 1 | ARS_df$TXN_Direction_Downstream == -1))] <- "Downstream"
ARS_df$ARS_TXN_type[which(ARS_df$TXN_Direction_Upstream == 2 & ARS_df$TXN_Direction_Downstream == 2)] <- "None_Flanking"

rm(ars, ARS_array, ARS_end_1000,ARS_str_1000, Genome_End,Genome_St,gene_holder_downstream,gene_holder_upstream)
#At this point, I verified a few of the ARS to see if info was correct. The few I saw checked out. 

ARS_df[ARS_df == -1] <- 0 #change negative ones for strand directionality to 0 so that I can use mean to count proportions

total_ARS = length(ARS_df$DBID)
ARS_gene_none <- sum(is.na(ARS_df$TXN_Direction_Downstream) & is.na(ARS_df$TXN_Direction_Upstream))/total_ARS
ARS_gene_upstream <- sum(!is.na(ARS_df$TXN_Direction_Upstream) & is.na(ARS_df$TXN_Direction_Downstream))/total_ARS
ARS_gene_downstream <- sum(!is.na(ARS_df$TXN_Direction_Downstream) & is.na(ARS_df$TXN_Direction_Upstream))/total_ARS

ARS_df_2 <- ARS_df
ARS_df_2[is.na(ARS_df_2)] <- 2 #change to 2 so that i can compare using logical operators
convergent <- sum(ARS_df_2$TXN_Direction_Downstream == 0  & ARS_df_2$TXN_Direction_Upstream == 1)/total_ARS
divergent <- sum(ARS_df_2$TXN_Direction_Downstream == 1  & ARS_df_2$TXN_Direction_Upstream == 0)/total_ARS
co_directional_1 <- sum(ARS_df_2$TXN_Direction_Downstream == 1  & ARS_df_2$TXN_Direction_Upstream == 1)/total_ARS
co_directional_2 <- sum(ARS_df_2$TXN_Direction_Downstream == 0  & ARS_df_2$TXN_Direction_Upstream == 0)/total_ARS

#should add up to 1 if I have accounted for all possibilities
ARS_accountedfor <- ARS_gene_none + ARS_gene_downstream + ARS_gene_upstream + convergent + 
                    divergent + co_directional_1 + co_directional_2

#Create dataframe with Data for each category to be able to plot 
ARS_gene_df <- data.frame(matrix(ncol = 5, nrow = 6))
colnames(ARS_gene_df) <- c("Type_of_ARS", "Percent_of_Total_ARS", "ARS_List", "Gene_List_Up", "Gene_List_Down")
ARS_gene_df$Type_of_ARS <- c("None flanking", "Upstream", "Downstream", "Co-Directional", "Divergent", "Convergent")
ARS_gene_df$Percent_of_Total_ARS <- c(ARS_gene_none, ARS_gene_upstream, ARS_gene_downstream, (co_directional_1 + co_directional_2), divergent, convergent)
ARS_gene_df$Type_of_ARS <- as.factor(ARS_gene_df$Type_of_ARS)

#Creating list for each ARS-type with corresponding ARS names, ARS sequences and genes


ARS_gene_df$ARS_List[1] <- ARS_df_2 %>% filter(TXN_Direction_Upstream == 2 & TXN_Direction_Downstream == 2) %>% pull(SystematicName) %>% list() %>% unlist() %>% paste(collapse=', ')
ARS_gene_df$Gene_List_Up[1] <- ARS_df_2 %>% filter(TXN_Direction_Upstream == 2 & TXN_Direction_Downstream == 2) %>% pull(Gene_Upstream) %>% list() %>% unlist() %>% paste(collapse=', ')
ARS_gene_df$Gene_List_Down[1] <- ARS_df_2 %>% filter(TXN_Direction_Upstream == 2 & TXN_Direction_Downstream == 2) %>% pull(Gene_Downstream) %>% list() %>% unlist() %>% paste(collapse=', ')


ARS_gene_df$ARS_List[2] <- ARS_df_2 %>% filter(TXN_Direction_Upstream != 2 & TXN_Direction_Downstream == 2) %>% pull(SystematicName) %>% list() %>% unlist() %>% paste(collapse=', ')
ARS_gene_df$Gene_List_Up[2] <- ARS_df_2 %>% filter(TXN_Direction_Upstream != 2 & TXN_Direction_Downstream == 2) %>% pull(Gene_Upstream) %>% list() %>% unlist() %>% paste(collapse=', ')
ARS_gene_df$Gene_List_Down[2] <- ARS_df_2 %>% filter(TXN_Direction_Upstream != 2 & TXN_Direction_Downstream == 2) %>% pull(Gene_Downstream) %>% list() %>% unlist() %>% paste(collapse=', ')


ARS_gene_df$ARS_List[3] <- ARS_df_2 %>% filter(TXN_Direction_Upstream == 2 & TXN_Direction_Downstream != 2) %>% pull(SystematicName) %>% list() %>% unlist() %>% paste(collapse=', ')
ARS_gene_df$Gene_List_Up[3] <- ARS_df_2 %>% filter(TXN_Direction_Upstream == 2 & TXN_Direction_Downstream != 2) %>% pull(Gene_Upstream) %>% list()%>% unlist() %>% paste(collapse=', ')
ARS_gene_df$Gene_List_Down[3] <- ARS_df_2 %>% filter(TXN_Direction_Upstream == 2 & TXN_Direction_Downstream != 2) %>% pull(Gene_Downstream) %>% list() %>% unlist() %>% paste(collapse=', ')


ARS_gene_df$ARS_List[4] <- ARS_df_2 %>% filter((TXN_Direction_Upstream == TXN_Direction_Downstream) & (TXN_Direction_Upstream != 2)) %>% pull(SystematicName) %>% list() %>% unlist() %>% paste(collapse=', ')
ARS_gene_df$Gene_List_Up[4] <- ARS_df_2 %>% filter((TXN_Direction_Upstream == TXN_Direction_Downstream) & (TXN_Direction_Upstream != 2)) %>% pull(Gene_Upstream) %>% list() %>% unlist() %>% paste(collapse=', ')
ARS_gene_df$Gene_List_Down[4] <- ARS_df_2 %>% filter((TXN_Direction_Upstream == TXN_Direction_Downstream) & (TXN_Direction_Upstream != 2)) %>% pull(Gene_Downstream) %>% list() %>% unlist() %>% paste(collapse=', ')
 

ARS_gene_df$ARS_List[5] <- ARS_df_2 %>% filter(TXN_Direction_Upstream == 0 & TXN_Direction_Downstream == 1) %>% pull(SystematicName) %>% list() %>% unlist() %>% paste(collapse=', ')
ARS_gene_df$Gene_List_Up[5] <- ARS_df_2 %>% filter(TXN_Direction_Upstream == 0 & TXN_Direction_Downstream == 1) %>% pull(Gene_Upstream) %>% list() %>% unlist() %>% paste(collapse=', ')
ARS_gene_df$Gene_List_Down[5] <- ARS_df_2 %>% filter(TXN_Direction_Upstream == 0 & TXN_Direction_Downstream == 1) %>% pull(Gene_Downstream) %>% list() %>% unlist() %>% paste(collapse=', ')


ARS_gene_df$ARS_List[6] <- ARS_df_2 %>% filter(TXN_Direction_Upstream == 1 & TXN_Direction_Downstream == 0) %>% pull(SystematicName) %>% list() %>% unlist() %>% paste(collapse=', ')
ARS_gene_df$Gene_List_Up[6] <- ARS_df_2 %>% filter(TXN_Direction_Upstream == 1 & TXN_Direction_Downstream == 0) %>% pull(Gene_Upstream) %>% list()%>% unlist() %>% paste(collapse=', ')
ARS_gene_df$Gene_List_Down[6] <- ARS_df_2 %>% filter(TXN_Direction_Upstream == 1 & TXN_Direction_Downstream == 0) %>% pull(Gene_Downstream) %>% list() %>% unlist() %>% paste(collapse=', ')

write.xlsx(ARS_gene_df, "ARS_transcription_overlap.xlsx")









convergent <- sum(ARS_df_2$TXN_Direction_Downstream == 0  & ARS_df_2$TXN_Direction_Upstream == 1)/total_ARS
divergent <- sum(ARS_df_2$TXN_Direction_Downstream == 1  & ARS_df_2$TXN_Direction_Upstream == 0)/total_ARS
co_directional_1 <- sum(ARS_df_2$TXN_Direction_Downstream == 1  & ARS_df_2$TXN_Direction_Upstream == 1)/total_ARS
co_directional_2 <- sum(ARS_df_2$TXN_Direction_Downstream == 0  & ARS_df_2$TXN_Direction_Upstream == 0)/total_ARS

#Plotting the %s of each type of gene-flanked ARS
ARS_gene_df %>% ggplot(aes(x = Type_of_ARS, Percent_of_Total_ARS)) + geom_col() + 
                geom_text(aes(label = round(Percent_of_Total_ARS, 3)), size = 4, position = position_stack(vjust = .5))+
                theme_classic()


 
#holder <- ARS_df_2 %>% filter(TXN_Direction_Upstream != 2) %>% pull(SystematicName)
#holder_1 <- ARS_df_2 %>% filter(TXN_Direction_Upstream != 2) %>% select(SystematicName, Gene_Upstream) %>% pull(SystematicName)
#colnames(dataframe)[which(names(dataframe) == "columnName")] <- "newColumnName"
# <- (sum(!is.na(ARS_df$TXN_Direction_Upstream) | !is.na(ARS_df$TXN_Direction_Downstream)) - 
#                      sum(!is.na(ARS_df$TXN_Direction_Upstream) & !is.na(ARS_df$TXN_Direction_Downstream)))/total_ARS
#values = c(ARS_gene_none, ARS_gene_upstream, ARS_gene_downstream, (co_directional_1 + co_directional_2), divergent, convergent)
#names_col = c("None flanking", "Upstream", "Downstream", "Co-Directional 1", "Co-Directional 2", "Divergent", "Convergent")
#ARS_gene_UpAndDown <- (sum(!is.na(ARS_df$TXN_Direction_Upstream) & !is.na(ARS_df$TXN_Direction_Downstream)))/total_ARS
#count_ARS_gene_UpAndDown <- sum(!is.na(ARS_df$TXN_Direction_Upstream) & !is.na(ARS_df$TXN_Direction_Downstream))
#Ignore these lines, they do the same as above but manually and without keeping it as data frame
# chromosomes <- ARS_df$ChrNum
# chromosomes <- chromosomes[!duplicated(chromosomes)]
# chr_num <- 1:length(chromosomes)
# ARSperChr <- integer(length(chromosomes))
# 
# for (Chr in chr_num){
#   ARSperChr[Chr] <- sum(ARS_df$ChrNum == chromosomes[Chr])
# }
# 
# plot(chr_num, ARSperChr)

#Lines used to figure out code logic and that it works
#Logic for gene being downstream or upstream from ARS also taking into account directionality of transcription. Dont think I needed to do this.
#Since information is encoded by strand parameter. 
# is_upstream_plus = Genome_St < ARS_str_1000[1] & Genome_End > ARS_str_1000[1] & Genome_St < ARS_end_1000[1] & Genome_End < ARS_end_1000[1]
# is_upstream_minus =  Genome_St > ARS_str_1000[1] & Genome_End < ARS_str_1000[1] & Genome_St < ARS_end_1000[1] & Genome_End < ARS_end_1000[1]
# is_downstream_plus = Genome_St > ARS_str_1000[1] & Genome_End > ARS_str_1000[1] & Genome_St < ARS_end_1000[1] & Genome_End > ARS_end_1000[1]
# is_downstream_minus = Genome_St > ARS_str_1000[1] & Genome_End > ARS_str_1000[1] & Genome_St > ARS_end_1000[1] & Genome_End < ARS_end_1000[1]
# 
# same_chromosome <- ORF_df$ChrNum == ARS_df$ChrNum[1]
# gene_holder <- ORF_df %>% filter((same_chromosome) & (is_upstream_plus | is_upstream_minus | is_downstream_plus | is_downstream_minus)) %>%
#                       select(GeneName,Strand)
# ARS_df['Gene_Upstream'] <- NA
# ARS_df$Gene_Upstream[[1]] <- gene_holder[[1]]
#require(ggplot2)
#qplot(Type_of_ARS, Percent_of_Total_ARS, data = ARS_gene_df, geom = "col")