library(tidyverse)
library(ggplot2)
library(xlsx)
library(gridExtra)
library(ggpubr)
library(plyr)
library(raster)
source("C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/R Analysis Scripts/LEMR_plotting_functions.R")
#Summary of Results for the analysis done on this script (2020_6_2)
#Code to generate plots that lead to the following conclusions is at the bottom. 
#No relationship between distance to ARS or Intergenic Upstream Length to levels of IIA,B,D,TBP or Pol.II
#Increasing levels of IIA,IIB correlate with increasing TBP and Pol.II (for IIB)
#No relationship between mismatch to TATA o ARS_Type to levels of GTFs
#Genes that form divergent transcription structure around ARS seem to have bigger upstream intergenic regions
#Levels of IIA,IIB,IID,Pol.II,TBP may be lower (boxplot is narrower and lower) in genes facing away from ARS
#No relationship detected to between being associated to ARS or downstream/upstream.
#Nothing special gained by analyzing Rhee 2012 data set. If anything, suggests that suppression is due to transcription of target or targets
#Could be that has to be deeper analysis using other data to see how being close to ARS affects transcription


##Analysis requires dataframe produced by ARS_ORF_Analysis.R file!!!!!!!!
### Need this to adjust the levels of Chr factor for Hawkins Dataframe Chr column.
values <- c("chrI",    "chrII",  "chrIII",  "chrIV",   "chrV",    "chrVI",   "chrVII",  "chrVIII", "chrIX",   "chrX",   
            "chrXI",   "chrXII",  "chrXIII", "chrXIV",  "chrXV",   "chrXVI")
#Import Hawkins_df with info of closest ARS to position determined by fitting en Hawkins 2013 paper.
#Map values so Chr factor order is same order as chromosome number
path_to_Hawkins = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/Hawkins_Data/Hawkins_df_modified.xlsx"
read.xlsx(path_to_Hawkins, header = TRUE, sheetIndex = 1) %>% mutate_if(is.factor, as.character) -> Hawkins_df
Hawkins_df$Chr <- as.factor(Hawkins_df$Chr)
Hawkins_df$Chr <- factor(Hawkins_df$Chr, levels = values)
Hawkins_columns <- colnames(Hawkins_df)

#Import yeastmine query output for uncharacterized_verified ORFs
path_to_ORFs = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/Query Output/Verified_ORF_info_2020-06-01.xlsx"
read.xlsx(path_to_ORFs, header = TRUE, sheetIndex=1, rowIndex = NULL,endRow=5287) %>% mutate_if(is.factor, as.character) -> ORF_df
ORF_df$ChrNum <- as.factor(ORF_df$ChrNum)
ORF_df <- filter(ORF_df, ORF_df$ChrNum != "chrmt")
names(ORF_df)[names(ORF_df)=="ChrNum"] <- "Chr"
ORF_df <- arrange(ORF_df, Chr, Genome_Start)
ORF_df$Chr <- factor(ORF_df$Chr, levels = values)

#Read Rhee data from Rhee 2012 nature paper. Contains top 6045 genes with PIC signatures also has data 
# Retain rows and columns with data since read.xlsx imports NA columns and rows.
#Change some of the columns to factors.
#Change name of chrom columns to Chr so that access is same between df. 
#mapvalues of Chr column to Hawkins to reorder factor levels.
path_to_Rhee = 'C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/Rhee_Data/Top_PIC_genes_modified.xls'
read.xlsx(path_to_Rhee, header = TRUE, sheetIndex = 1) %>% mutate_if(is.factor, as.character) -> Rhee_df
Rhee_df <- Rhee_df[1:6045,1:25]
names(Rhee_df)[names(Rhee_df)=="chrom"] <- "Chr"
Rhee_columns <- colnames(Rhee_df)
Rhee_df %>% mutate_at(Rhee_columns[c(2,3,6,8,18,19)], factor) -> Rhee_df
Rhee_df$Chr <- mapvalues(Rhee_df$Chr, from = levels(Rhee_df$Chr), to = levels(Hawkins_df$Chr))
Rhee_df %>% mutate_at(Rhee_columns[c(5,24,25)], as.integer) -> Rhee_df


#Create empty array to fill Distance to ARS column (or other columns to be created and filled by for loop)
empty_array <- array(data = NA, length(Rhee_df$Chr))
Rhee_df['Dist_to_ARS'] <- empty_array
#For every gene, determine the minimum distance between its IIB_coordinate and the ARS that are in the same chromosome
#in Hawkins. 
for (gene in 1:length(Rhee_df$Chr)){
  #Filter df to ARS that are on same chromosome and store in df_holder
  df_holder <- Hawkins_df %>% filter(Chr == Rhee_df$Chr[gene]) 
  #Determine absolute bp distance between position of IIB and ARS
  dist_holder <- abs(df_holder$Position - Rhee_df$IIB_coor[gene])
  #Choose the minimum value in vector holding distances
  Rhee_df$Dist_to_ARS[gene] <- min(dist_holder)
}
rm(df_holder,dist_holder)
#Determine genes in Rhee_df that are flanking ARSs in Hawkins_df and assign the Architecture that they are part of
#Do first for genes that are upstream of ARS then for genes downstream.
Rhee_df['ARS_Type'] <-empty_array
Rhee_df['Location_to_ARS'] <- empty_array
idx_genes_flanking_up <- which(Rhee_df$gene_id %in% Hawkins_df$SystematicName_Upstream)
#Filter then sort since Rhee_df is already sorted by mRNA then gene_id.
df_holder <- Hawkins_df %>% filter(SystematicName_Upstream %in% Rhee_df[idx_genes_flanking_up,'gene_id']) %>% arrange(SystematicName_Upstream) %>% dplyr::select(SystematicName_Upstream,ARS_TXN_type) 
Rhee_df[idx_genes_flanking_up,27] <- df_holder[,2]
Rhee_df[idx_genes_flanking_up,28] <- "Upstream"

idx_genes_flanking_down <- which(Rhee_df$gene_id %in% Hawkins_df$SystematicName_Downstream)
df_holder <- Hawkins_df %>% filter(SystematicName_Downstream %in% Rhee_df[idx_genes_flanking_down,'gene_id']) %>% arrange(SystematicName_Downstream) %>% dplyr::select(SystematicName_Downstream,ARS_TXN_type) 
Rhee_df[idx_genes_flanking_down,27] <- df_holder[,2]
Rhee_df[idx_genes_flanking_down,28] <- "Downstream"
Rhee_df %>% filter(!is.na(ARS_Type)) %>% nrow()
rm(df_holder)
#Label the ones that aren't flanking the genes
Rhee_df[which(is.na(Rhee_df$ARS_Type)),27] <-"Not_Flanking"
Rhee_df[which(is.na(Rhee_df$Location_to_ARS)),28] <- "Not_Flanking"

Rhee_df$ARS_Type <- as.factor(Rhee_df$ARS_Type)
Rhee_df %>% mutate(ARS_flag = ifelse(ARS_Type == "Not_Flanking", "Not_Flanking", "Flanking")) -> Rhee_df


Rhee_df['Direction_to_ARS'] <- empty_array
idx_genes_up_to_ARS <- which(Rhee_df$gene_id %in% Hawkins_df$SystematicName_Upstream & Rhee_df$strand == "W")
Rhee_df[idx_genes_up_to_ARS,30] <- "Towards_ARS"

idx_genes_up_away_ARS <- which(Rhee_df$gene_id %in% Hawkins_df$SystematicName_Upstream & Rhee_df$strand == "C")
Rhee_df[idx_genes_up_away_ARS,30] <- "Away_ARS"

idx_genes_down_to_ARS <- which(Rhee_df$gene_id %in% Hawkins_df$SystematicName_Downstream & Rhee_df$strand == "C")
Rhee_df[idx_genes_down_to_ARS,30] <- "Towards_ARS"

idx_genes_down_away_ARS <- which(Rhee_df$gene_id %in% Hawkins_df$SystematicName_Downstream & Rhee_df$strand == "W")
Rhee_df[idx_genes_away_to_ARS,30] <- "Away_ARS"
Rhee_df[which(is.na(Rhee_df$Direction_to_ARS)),30] <-"Not_Flanking"

#For some reason, some genes arent assigned correctly their direction. Manually determined which ones. They were all downstream and on the watson strand. 
#assigned them to Away_ARS and made sure all were assinged.
Rhee_df %>% filter(ARS_flag == "Flanking" & Direction_to_ARS == "Not_Flanking") %>% dplyr::select(1,6,26,27,28,29,30)
Rhee_df[which(Rhee_df$ARS_flag == "Flanking" & Rhee_df$Direction_to_ARS == "Not_Flanking"),30] <- "Away_ARS"
Rhee_df %>% filter(ARS_flag == "Flanking" & Direction_to_ARS == "Not_Flanking") %>% nrow()

#Convert column type to factor for new columns
Rhee_columns <- colnames(Rhee_df)
Rhee_df %>% mutate_at(Rhee_columns[c(28:30)], factor) -> Rhee_df

###Generate plots for analysis###
Rhee_df %>% ggplot(aes(x = IIA, y=TBP)) + geom_point() + xlim(0, 4000)
Rhee_df %>% ggplot(aes(x = IIA, y=Pol.II)) + geom_point() + xlim(0, 4000)

Rhee_df %>% ggplot(aes(x = IIB, y=TBP)) + geom_point() 
Rhee_df %>% ggplot(aes(x = IIB, y=Pol.II)) + geom_point() 

Rhee_df %>% ggplot(aes(x = IID, y=TBP)) + geom_point() 
Rhee_df %>% ggplot(aes(x = IID, y=Pol.II)) + geom_point() 

Rhee_df %>% ggplot(aes(x = TBP, y=Pol.II)) + geom_point() 

Rhee_df %>% filter(!is.na(Intergenic_Up_Length)) %>% ggplot(aes(x = Intergenic_Up_Length, y=IIA)) + geom_point()
Rhee_df %>% filter(!is.na(Intergenic_Up_Length)) %>% ggplot(aes(x = Intergenic_Up_Length, y=IIB)) + geom_point()
Rhee_df %>% filter(!is.na(Intergenic_Up_Length)) %>% ggplot(aes(x = Intergenic_Up_Length, y=IID)) + geom_point()
Rhee_df %>% filter(!is.na(Intergenic_Up_Length)) %>% ggplot(aes(x = Intergenic_Up_Length, y=TBP)) + geom_point()
Rhee_df %>% filter(!is.na(Intergenic_Up_Length)) %>% ggplot(aes(x = Intergenic_Up_Length, y=Pol.II)) + geom_point()

Rhee_df %>%  ggplot(aes(x = Dist_to_ARS, y=IIA)) + geom_point()
Rhee_df %>%  ggplot(aes(x = Dist_to_ARS, y=IIB)) + geom_point()
Rhee_df %>%  ggplot(aes(x = Dist_to_ARS, y=IID)) + geom_point()
Rhee_df %>%  ggplot(aes(x = Dist_to_ARS, y=TBP)) + geom_point()
Rhee_df %>%  ggplot(aes(x = Dist_to_ARS, y=Pol.II)) + geom_point()

Rhee_df %>% filter(!is.na(Intergenic_Up_Length)) %>% ggplot(aes(x = Intergenic_Down_Length, y=TBP)) + geom_point()
Rhee_df %>% filter(!is.na(Intergenic_Up_Length)) %>% ggplot(aes(x = Intergenic_Down_Length, y=Pol.II)) + geom_point()

#Facetting by RNA
Rhee_df %>% ggplot(aes(x = IIA, y=TBP, color = RNA.class)) + geom_point() + xlim(0,1500) + facet_wrap(RNA.class~., nrow=2) + theme(legend.position = "none")
Rhee_df %>% ggplot(aes(x = IIA, y=Pol.II, color = RNA.class)) + geom_point() + xlim(0,1500) +facet_wrap(RNA.class~., nrow=2) + theme(legend.position = "none")

Rhee_df %>% ggplot(aes(x = IIB, y=TBP, color = RNA.class)) + geom_point() + facet_wrap(RNA.class~., nrow=2) + theme(legend.position = "none")
Rhee_df %>% ggplot(aes(x = IIB, y=Pol.II, color = RNA.class)) + geom_point() + facet_wrap(RNA.class~., nrow=2) + theme(legend.position = "none")

Rhee_df %>% ggplot(aes(x = IID, y=TBP, color = RNA.class)) + geom_point() + facet_wrap(RNA.class~., nrow=2) + theme(legend.position = "none")
Rhee_df %>% ggplot(aes(x = IID, y=Pol.II, color = RNA.class)) + geom_point() + facet_wrap(RNA.class~., nrow=2) + theme(legend.position = "none")

Rhee_df %>% ggplot(aes(x = TBP, y=Pol.II, color = RNA.class)) + geom_point() + facet_wrap(RNA.class~., nrow=2) + theme(legend.position = "none")

Rhee_df %>%  ggplot(aes(x = Dist_to_ARS, y=Pol.II)) + geom_point() + facet_wrap(RNA.class~., nrow = 2)
Rhee_df %>%  ggplot(aes(x = Dist_to_ARS, y=TBP)) + geom_point() + facet_wrap(RNA.class~., nrow = 2)

#Facetting by TATA type
Rhee_df %>% ggplot(aes(x = IIA, y=TBP, color = mismatch)) + geom_point() + xlim(0,2500)+facet_wrap(mismatch~., nrow=2) + theme(legend.position = "none")
Rhee_df %>% ggplot(aes(x = IIA, y=Pol.II, color = mismatch)) + geom_point() + xlim(0,2500)+facet_wrap(mismatch~., nrow=2) + theme(legend.position = "none")

Rhee_df %>% ggplot(aes(x = IIB, y=TBP, color = mismatch)) + geom_point() + facet_wrap(mismatch~., nrow=2) + theme(legend.position = "none")
Rhee_df %>% ggplot(aes(x = IIB, y=Pol.II, color = mismatch)) + geom_point() + facet_wrap(mismatch~., nrow=2) + theme(legend.position = "none")

Rhee_df %>% ggplot(aes(x = IID, y=TBP, color = mismatch)) + geom_point() + facet_wrap(mismatch~., nrow=2) + theme(legend.position = "none")
Rhee_df %>% ggplot(aes(x = IID, y=Pol.II, color = mismatch)) + geom_point() + facet_wrap(mismatch~., nrow=2) + theme(legend.position = "none")

Rhee_df %>% ggplot(aes(x = TBP, y=Pol.II, color = mismatch)) + geom_point() + facet_wrap(mismatch~., nrow=2) + theme(legend.position = "none")

Rhee_df %>% filter(!is.na(Intergenic_Up_Length)) %>% ggplot(aes(x = Intergenic_Up_Length, y=TBP, color = mismatch)) + geom_point() + facet_wrap(mismatch~., nrow=2) + theme(legend.position = "none")
Rhee_df %>% filter(!is.na(Intergenic_Up_Length)) %>% ggplot(aes(x = Intergenic_Up_Length, y=Pol.II, color = mismatch)) + geom_point() + facet_wrap(mismatch~., nrow=2) + theme(legend.position = "none")
#Facetting by ARS_Type
Rhee_df %>% ggplot(aes(x = IIA, y=TBP, color = ARS_Type)) + geom_point() +xlim(0,2500)+facet_wrap(ARS_Type~., nrow=2) + theme(legend.position = "none")

Rhee_df %>% ggplot(aes(x = ARS_Type, y = IIA, color = ARS_Type)) + 
  geom_boxplot() + geom_jitter(alpha= .4)+ ylim(0,1000)+theme_classic() + theme(legend.position = "none") 

Rhee_df %>% ggplot(aes(x = ARS_Type, y = IIB, color = ARS_Type)) + 
  geom_boxplot() + geom_jitter(alpha= .4)+ ylim(0,5000)+theme_classic() + theme(legend.position = "none") 

Rhee_df %>% ggplot(aes(x = ARS_Type, y = IID, color = ARS_Type)) + 
  geom_boxplot() + geom_jitter(alpha= .4)+ ylim(0,1800)+theme_classic() + theme(legend.position = "none") 

Rhee_df %>% ggplot(aes(x = ARS_Type, y = TBP, color = ARS_Type)) + 
  geom_boxplot() + geom_jitter(alpha= .4)+ ylim(0,600)+ theme_classic() + theme(legend.position = "none") 

Rhee_df %>% ggplot(aes(x = ARS_Type, y = Pol.II, color = ARS_Type)) + 
  geom_boxplot() + geom_jitter(alpha= .4)+ ylim(0,500)+ theme_classic() + theme(legend.position = "none") 

Rhee_df %>% filter(!is.na(Intergenic_Up_Length)) %>% ggplot(aes(x = ARS_Type, y = Intergenic_Up_Length, color = ARS_Type)) + 
  geom_boxplot() + geom_jitter(alpha= .4)+ theme_classic() + theme(legend.position = "none") 

Rhee_df %>% filter(RNA.class != "rRNA")%>%ggplot(aes(x = RNA.class, y = IIA, color = RNA.class)) + 
  geom_boxplot() + geom_jitter(alpha= .4)+ theme_classic() + theme(legend.position = "none") 

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = mismatch, y = IIA, color = mismatch)) + 
  geom_boxplot() + ylim(0,2500)+theme_classic() + theme(legend.position = "none") 


#Looking at broader classifications
#Whether genes is upstream or downstream
#Whether gene is facing towards or away from ARS
#Only look at mRNA
Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = Location_to_ARS, y = IIA, color = Location_to_ARS)) + 
  geom_boxplot() + ylim(0,500)+theme_classic() + theme(legend.position = "none")  

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = Location_to_ARS, y = TBP, color = Location_to_ARS)) + 
  geom_boxplot() +ylim(0,200)+theme_classic() + theme(legend.position = "none") 

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = Location_to_ARS, y = Pol.II, color = Location_to_ARS)) + 
  geom_boxplot() + ylim(0,200)+theme_classic() + theme(legend.position = "none") 

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = Location_to_ARS, y = IIB, color = Location_to_ARS)) + 
  geom_boxplot() + ylim(0,500)+theme_classic() + theme(legend.position = "none") 

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = Location_to_ARS, y = IID, color = Location_to_ARS)) + 
  geom_boxplot() + ylim(0,500)+theme_classic() + theme(legend.position = "none") 

#Next to origin or not
Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = ARS_flag, y = IIA, color = ARS_flag)) + 
  geom_boxplot() + ylim(0,500)+theme_classic() + theme(legend.position = "none")  

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = ARS_flag, y = TBP, color = ARS_flag)) + 
  geom_boxplot() +ylim(0,200)+theme_classic() + theme(legend.position = "none") 

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = ARS_flag, y = Pol.II, color = ARS_flag)) + 
  geom_boxplot() + ylim(0,200)+theme_classic() + theme(legend.position = "none") 

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = ARS_flag, y = IIB, color = ARS_flag)) + 
  geom_boxplot() + ylim(0,500)+theme_classic() + theme(legend.position = "none") 

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = ARS_flag, y = IID, color = ARS_flag)) + 
  geom_boxplot() + ylim(0,500)+theme_classic() + theme(legend.position = "none") 

#Towards or away from ARS
Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = Direction_to_ARS, y = IIA, color = Direction_to_ARS)) + 
  geom_boxplot() + ylim(0,200)+theme_classic() + theme(legend.position = "none")  

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = Direction_to_ARS, y = TBP, color = Direction_to_ARS)) + 
  geom_boxplot() +ylim(0,200)+theme_classic() + theme(legend.position = "none") 

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = Direction_to_ARS, y = Pol.II, color = Direction_to_ARS)) + 
  geom_boxplot() + ylim(0,200)+theme_classic() + theme(legend.position = "none") 

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = Direction_to_ARS, y = IIB, color = Direction_to_ARS)) + 
  geom_boxplot() + ylim(0,500)+theme_classic() + theme(legend.position = "none") 

Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = Direction_to_ARS, y = IID, color = Direction_to_ARS)) + 
  geom_boxplot() + ylim(0,500)+theme_classic() + theme(legend.position = "none")


Rhee_df %>% filter(RNA.class == "mRNA") %>%ggplot(aes(x = IIA, y=TBP, color = Direction_to_ARS)) + geom_point() +
  xlim(0,2500)+facet_wrap(Direction_to_ARS~.) + theme(legend.position = "none")
###End###
#Code for other analysis and figuring out how to do things

# Rhee_df$gene_id[idx_genes_flanking_up]
# idx_ars_flanking_up <- which(Hawkins_df$SystematicName_Upstream %in% Rhee_df[idx_genes_flanking,'gene_id'])
# idx_ars_flanking_down <- which(Hawkins_df$SystematicName_Downstream %in% Rhee_df[idx_genes_flanking, 'gene_id'])
# sort(Hawkins_df[idx_ars_flanking_up,14])
# order(Hawkins_df[idx_ars_flanking_up,14])
# length(idx_ars_flanking_up) + length(idx_ars_flanking_down)
# Hawkins_df[order(Hawkins_df[idx_ars_flanking_up,14]),24]
# Hawkins_df[order(Hawkins_df[idx_ars_flanking_up,14]),14]

#Import wang df to also look at relationship of PIC to type of ARS near gene.
#Retain columns with data and rename df columns
#Convert certain columns to factor, chatacters and reorder Chr column.

#Create a list that contains genes adjacent to ARSs for later analysis of the relationship between ARS adjacency 
#and PIC levels.
# ARS_adjacent = c()
# for (i in 1:length(Wang_df$ARS)){
#   #Processing that has to be done to create list. Genes have ; that separates genes next to ARS.
#   #Some just have one gene nearby. They have to processed differently.
#   #If they have ;, split string by that, append the first one after removing * character.
#   #then remove space character fpr the second one and append to list
#   if(grepl(";", Wang_df$Adjancent_Genes[i])){
#     split_string <- strsplit(Wang_df$Adjancent_Genes[i], ";")
#     ARS_adjacent <- append(ARS_adjacent, str_remove_all(split_string[[1]][1], fixed("*")))
#     holder <- str_remove_all(split_string[[1]][2], fixed("*"))
#     ARS_adjacent <- append(ARS_adjacent, str_remove_all(holder, fixed(" ")))
#   } else{ 
#     #If gene doesnt have ;, just remove * character and append to list. 
#     ARS_adjacent <- append(ARS_adjacent, str_remove_all(Wang_df$Adjancent_Genes[i], fixed("*")))
#   }
# }
# ARS_adjacent <- ARS_adjacent[1:699] #retain Genes only. Due to xlsx import, some NAs are added to list. 

#Import ARS_dataframe file that has associated genesto ARS based on my overlap analysis using ARS and ORFs from Intermine
#Used for same reason as Wang dataframe. 
#Retain columns with data and rename df columns
#Convert certain columns to factor, chatacters and reorder Chr column.
# path_to_ARS_Intermine = 'C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/ARS_dataframe.xlsx'
# read.xlsx(path_to_ARS_Intermine, header = TRUE, sheetIndex = 1) %>% mutate_if(is.factor, as.character) -> ARS_df
# ARS_df <- ARS_df[,2:13]
# names(ARS_df)[names(ARS_df)=="ChrNum"] <- "Chr"
# ARS_df_columns <- colnames(ARS_df)
# ARS_df %<>% mutate_at(ARS_df_columns[c(4,12)], factor)
# ARS_df$Chr <- mapvalues(ARS_df$Chr, from = levels(ARS_df$Chr), to = levels(Hawkins_df$Chr))

#Create column that represents whether gene is ARS_adjacent by checking whether gene_id or standard name
#is in list xtracted from Wang Analysis. 
# Rhee_df <- Rhee_df %>% mutate(ARS_Adjacent = ifelse(gene_id %in% ARS_adjacent | Standard_Name %in% ARS_adjacent, T, F))

#Convert all instances of "None" in Standard_Name column to "Not Available since None is in list for genes that I determined
#and would appear as positive. Also convert NA to "Not Available"
# Rhee_df$Standard_Name[Rhee_df$Standard_Name == "None"] <- "Not Available"
# Rhee_df$Standard_Name[is.na(Rhee_df$Standard_Name)] <- "Not Available"

#Create column that represents whether gene is ARS adjacent by checking whether gene_id or standard name
#is in Gene_Upstream or Gene_Downstream column of ARS_df that I made using Intermine comparisons. 
#THIS requires the dataframe created by ARS_and_ORF_data_analysis.R script.
# Rhee_df <- Rhee_df %>% mutate(ARS_Adjacent_2 = ifelse(Standard_Name %in% ARS_df$Gene_Upstream | Standard_Name %in% ARS_df$Gene_Downstream, T, F))

#Create column to label whether gene is adjacent to a certain type of ARS based on the number and ways gene surrounds it. 
# Rhee_df['ARS_Type'] <- empty_array

#For all rows in Rhee_df dataframe, see if it is in Gene_Upstream or Downstream,
#access the ARS_type column in ARS_df for the gene that is in one of the columns and
#put it in ARS_type of the Rhee_df. Gives error of length of replacement. 
#Have checked it 
# count= 0
# for (i in 1:length(Rhee_df$Standard_Name)){
#   if(Rhee_df$Standard_Name[i] %in% ARS_df$Gene_Upstream){
#     count <- count +1
#     Rhee_df$ARS_Type[i] <- as.character(ARS_df$ARS_TXN_type[which(Rhee_df$Standard_Name[i] == ARS_df$Gene_Upstream)])
#     
#   } else if (Rhee_df$Standard_Name[i] %in% ARS_df$Gene_Downstream){
#     Rhee_df$ARS_Type[i] <- as.character(ARS_df$ARS_TXN_type[which(Rhee_df$Standard_Name[i] == ARS_df$Gene_Downstream)])
#     count <- count +1
#   } else{
#     Rhee_df$ARS_Type[i] <- "Not_Adjacent"
#   }
# }
# #Turn to column to factor. 
# Rhee_df$ARS_Type <- as.factor(Rhee_df$ARS_Type)

# Rhee_df_1 <- Rhee_df
# Rhee_1_columns <- colnames(Rhee_df_1)
#  
# x_list_for_normalization = c(9:17,24)
# for (i in 1:length(x_list_for_normalization)){
#   name_for_column <- paste(Rhee_1_columns[x_list_for_normalization[i]], "Norm", sep="") 
#   normalized_data <- (Rhee_df_1[,x_list_for_normalization[i]] - min(Rhee_df_1[,x_list_for_normalization[i]]))/(max(Rhee_df_1[,x_list_for_normalization[i]]) - min(Rhee_df_1[,x_list_for_normalization[i]]))
#   Rhee_df_1[name_for_column] <- normalized_data
# }
# 
# x_list_for_standarization = c(9:17,24)
# for (i in 1:length(x_list_for_standarization)){
#   name_for_column <- paste(Rhee_1_columns[x_list_for_standarization[i]], "ZStd", sep="")
#   standarized_data <- (Rhee_df_1[,x_list_for_standarization[i]] - mean(Rhee_df_1[,x_list_for_standarization[i]], na.rm = TRUE))/sd(Rhee_df_1[,x_list_for_standarization[i]], na.rm = TRUE)
#   Rhee_df_1[name_for_column] <- standarized_data
# }
# rm(name_for_column, normalized_data, standarized_data)


#This part of the script produces plots that attempt to explore the relationship between ARS and Transcription,
#from the transcription stand-point ie are PIC levels related to being close to an ARS element. \
#This section uses functions found in LEMR_plotting_functions.R to make scatter, density, and box plot
#To do this, define x,y variables and index to facet (only works on factor columns using strings) 
#and filter by and finally index to color data. 
#Better to do targeted plots

# #View the distribution of PIC levels, distance to closest ARS across dataset
# x_list = c(9:17,24)
# plot_density_and_save_nocolor(Rhee_df_1, x_list)
# 
# #Look at relationship between PIC levels and to distance to ARS
# x_list = c(9:17,24)
# y_list = x_list
# plot_point_and_save_nocolor(Rhee_df_1, x_list, y_list)
# 
# #Look at density and scatter faceted by RNA class
# plot_density_and_save(Rhee_df_1, x_list, facet_index = 2, color_index = 2)
# plot_point_and_save(Rhee_df_1, x_list, y_list, facet_index = 2, color_index = 2)
# 
# #Look at density and scatter faceted by mismatch to TATA box
# plot_density_and_save(Rhee_df_1, x_list, facet_index = 8, color_index = 8)
# plot_point_and_save(Rhee_df_1, x_list, y_list, facet_index = 8, color_index = 8)
# 
# #By chromosome 
# plot_density_and_save_nocolor(Rhee_df_1, x_list, facet_index = 3)
# plot_point_and_save_nocolor(Rhee_df_1, x_list, y_list, facet_index = 3)
# 
# #By chromosome filtered to mRNA and colored by ARS_Type
# plot_density_and_save(Rhee_df_1, x_list, facet_index = 3, color_index = 27, filter_factor = 2, filter_factor ="mRNA")
# plot_point_and_save(Rhee_df_1, x_list, y_list, facet_index = 3, color_index = 27, filter_factor = 2, filter_factor ="mRNA")
# 
# #Look at density and scatter filtered to RNA y facet by ARS Adjancency 1 and 2, and ARS_type
# #ARS_Adjancency from Wang data
# plot_density_and_save(Rhee_df_1, x_list, facet_index = 25, color_index = 25, filter_index = 2, filter_factor = "mRNA")
# plot_point_and_save(Rhee_df_1, x_list, y_list, facet_index = 25, color_index = 25, filter_index = 2, filter_factor = "mRNA")
# 
# #ARS_Adjacency from Intermine ARS and ORF comparisons
# plot_density_and_save(Rhee_df_1, x_list, facet_index = 26, color_index = 26, filter_index = 2, filter_factor = "mRNA")
# plot_point_and_save(Rhee_df_1, x_list, y_list, facet_index = 26, color_index = 26, filter_index = 2, filter_factor = "mRNA")
# 
# #ARS_Type from Intermine ARS and ORF comparisons
# plot_density_and_save(Rhee_df_1, x_list, facet_index = 27, color_index = 26, filter_index = 2, filter_factor = "mRNA")
# plot_point_and_save(Rhee_df_1, x_list, y_list, facet_index = 27, color_index = 26, filter_index = 2, filter_factor = "mRNA")
# 
# #By chromosome and color by ARS_Type, filter to mRNA
# plot_density_and_save(Rhee_df_1, x_list, facet_index = 3, color_index = 27, filter_index = 2, filter_factor = "mRNA")
# plot_point_and_save(Rhee_df_1, x_list, y_list, facet_index = 3, color_index = 27, filter_index = 2, filter_factor = "mRNA")
# 
# #Compare distributions of PIC components for ARS_Adjacency columns and ARS_Type
# x_list = c(25,26,27)
# y_list = c(9:17,24)
# plot_box_and_save_nocolor(Rhee_df_1, x_list, y_list, filter_index = 2, filter_factor = "mRNA")
# 
# #Compare distributions and facet by mismatch to TATA
# plot_box_and_save_nocolor(Rhee_df_1, x_list, y_list, facet_index = 8,  filter_index = 2, filter_factor = "mRNA")
# 
# #Compare distributions and facet by Chr 
# plot_box_and_save_nocolor(Rhee_df_1, x_list, y_list, facet_index = 3,  filter_index = 2, filter_factor = "mRNA")

# x_list = c(25)
# y_list = c(9:17)
# facet_idx = 8
# plot_box_and_save(Rhee_df, x_list, y_list, facet_idx)
# 
# x_list = c(9:17,24)
# facet_idx = 8
# Rhee_df <- Rhee_df %>% arrange(ARS_Adjacent)
# plot_point_and_save(Rhee_df, x_list, y_list, facet_idx, color_index = 25)
# 
# 
#   
#   
# Rhee_columns <- colnames(Rhee_df)
# 
# x_list = c(10)
# y_list = c(9, 11:17)
# facet_idx = 24
# Rhee_df <- Rhee_df %>% arrange(ARS_Adjacent_2)
# 
# plot_point_and_save(Rhee_df, x_list, y_list, facet_idx, color_index = 25)
# 
# #These produce files with same name. Careful!
# plot_density_and_save(Rhee_df, x_list, facet_index = 0, color_index = 24)
# 
# plot_density_and_save(Rhee_df, x_list, facet_index = 0, color_index = 25)
# 
# Rhee_df %>% ggplot(aes(x = IIA, color = ARS_Type)) + geom_density(alpha =.4) + facet_wrap(ARS_Type~., nrow=2) + xlim(0, 750) +theme_classic()
# ggsave("Gene_IIA_density_by_ARS_Type_50rows_not_visualized.png")
# 
# Rhee_df %>% ggplot(aes(x = IIA, color = ARS_Adjacent_2)) + geom_density(alpha =.4) + xlim(0, 750) +theme_classic()
# ggsave("Gene_IIA_density_by_ARS_Adjacency_50rows_not_visualized.png")
# 
# x_list = c(24,25)
# y_list = c(9:17)
# facet_idx = 24
# filter_index = 0 
# plot_box_and_save_w_filter(Rhee_df, x_list, y_list)
# 
# plot_box_and_save_w_filter(Rhee_df, x_list, y_list, filter_index = 2, filter_factor = "mRNA")
# Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = IIA, color = ARS_Type)) + geom_density(alpha =.4) + facet_wrap(ARS_Type~., nrow=2) + xlim(0, 750) +theme_classic()
# ggsave("Gene_IIA_density_by_ARS_Adjacency_fitler_mRNA_24rows_not_visualized.png")
# 
# 
# 
# Rhee_df %>% ggplot(aes(x = IIA, color = numeric(0))) + geom_density(alpha =.4) + facet_wrap(ARS_Type~., nrow=2) + xlim(0, 750) +theme_classic()
# x_list = c(9:17)
# plot_density_and_save_w_filter(Rhee_df, x_list, facet_index = 24, color_index = 24)
# plot_density_and_save_w_filter(Rhee_df, x_list, facet_index = 25, color_index = 25)
# 
# plot_density_and_save_w_filter(Rhee_df, x_list, facet_index = 24, filter_index = 2, filter_factor ="mRNA",color_index = 24)
# plot_density_and_save_w_filter(Rhee_df, x_list, facet_index = 25, filter_index = 2, filter_factor ="mRNA",color_index = 25)
# 
# plot_density_and_save_w_filter(Rhee_df, x_list, facet_index = 24, color_index = 2)
# plot_density_and_save_w_filter(Rhee_df, x_list, facet_index = 25, color_index = 2)
# 
# 
# plot_point_and_save(Rhee_df, x_list, y_list, facet_idx)
# plot_density_and_save(Rhee_df, x_list, facet_idx)
# 

#Doing the same plot but for normalized and standarized data.
#Looks same. Just wanted to make sure. 
# x_list = c(28:37)
# y_list = x_list
# plot_point_and_save_nocolor(Rhee_df_1, x_list, y_list)
# 
# x_list = c(38:47)
# y_list = x_list
# plot_point_and_save_nocolor(Rhee_df_1, x_list, y_list)
# Rhee_df_1 <- subset(Rhee_df_1, select = -28)
# ARS_df %>% filter(Gene_Upstream == Gene_Downstream) %>% nrow()
# 
# Rhee_df <- Rhee_df %>% mutate(ARS_Type = ifelse(Standard_Name %in% ARS_df$Gene_Downstream, ARS_df$ARS_TXN_type[which(ARS_df$Gene_Downstream %in% Standard_Name[which(Standard_Name %in% ARS_df$Gene_Downstream)])], 
#                       ifelse(Standard_Name %in% ARS_df$Gene_Upstream, ARS_df$ARS_TXN_type[which(ARS_df$Gene_Upstream %in% Standard_Name[which(Standard_Name %in% ARS_df$Gene_Upstream)])], "None")))
# 
# ARS_df$ARS_TXN_type[which(ARS_df$Gene_Upstream %in% Rhee_df$Standard_Name[which(Rhee_df$Standard_Name %in% ARS_df$Gene_Downstream)])]
# 
# ARS_df$ARS_TXN_type[which(ARS_df$Gene_Downstream %in% Rhee_df$Standard_Name[which(Rhee_df$Standard_Name %in% ARS_df$Gene_Downstream)])]
# which(Rhee_df$Standard_Name[5716] %in% ARS_df$Gene_Upstream)
# path_to_Intermine <- 'C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/ARS_transcription_overlap.xlsx'
# read.xlsx(path_to_Intermine, header = TRUE, sheetIndex =1, ) %>% mutate_if(is.factor, as.character) -> TXN_Overlap_df
# TXN_Overlap_df <- TXN_Overlap_df[,2:6]
# 
# empty_array <- array(data = NA, length(TXN_Overlap_df$Type_of_ARS))
# TXN_Overlap_df['Gene_Up'] <- empty_array
# TXN_Overlap_df['Gene_Down'] <- empty_array

# 
# for (i in 1:length(c(10,24))){
#       Rhee_df %>% ggplot(aes_string(x = Rhee_columns[c(10,24)[i]], color = Rhee_columns[2])) + geom_density(alpha = .4, fill = "black") + facet_wrap(as.formula(paste(Rhee_columns[2],"~", ".", sep = "")), nrow = 2)
#       name_for_file <- paste(Rhee_columns[c(10,24)[i]], "density", ".png", sep = "_")
#       ggsave(name_for_file)
#     }
#   
# Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x=IIA, y = Pol.II, color = ARS_Adjacent)) + geom_point() 
# Rhee_df %>% filter(RNA.class == "mRNA") %>% ggplot(aes(x = ARS_Adjacent, y = IIA)) +
#   geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +facet_wrap(mismatch~.)
# 
# 
# Rhee_df['ARS_Adjacent_2'] <- empty_array
# Rhee_df['ARS_TXN_type'] <- empty_array
# Rhee_df['Loc_to_ARS'] <- empty_array
# Rhee_df$Standard_Name[Rhee_df$Standard_Name == "None"] <- "Not Avalaible"
# Rhee_df$Standard_Name[is.na(Rhee_df$Standard_Name)] <- "Not Available"
# #M
# for (i in 1:length(Rhee_df$gene_id)){
#   
#   for (j in 1:length(TXN_Overlap_df$Type_of_ARS)){
#     
#     if (grepl(Rhee_df$Standard_Name[i], TXN_Overlap_df$Gene_List_Up[j], fixed = TRUE)){
#       Rhee_df$ARS_Adjacent_2[i] <- TRUE
#       Rhee_df$ARS_TXN_type[i] <- TXN_Overlap_df$Type_of_ARS[j]
#       Rhee_df$Loc_to_ARS[i] <- "Up"
#       
#     } else if (grepl(Rhee_df$Standard_Name[i], TXN_Overlap_df$Gene_List_Down[j], fixed = TRUE)){
#       Rhee_df$ARS_Adjacent_2[i] <- TRUE
#       Rhee_df$ARS_TXN_type[i] <- TXN_Overlap_df$Type_of_ARS[j]
#       Rhee_df$Loc_to_ARS[i] <- "Down"
#     } 
#   }
# }
# Rhee_df$ARS_Adjacent_2[is.na(Rhee_df$ARS_Adjacent_2)] <- FALSE
# Rhee_df$ARS_TXN_type[is.na(Rhee_df$ARS_TXN_type)] <- "None"
# Rhee_df$Loc_to_ARS[is.na(Rhee_df$Loc_to_ARS)] <- "None"
# 
# # for (i in 1:length(Rhee_df$Standard_Name)){
# #    
# #   if( Rhee_df$gene_id[i] %in% ARS_adjacent | Rhee_df$Standard_Name[i] %in% ARS_adjacent){
# #     Rhee_df$ARS_Adjacent[i] <- TRUE
# #   } else{
# #     Rhee_df$ARS_Adjacent[i] <- FALSE
# #   }
# # }
# 
# # ARS_adjacent[1]
# # which(ARS_adjacent == ARS_adjacent[1])
# # 
# # b[[1]][1]
# # b <- strsplit(a, ";")
# # str_remove_all(a, ";")
# # c <- append(c, b)
# # 
# # for (i in 1:length(c(9:10))){
# #   Rhee_df %>% ggplot(aes_string(x = "Dist_to_ARS", y = Rhee_columns[c(9:10)[i]], color = "RNA.class")) + geom_point() + facet_wrap(as.formula(paste("RNA.class","~", ".", sep = "")), nrow = 2)
# #   ggsave(paste("Trial", as.character(i), ".png", sep = "_"), width = 13.7, height = 8.92)
# # }
# f# ARS_array <- array(data = NA, length(Hawkins_df$Chr))
# # Hawkins_df['ClosestARS'] <- ARS_array
# # Hawkins_df['ClosestARS_Start'] <- ARS_array
# # Hawkins_df['ClosestARS_End'] <-ARS_array
# # Hawkins_df['Distance_to_closest_ARS'] <- ARS_array
# # for (index in 1:length(Hawkins_df$Chr)){
# #   df_holder <- ARS_df %>% filter(Chr == Hawkins_df$Chr[index])
# #   start_holder <- abs(df_holder$Start - Hawkins_df$Position[index])
# #   end_holder <- abs(df_holder$End - Hawkins_df$Position[index])
# #   Hawkins_df$Distance_to_closest_ARS[index] <- min(c(min(start_holder), min(end_holder)))
# #   if (which.min(start_holder) == which.min(end_holder)){
# #     Hawkins_df$ClosestARS[index] <- df_holder$Standard_Name[which.min(start_holder)]
# #     Hawkins_df$ClosestARS_Start[index] <- df_holder$Start[which.min(start_holder)]
# #     Hawkins_df$ClosestARS_End[index] <- df_holder$End[which.min(start_holder)]
# #   } else{
# #     if(which.min(c(start_holder[which.min(start_holder)], end_holder[which.min(end_holder)])) == 1){
# #       Hawkins_df$ClosestARS[index] <- df_holder$Standard_Name[which.min(start_holder)]
# #       Hawkins_df$ClosestARS_Start[index] <- df_holder$Start[which.min(start_holder)]
# #       Hawkins_df$ClosestARS_End[index] <- df_holder$End[which.min(start_holder)]
# #     } else{
# #       Hawkins_df$ClosestARS[index] <- df_holder$Standard_Name[which.min(end_holder)]
# #       Hawkins_df$ClosestARS_Start[index] <- df_holder$Start[which.min(end_holder)]
# #       Hawkins_df$ClosestARS_End[index] <- df_holder$End[which.min(end_holder)]
# #     }
# #   }
# #   
# # }
# 
# 
# # c(9:17)
# 


