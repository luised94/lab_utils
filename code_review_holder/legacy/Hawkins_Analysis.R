library(tidyverse)
library(ggplot2)
library(xlsx)
library(gridExtra)
library(ggpubr)
library(plyr)
library(raster)
source("C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/R Analysis Scripts/LEMR_plotting_functions.R")
### Need this to adjust the levels of Chr factor for Hawkins Dataframe Chr column.
values <- c("chrI",    "chrII",  "chrIII",  "chrIV",   "chrV",    "chrVI",   "chrVII",  "chrVIII", "chrIX",   "chrX",   
            "chrXI",   "chrXII",  "chrXIII", "chrXIV",  "chrXV",   "chrXVI")

path_1 = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/Hawkins_Data/Hawkins_df_with_closest_ARS.xlsx"
read.xlsx(path_1, header = TRUE, sheetIndex = 1) %>% mutate_if(is.factor, as.character) -> Hawkins_df
# Hawkins_df <- Hawkins_df[0:459,c(0:6)]
Hawkins_df$Chr <- as.factor(Hawkins_df$Chr)
Hawkins_df$Chr <- factor(Hawkins_df$Chr, levels = values)
# Hawkins_df$Chr <- mapvalues(Hawkins_df$Chr, from = levels(Hawkins_df$Chr), to = levels(ARS_df$Chr))

path_2 = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/Query Output/Chr_info_2020-04-06.csv"
read.csv(path_2, header = TRUE) -> CHR_df
CHR_df$Chr <- factor(CHR_df$Chr, levels = values)

path_to_ORFs = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/Query Output/Verified_ORF_info_2020-06-01.xlsx"
read.xlsx(path_to_ORFs, header = TRUE, sheetIndex=1, rowIndex = NULL,endRow=5287) %>% mutate_if(is.factor, as.character) -> ORF_df
ORF_df$ChrNum <- as.factor(ORF_df$ChrNum)
ORF_df <- filter(ORF_df, ORF_df$ChrNum != "chrmt")
names(ORF_df)[names(ORF_df)=="ChrNum"] <- "Chr"
ORF_df <- arrange(ORF_df, Chr, Genome_Start)
ORF_df$Chr <- factor(ORF_df$Chr, levels = values)

path_to_Rhee = 'C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/Rhee_Data/Top_PIC_genes_modified.xls'
read.xlsx(path_to_Rhee, header = TRUE, sheetIndex = 1) %>% mutate_if(is.factor, as.character) -> Rhee_df
Rhee_df <- Rhee_df[1:6045,1:23]
names(Rhee_df)[names(Rhee_df)=="chrom"] <- "Chr"
Rhee_columns <- colnames(Rhee_df)
Rhee_df %>% mutate_at(Rhee_columns[c(2,3,6,8,18,19)], factor) -> Rhee_df
Rhee_df$Chr <- mapvalues(Rhee_df$Chr, from = levels(Rhee_df$Chr), to = levels(Hawkins_df$Chr))
Rhee_df$TSS <- as.integer(Rhee_df$TSS)


#Initialize array with # of elements = to # of ARS and then create a column to put gene hits
Hawkins_array <- array(data = NA, length(Hawkins_df$Chr))
Hawkins_df['Gene_Upstream'] <- Hawkins_array
Hawkins_df['SystematicName_Upstream'] <- Hawkins_array
Hawkins_df['TXN_Direction_Upstream'] <- Hawkins_array
Hawkins_df['Phenotype_Upstream']  <- Hawkins_array
Hawkins_df['Gene_Downstream'] <- Hawkins_array
Hawkins_df['SystematicName_Downstream'] <- Hawkins_array
Hawkins_df['TXN_Direction_Downstream'] <- Hawkins_array
Hawkins_df['Phenotype_Downstream']  <- Hawkins_array
Hawkins_df['Intergenic_Up_Length'] <- Hawkins_array
Hawkins_df['Intergenic_Down_Length'] <- Hawkins_array
Hawkins_df['Intergenic_Length'] <- Hawkins_array
#Create variables to use in for loop (Update: 2020_5_27 Distance updated to 500 after reviewing A genome view of DNA replication paper by macAlpine)
Hawkins_pos_minus1000 = Hawkins_df$Position - 550
Hawkins_pos_plus1000 = Hawkins_df$Position + 550
Genome_St = ORF_df$Genome_Start
Genome_End = ORF_df$Genome_End

#For all of the ARS, determine genes that are in the same chromosome, then compare Start and End positions with ARS gene start and 
#end positions then filter ORF dataframe and select Gene Name and Strand for further quantification.
for (ars in 1:length(Hawkins_array)){
  same_chromosome <- ORF_df$Chr == Hawkins_df$Chr[ars] #determine genes in same chromosome to current ARS
  #logic for determining whether gene is upstream accounting for directionality
  is_upstream_plus = Genome_St < Hawkins_pos_minus1000[ars] & Genome_End > Hawkins_pos_minus1000[ars] & Genome_St < Hawkins_pos_plus1000[ars] & 
    Genome_End < Hawkins_pos_plus1000[ars] 
  is_upstream_minus =  Genome_St > Hawkins_pos_minus1000[ars] & Genome_End <Hawkins_pos_minus1000[ars] & Genome_St <Hawkins_pos_plus1000[ars] & 
    Genome_End < Hawkins_pos_plus1000[ars]
  #filter dataframe to gene that is in front of ARS
  gene_holder_upstream <- ORF_df %>% filter((same_chromosome) & (is_upstream_plus | is_upstream_minus)) %>% dplyr::select(SystematicName,GeneName,Strand,Phenotype)
  #if it is not empty, then add the value to the column of ARS_df
  if( !is_empty(gene_holder_upstream[[1]])){
    Hawkins_df$SystematicName_Upstream[[ars]] <- gene_holder_upstream['SystematicName'][[1]]
    Hawkins_df$Gene_Upstream[[ars]] <- gene_holder_upstream['GeneName'][[1]]
    Hawkins_df$TXN_Direction_Upstream[[ars]] <- gene_holder_upstream['Strand'][[1]]
    Hawkins_df$Phenotype_Upstream[[ars]] <- strsplit(gene_holder_upstream['Phenotype'][[1]], ";")[[1]][1]
  }
  
  #logic for determining whether gene is upstream accounting for directionality
  is_downstream_plus = Genome_St > Hawkins_pos_minus1000[ars] & Genome_End > Hawkins_pos_minus1000[ars] & Genome_St < Hawkins_pos_plus1000[ars] & 
    Genome_End > Hawkins_pos_plus1000[ars]
  is_downstream_minus = Genome_St > Hawkins_pos_minus1000[ars] & Genome_End > Hawkins_pos_minus1000[ars] & Genome_St > Hawkins_pos_plus1000[ars] &
    Genome_End < Hawkins_pos_plus1000[ars]
  #filter dataframe to gene that is after the ARS 
  gene_holder_downstream <- ORF_df %>% filter((same_chromosome) & (is_downstream_plus | is_downstream_minus)) %>% dplyr::select(SystematicName,GeneName,Strand,Phenotype)
  #if it is not empty, then add the value to the column of ARS_df
  if( !is_empty(gene_holder_downstream[[1]])){
    Hawkins_df$SystematicName_Downstream[[ars]] <- gene_holder_downstream['SystematicName'][[1]]
    Hawkins_df$Gene_Downstream[[ars]] <- gene_holder_downstream['GeneName'][[1]]
    Hawkins_df$TXN_Direction_Downstream[[ars]] <- gene_holder_downstream['Strand'][[1]]
    Hawkins_df$Phenotype_Downstream[[ars]] <- strsplit(gene_holder_downstream['Phenotype'][[1]], ";")[[1]][1]
  }
  #Determine in what intergenic region the ars is. 
  intergenic_down <- which(same_chromosome & Hawkins_df$Position[ars] < ORF_df$Downstream_Intergenic_End & Hawkins_df$Position[ars] > ORF_df$Downstream_Intergenic_Start)
  intergenic_up <- which(same_chromosome & Hawkins_df$Position[ars] < ORF_df$Upstream_Intergenic_End & Hawkins_df$Position[ars] > ORF_df$Upstream_Intergenic_Start)
  #If it is between an intergenic region of ORF_df then assign the minimum length (minimum is precaution)
  if (!is_empty(intergenic_down)){
    Hawkins_df$Intergenic_Down_Length[[ars]] <- min(ORF_df[intergenic_down,11] - ORF_df[intergenic_down,10])
  }
  if(!is_empty(intergenic_up)){
    Hawkins_df$Intergenic_Up_Length[[ars]] <- min(ORF_df[intergenic_up,13] - ORF_df[intergenic_up,12])
  }
  #If it is found in an intergenic region, assign minimum value (also precaution). Else assign 0 (ARS is in a gene)
  if(!(is.na(Hawkins_df$Intergenic_Down_Length[[ars]]) & is.na(Hawkins_df$Intergenic_Up_Length[[ars]]))){
    Hawkins_df$Intergenic_Length[[ars]] <- min(Hawkins_df$Intergenic_Down_Length[[ars]], Hawkins_df$Intergenic_Up_Length[[ars]], na.rm = TRUE)
  }else{
    Hawkins_df$Intergenic_Length[[ars]] <- 0
  }
}
rm(ars,Genome_End,Genome_St,intergenic_down,intergenic_up,is_downstream_minus,is_upstream_plus,is_downstream_plus,is_upstream_minus,same_chromosome)

Hawkins_df$TXN_Direction_Upstream[is.na(Hawkins_df$TXN_Direction_Upstream)] <- 2
Hawkins_df$TXN_Direction_Downstream[is.na(Hawkins_df$TXN_Direction_Downstream)] <- 2

Hawkins_df$Gene_Upstream[is.na(Hawkins_df$Gene_Upstream)] <- 2
Hawkins_df$Gene_Downstream[is.na(Hawkins_df$Gene_Downstream)] <- 2

Hawkins_df$Phenotype_Upstream[is.na(Hawkins_df$Phenotype_Upstream)] <- 2
Hawkins_df$Phenotype_Downstream[is.na(Hawkins_df$Phenotype_Downstream)] <- 2

for (i in 1:length(Hawkins_df$Phenotype_Upstream)){
  Hawkins_df$Phenotype_Upstream[[i]] <- strsplit(Hawkins_df$Phenotype_Upstream[i], " ")[[1]][1]
  Hawkins_df$Phenotype_Downstream[[i]] <- strsplit(Hawkins_df$Phenotype_Downstream[i], " ")[[1]][1]
}

Hawkins_df['ARS_TXN_type'] <- Hawkins_array

Hawkins_df$ARS_TXN_type[which(Hawkins_df$TXN_Direction_Upstream == 1 & Hawkins_df$TXN_Direction_Downstream == 1)] <- "Codirectional"
Hawkins_df$ARS_TXN_type[which(Hawkins_df$TXN_Direction_Upstream == -1 & Hawkins_df$TXN_Direction_Downstream == -1)] <- "Codirectional"
Hawkins_df$ARS_TXN_type[which(Hawkins_df$TXN_Direction_Upstream == 1 & Hawkins_df$TXN_Direction_Downstream == -1)] <- "Convergent"
Hawkins_df$ARS_TXN_type[which(Hawkins_df$TXN_Direction_Upstream == -1 & Hawkins_df$TXN_Direction_Downstream == 1)] <-  "Divergent"
Hawkins_df$ARS_TXN_type[which(Hawkins_df$TXN_Direction_Upstream == 1  & Hawkins_df$TXN_Direction_Downstream == 2)] <- "Upstream_Towards"
Hawkins_df$ARS_TXN_type[which(Hawkins_df$TXN_Direction_Upstream == -1 & Hawkins_df$TXN_Direction_Downstream == 2)] <- "Upstream_Away"
Hawkins_df$ARS_TXN_type[which(Hawkins_df$TXN_Direction_Upstream == 2 & Hawkins_df$TXN_Direction_Downstream == -1)] <- "Downstream_Towards"
Hawkins_df$ARS_TXN_type[which(Hawkins_df$TXN_Direction_Upstream == 2 & Hawkins_df$TXN_Direction_Downstream == 1)] <- "Downstream_Away"
Hawkins_df$ARS_TXN_type[which(Hawkins_df$TXN_Direction_Upstream == 2 & Hawkins_df$TXN_Direction_Downstream == 2)] <- "None_Flanking"

Hawkins_df['Dist_to_TSS'] <- Hawkins_array
for (index in 1:length(Hawkins_df$Chr)){
  #Filter df to PIC locations that are on same chromosome and store in df_holder
  df_holder <-  Rhee_df %>% filter(!is.na(TSS) & TSS != 0 & Chr == Hawkins_df$Chr[index]) %>% dplyr::select(TSS) 
  #Determine absolute bp distance between position of TSS and ARS
  dist_holder <- abs(df_holder$TSS - Hawkins_df$Position[index])
  #Choose the minimum value in vector holding distances
  Hawkins_df$Dist_to_TSS[index] <- min(dist_holder)
}
#Plot the count for each Transcription directionality around ARS 
p <- count(Hawkins_df$ARS_TXN_type) %>% ggplot(aes_string(x='x', y='freq', fill = 'x', label = 'freq')) + 
  geom_col() +  geom_text(nudge_y = 3) +
  ggtitle("Count of Transcription Architecture around ARS (1100bp) ") + 
  xlab("Transcription around ARS") + ylab("Count for each Structure") +
  theme_classic() + theme(legend.position = "none")
ggsave("Count of Transcription Architecture around ARS_1100bp.png", plot = p, device = "png")
rm(p)
#Scatter plot for Intergenic Length and ARS Efficiency and Competence 
p1 <- Hawkins_df %>% ggplot(aes(x= Intergenic_Length, y=Efficiency)) + geom_point()
p2 <- Hawkins_df %>% ggplot(aes(x= Intergenic_Length, y=Competence)) + geom_point()
ggsave("Intergenic Length vs Efficiency and Competence.png", plot = ggarrange(p1, p2), device = "png")
rm(p1,p2)

#Scatter plot for Dist to TSS and ARS Efficiency and Competence
p1 <- Hawkins_df %>% ggplot(aes(x= Dist_to_TSS, y=Efficiency)) + geom_point()
p2 <- Hawkins_df %>% ggplot(aes(x= Dist_to_TSS, y=Competence)) + geom_point()
ggsave("Distance to TSS vsEfficiency and Competence.png", plot = ggarrange(p1, p2), device = "png")
rm(p1,p2)

#Scatter plots for Dist to TSS or Intergenic Length vs ARS_Efficiency and Competence faceting by ARS Transcription
Hawkins_df %>% ggplot(aes(x = Dist_to_TSS, y = Efficiency, color = ARS_TXN_type)) + geom_point() + 
  facet_wrap(ARS_TXN_type~., nrow=2) + theme_classic() + theme(legend.position = "none")
Hawkins_df %>% ggplot(aes(x = Dist_to_TSS, y = Competence, color = ARS_TXN_type)) + geom_point() + 
  facet_wrap(ARS_TXN_type~., nrow=2) + theme_classic() + theme(legend.position = "none")
Hawkins_df %>% ggplot(aes(x = Intergenic_Length, y = Efficiency, color = ARS_TXN_type)) + geom_point() + 
  facet_wrap(ARS_TXN_type~., nrow=2) + theme_classic() + theme(legend.position = "none")
Hawkins_df %>% ggplot(aes(x = Intergenic_Length, y = Competence, color = ARS_TXN_type)) + geom_point() + 
  facet_wrap(ARS_TXN_type~., nrow=2) + theme_classic() + theme(legend.position = "none")

#Boxplots for Distance to TSS, Intergenic Length, Efficiency and Competence for Each Transcription Structure 
p <- Hawkins_df %>% ggplot(aes(x = ARS_TXN_type, y = Efficiency, color = ARS_TXN_type)) + 
  geom_boxplot() + geom_jitter(alpha= .4) + theme_classic() + theme(legend.position = "none") +
  ggtitle("Distribution of Efficiency for each Transcription Structure around ARS") + 
  xlab("Transcription around ARS") + ylab("Efficiency") +
  theme_classic() + theme(legend.position = "none")
ggsave("Boxplot_Transcription around ARS vs Efficiency.png", plot = p, device = "png")
rm(p)

p <- Hawkins_df %>% ggplot(aes(x = ARS_TXN_type, y = Competence, color = ARS_TXN_type)) + 
  geom_boxplot() + geom_jitter(alpha= .4)+ theme_classic() + theme(legend.position = "none") +
  ggtitle("Distribution of Competence for each Transcription Structure around ARS") + 
  xlab("Transcription around ARS") + ylab("Competence") +
  theme_classic() + theme(legend.position = "none")
ggsave("Boxplot_Transcription around ARS vs Competence.png", plot = p, device = "png")
rm(p)

p <- Hawkins_df %>% ggplot(aes(x = ARS_TXN_type, y = Intergenic_Length, color = ARS_TXN_type)) + 
  geom_boxplot() + geom_jitter(alpha= .4) + theme_classic() + theme(legend.position = "none") +
  ggtitle("Distribution of Intergenic Length for each Transcription Structure around ARS") + 
  xlab("Transcription around ARS") + ylab("Intergenic Length") +
  theme_classic() + theme(legend.position = "none")
ggsave("Boxplot_Transcription around ARS vs Intergenic Length.png", plot = p, device = "png")
rm(p)

p <- Hawkins_df %>% ggplot(aes(x = ARS_TXN_type, y = Dist_to_TSS, color = ARS_TXN_type)) + 
  geom_boxplot() + geom_jitter(alpha= .4) + theme_classic() + theme(legend.position = "none")
  ggtitle("Distribution of Dist_to_TSS for each Transcription Structure around ARS") + 
  xlab("Transcription around ARS") + ylab("Dist_to_TSS") +
  theme_classic() + theme(legend.position = "none")
ggsave("Boxplot_Transcription around ARS vs Dist_to_TSS.png", plot = p, device = "png")
rm(p)

write.xlsx(Hawkins_df, "Hawkins_df_modified.xlsx", row.names = FALSE)

# Hawkins_df %>% filter(Intergenic_Length != 0) %>% ggplot(aes(x= Intergenic_Length, y=Efficiency)) + geom_point()
# Hawkins_df %>% filter(Intergenic_Length != 0) %>% ggplot(aes(x= Intergenic_Length, y=Competence)) + geom_point()
# 
# Hawkins_df %>% filter(ARS_TXN_type == "Convergent") %>% ggplot(aes(x = Chr, y = ..count..)) + geom_bar()
# Hawkins_df %>% filter(ARS_TXN_type == "Divergent") %>% ggplot(aes(x = Chr, y = ..count..)) + geom_bar()
# 
# Hawkins_df %>% filter(ARS_TXN_type == "Convergent") %>% ggplot(aes(x = Efficiency, y = IIAHawkins)) + geom_point()
# Hawkins_df %>% filter(ARS_TXN_type == "Convergent") %>% ggplot(aes(x = Competence, y = IIAHawkins)) + geom_point()
# Hawkins_df %>% filter(ARS_TXN_type == "Convergent") %>% ggplot(aes(x = Efficiency, y = Pol.IIHawkins)) + geom_point()
# Hawkins_df %>% filter(ARS_TXN_type == "Convergent") %>% ggplot(aes(x = Competence, y = Pol.IIHawkins)) + geom_point()
# 
# Hawkins_df %>% filter(ARS_TXN_type == "Divergent") %>% ggplot(aes(x = IIAHawkins, y = Efficiency)) + geom_point()
# Hawkins_df %>% filter(ARS_TXN_type == "Divergent") %>% ggplot(aes(x = IIAHawkins, y = Competence)) + geom_point()
# Hawkins_df %>% filter(ARS_TXN_type == "Divergent") %>% ggplot(aes(x = Pol.IIHawkins, y = Efficiency)) + geom_point()
# Hawkins_df %>% filter(ARS_TXN_type == "Divergent") %>% ggplot(aes(x = Pol.IIHawkins, y = Competence)) + geom_point()

# Hawkins_df['ARS_Phen_type'] <- Hawkins_array
# 
# for (i in 1:length(Hawkins_df$Chr)){
#   phen_string <- ""
#   
#   if (Hawkins_df$Phenotype_Upstream[i] == 2){
#   phen_string  <- "NA"
#   }else if (Hawkins_df$Phenotype_Upstream == "None"){
#     phen_string <- "No_Test"
#   }else{
#     phen_string <- substring(Hawkins_df$Phenotype_Upstream[i], 1, 4)
#   }
#   if (Hawkins_df$Phenotype_Downstream[i] == "2"){
#     phen_string  <- paste(phen_string, "NA", sep = "&")
#   }else if (Hawkins_df$Phenotype_Downstream == "None"){
#     phen_string <- paste(phen_string, "No_Test", sep = "&")
#   }else{
#     phen_string <- paste(phen_string, substring(Hawkins_df$Phenotype_Downstream[i], 1, 4), sep = "&")
#   }
#   
#   Hawkins_df$ARS_Phen_type[[i]] <- str_remove_all(phen_string, "-")
# }
# 
# 
# PIC_columns = Rhee_columns[9:17]
# for (i in 1:length(PIC_columns)) {
#   Hawkins_df[paste(PIC_columns[i],"Hawkins", sep = "")] <- Hawkins_array
# }
# 
# count_both = 0
# count_up_only = 0
# count_down_only =0
# both_not_in_Rhee_df = 0
# up_not_in_Rhee_df = 0
# down_not_in_Rhee_df = 0
# not_available = 0 
# for (i in 1:length(Hawkins_df$Chr)){
#   sys_up_in_Rhee <- Hawkins_df$SystematicName_Upstream[i] %in% Rhee_df$gene_id 
#   sys_down_in_Rhee <- Hawkins_df$SystematicName_Downstream[i] %in% Rhee_df$gene_id
#   not_na_upstream <- !is.na(Hawkins_df$SystematicName_Upstream[i])
#   not_na_downstream <- !is.na(Hawkins_df$SystematicName_Downstream[i])
#   
#   if (sys_up_in_Rhee & sys_down_in_Rhee){
#     for (j in 1:length(PIC_columns)){
#       Hawkins_df[,paste(PIC_columns[j],"Hawkins", sep = "")][i] <- Rhee_df[,PIC_columns[j]][which(Rhee_df$gene_id == Hawkins_df$SystematicName_Upstream[i])] + Rhee_df[,PIC_columns[j]][which(Rhee_df$gene_id == Hawkins_df$SystematicName_Downstream[i])]
#          
#     }
#     count_both <- count_both + 1 
#   }else if (sys_up_in_Rhee){
#     for (j in 1:length(PIC_columns)){
#       Hawkins_df[,paste(PIC_columns[j],"Hawkins", sep = "")][i] <- Rhee_df[,PIC_columns[j]][which(Rhee_df$gene_id == Hawkins_df$SystematicName_Upstream[i])]
#        
#     }
#     count_up_only <- count_up_only + 1 
#   }else if (sys_down_in_Rhee){
#     for (j in 1:length(PIC_columns)){
#       Hawkins_df[,paste(PIC_columns[j],"Hawkins", sep = "")][i] <- Rhee_df[,PIC_columns[j]][which(Rhee_df$gene_id == Hawkins_df$SystematicName_Downstream[i])]
#         
#     }
#     count_down_only <- count_down_only + 1
#   }else if (not_na_upstream & not_na_downstream){
#     for (j in 1:length(PIC_columns)){
#       Hawkins_df[,paste(PIC_columns[j],"Hawkins", sep = "")][i] <- 0
#       
#     }
#     both_not_in_Rhee_df = both_not_in_Rhee_df + 1
#   }else if (not_na_upstream){
#     for (j in 1:length(PIC_columns)){
#       Hawkins_df[,paste(PIC_columns[j],"Hawkins", sep = "")][i] <- 0
#       
#     }
#     up_not_in_Rhee_df = up_not_in_Rhee_df + 1
#   }else if(not_na_downstream){
#     for (j in 1:length(PIC_columns)){
#       Hawkins_df[,paste(PIC_columns[j],"Hawkins", sep = "")][i] <- 0
#       
#     }
#     down_not_in_Rhee_df = down_not_in_Rhee_df + 1
#   }else{
#     for (j in 1:length(PIC_columns)){
#       Hawkins_df[,paste(PIC_columns[j],"Hawkins", sep = "")][i] <- 0
#       
#     }
#     not_available <- not_available + 1
#   }
# }
# print(both_not_in_Rhee_df + count_both + count_down_only + count_up_only + down_not_in_Rhee_df + up_not_in_Rhee_df + not_available)
# sum(!is.na(Hawkins_df$SystematicName_Upstream)) + sum(!is.na(Hawkins_df$SystematicName_Downstream))
#   
# 
# Hawkins_df$ARS_Phen_type[which(Hawkins_df$ARS_Phen_type == "Non&Gene")] <- "Non&Non"
# Hawkins_df$Phenotype_Downstream[which(Hawkins_df$Phenotype_Downstream == "Gene")] <- "Non-essential"



# Hawkins_columns <- colnames(Hawkins_df)
# Hawkins_df %<>% mutate_at(Hawkins_columns[c(15,16,19,20,21,22)], factor)

# x_list = c(3:8,10)
# y_list = x_list
# 
# plot_density_and_save_nocolor(Hawkins_df, x_list)
# plot_point_and_save_nocolor(Hawkins_df, x_list, y_list)
# 
# plot_density_and_save(Hawkins_df, x_list, facet_index = 1, color_index = 1)
# plot_density_and_save(Hawkins_df,x_list, facet_index = 21, color_index = 21)
# plot_density_and_save(Hawkins_df, x_list, facet_index = 22, color_index = 22)
# 
# plot_point_and_save(Hawkins_df, x_list, y_list, facet_index = 1, color_index = 1)
# plot_point_and_save(Hawkins_df, x_list, y_list, facet_index = 19, color_index = 19)
# plot_point_and_save(Hawkins_df, x_list, y_list, facet_index = 21, color_index = 21)
# 
# x_list = c(1,21,22)
# plot_box_and_save_nocolor(Hawkins_df, x_list, y_list)
# 
# 
# x_list = c(23:32)
# y_list = x_list
# plot_point_and_save_nocolor(Hawkins_df, x_list, y_list)
# 
# x_list = c(23:32)
# y_list = c(3:6)
# plot_point_and_save_nocolor(Hawkins_df, x_list, y_list)
# 
# x_list = c(23:32)
# y_list = x_list
# plot_point_and_save(Hawkins_df, x_list, y_list, facet_index= 21, color_index = 21)
# 

# x_list = c(23:31)
# plot_density_and_save_nocolor(Hawkins_df, x_list)
# plot_density_and_save(Hawkins_df, x_list, facet_index = 21, color_index = 21)
# 
# 
# Hawkins_df %>% filter(ARS_TXN_type != "None_Flanking") %>% ggplot(aes(x = Pol.IIHawkins, color = ARS_TXN_type)) + geom_density(alpha = .4)+ facet_wrap(ARS_TXN_type~., nrow = 2)
# Hawkins_df %>% filter(ARS_TXN_type != "None_Flanking") %>% ggplot(aes(x = IIAHawkins, color = ARS_TXN_type)) + geom_density(alpha = .4)+ facet_wrap(ARS_TXN_type~., nrow = 2)
# Hawkins_df %>% filter(ARS_TXN_type != "None_Flanking") %>% ggplot(aes(x = TBPHawkins, color = ARS_TXN_type)) + geom_density(alpha = .4)+ facet_wrap(ARS_TXN_type~., nrow = 2)

# df_columns = colnames(Hawkins_df)
# x_index = c(23,24,31)
# for (i in 1:length(x_index)){
#   plt <- Hawkins_df %>% filter(ARS_TXN_type != "None_Flanking" & TBPHawkins != 0) %>% ggplot(aes_string(x = df_columns[x_index[i]], color = "ARS_TXN_type")) + geom_density(alpha = .4) + 
#     facet_wrap(as.formula(paste("ARS_TXN_type","~", ".", sep = "")), nrow = 2)
#   name_for_file <- paste(df_columns[x_index[i]], "density", "by","ARS_TXN_type", "filter", "None_FlankingAndNot0","col", "ARS_TXN_type", ".png", sep = "_")
#   ggsave(name_for_file, plot = plt)
# }
# for (i in 1:length(x_index)){
#   plt <- Hawkins_df %>% filter(ARS_TXN_type != "None_Flanking" & TBPHawkins != 0) %>% ggplot(aes_string(x = df_columns[x_index[i]])) + geom_density(alpha = .4) 
#   name_for_file <- paste(df_columns[x_index[i]], "density", "by","ARS_TXN_type", "filter", "None_FlankingAndNot0","col", "ARS_TXN_type", ".png", sep = "_")
#   ggsave(name_for_file, plot = plt)
# }
# p2 <-  Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = Competence, fill = "red")) + geom_density(alpha = .4)
# p2
# ggsave("Origin Competence Distribution.png")
# p1 <- Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = Competence, fill = Chr)) + geom_density(alpha = .4) + facet_wrap(Chr~., nrow = 2)
# p1
# ggsave("Origin Competence Distribution per Chromosome.png")
# 
# p3 <- Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = Efficiency, fill = "red")) + geom_density(alpha = .4)
# p3
# ggsave("Origin Efficiency Distribution.png")
# 
# p4 <- Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = Efficiency, fill = Chr)) + geom_density(alpha = .4) + facet_wrap(Chr~., nrow = 2)
# p4
# ggsave("Origin Efficiency Distribution per Chromosome.png")
# 
# p5 <- Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = T_wdith, fill = "red")) + geom_density(alpha = .4)
# p5
# ggsave("Origin Time Width Distribution.png")
# 
# p6 <- Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = T_wdith, fill = Chr)) + geom_density(alpha = .4) + facet_wrap(Chr~., nrow = 2)
# p6
# ggsave("Origin Time Width Distribution per Chromosome.png")
# 
# p7 <- Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = T_half, fill = "red")) + geom_density(alpha = .4)
# p7
# ggsave("Origin Time Half Distribution.png")
# 
# p8 <- Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = T_half, fill = Chr)) + geom_density(alpha = .4) + facet_wrap(Chr~., nrow = 2)
# p8
# ggsave("Origin Time Half Distribution per Chromosome.png")

# Hawkins_cols <- colnames(Hawkins_df)
# summary_Hawkins <- Hawkins_df %>% summarize_at(vars(Hawkins_cols[c(3,4,5,6,7,8)]), funs(mean, median, sd, cv))
# summary_Hawkins_by_Chr <- Hawkins_df %>% group_by(Chr) %>% summarize_at(vars(Hawkins_cols[c(3,4,5,6,7,8)]), funs(mean, median, sd, cv)) 
# 
# summary_Hawkins_by_Chr['numARS'] <- count(Hawkins_df$Chr)[2]
# summary_Hawkins_by_Chr['Length'] <- CHR_df$Length

# p1 <- summary_Hawkins_by_Chr %>% ggplot(aes(Length, numARS)) + geom_point() + geom_smooth(method = 'lm')
# p2 <-- summary_Hawkins_by_Chr %>% ggplot(aes(Length, Competence_median)) + geom_point() + geom_smooth(method = 'lm')
# p3 <- summary_Hawkins_by_Chr %>% ggplot(aes(Length, Efficiency_median)) + geom_point() + geom_smooth(method = 'lm')
# p4 <- summary_Hawkins_by_Chr %>% ggplot(aes(Length, OriginDist_median)) + geom_point() + geom_smooth(method = 'lm')
# ggarrange(p1, p2, p3, p4)
# ggsave("Different metrics vs Length of Chromosome.png")

# setwd("C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/Hawkins_Data/")
# write.xlsx(as.data.frame(summary_Hawkins_by_Chr), "Hawkins_Origin_Metrics_Summary_By_Chr.xlsx", row.names = FALSE)
 # mean_by_chr <- Hawkins_df %>% select(Hawkins_cols[c(3,4,5,6,7,8)]) %>% aggregate(by = list(Hawkins_df$Chr), FUN = mean) 
# median_by_chr <- Hawkins_df %>% select(Hawkins_cols[c(3,4,5,6,7,8)]) %>% aggregate(by = list(Hawkins_df$Chr), FUN = median)
# sd_by_chr <- Hawkins_df %>% select(Hawkins_cols[c(3,4,5,6,7,8)]) %>% aggregate(by = list(Hawkins_df$Chr), FUN = sd)


# strsplit(substring(Hawkins_df$Phenotype_Downstream[4], 1,14), " ")[[1]]
# Rhee_df$Standard_Name[Rhee_df$Standard_Name == "None"] <- "Not Available"
# Rhee_df$Standard_Name[is.na(Rhee_df$Standard_Name)] <- "Not Available"

# first_index_up <- 1
# last_index_up <- 459
# first_index_down <- 1
# last_index_down <- 459
# 
# for (i in length(1:4)){
#   first_idx_up = first_index_up + (459*i-1)
#   last_idx_up = last_index_up + (459*i-1)
#   
#   for (j in length(1:4)){
#     first_idx_down = first_index_down + (459*(i-1))
#     last_idx_down = last_index_down + (459*(i-1))
#     
#     which(phenotype_upstream_values[first_idx_up:last_idx_up] & phenotype_downstream_values[first_idx_down:last_idx_down])
#     Hawkins_df$ARS_Phen_type[which(phenotype_upstream_values[first_index_up])] <- phenotype_categories[i+4*(i-1) + j]
#   }
# }

# upstream_essential = grepl("Essential", Hawkins_df$Phenotype_Upstream, fixed = TRUE)
# upstream_non_essential = grepl("Non-essential", Hawkins_df$Phenotype_Upstream, fixed = TRUE)
# upstream_nothing = grepl("2", Hawkins_df$Phenotype_Upstream, fixed = TRUE)
# upstream_not_available = grepl("None", Hawkins_df$Phenotype_Upstream, fixed = TRUE)
# 
# downstream_essential = grepl("Essential", Hawkins_df$Phenotype_Downstream, fixed = TRUE)
# downstream_non_essential = grepl("Non-essential", Hawkins_df$Phenotype_Downstream, fixed = TRUE)
# downstream_nothing = grepl("2", Hawkins_df$Phenotype_Downstream, fixed = TRUE)
# downstream_not_available = grepl("None", Hawkins_df$Phenotype_Downstream, fixed = TRUE)
# 
# phenotype_categories <- c("Ess&Ess","Ess&NonEss","Ess&NotTested","Ess&NA"
#                           ,"NonEss&Ess","NonEss&NonEss","NonEss&NotTested","NonEss&NA"
#                           ,"NotTested&Ess","NotTested&NonEss","NotTested&NotTested","NotTested&NA"
#                           ,"NA&Ess","NA&NonEss","NA&NotTested","NA&NA")
# 
# 
# 
# Hawkins_df$ARS_Phen_type[which(upstream_essential & downstream_essential)] <- phenotype_categories[1]
# Hawkins_df$ARS_Phen_type[which(upstream_essential & downstream_non_essential)]<- phenotype_categories[2]
# Hawkins_df$ARS_Phen_type[which(upstream_essential & downstream_nothing)]<- phenotype_categories[3]
# Hawkins_df$ARS_Phen_type[which(upstream_essential & downstream_not_available)]<- phenotype_categories[4]
# Hawkins_df$ARS_Phen_type[which(upstream_non_essential & downstream_essential)]<- phenotype_categories[5]
# Hawkins_df$ARS_Phen_type[which(upstream_non_essential & downstream_non_essential)]<- phenotype_categories[6]
# Hawkins_df$ARS_Phen_type[which(upstream_non_essential & downstream_nothing)]<- phenotype_categories[7]
# Hawkins_df$ARS_Phen_type[which(upstream_non_essential & downstream_not_available)]<- phenotype_categories[8]
# Hawkins_df$ARS_Phen_type[which(upstream_nothing & downstream_essential)]<- phenotype_categories[9]
# Hawkins_df$ARS_Phen_type[which(upstream_nothing & downstream_non_essential)]<- phenotype_categories[10]
# Hawkins_df$ARS_Phen_type[which(upstream_nothing & downstream_nothing)]<- phenotype_categories[11]
# Hawkins_df$ARS_Phen_type[which(upstream_nothing & downstream_not_available)]<- phenotype_categories[12]
# Hawkins_df$ARS_Phen_type[which(upstream_not_available & downstream_essential)]<- phenotype_categories[13]
# Hawkins_df$ARS_Phen_type[which(upstream_not_available & downstream_non_essential)]<- phenotype_categories[14]
# Hawkins_df$ARS_Phen_type[which(upstream_not_available & downstream_nothing)]<- phenotype_categories[15]
# Hawkins_df$ARS_Phen_type[which(upstream_not_available & downstream_not_available)]<- phenotype_categories[16]
