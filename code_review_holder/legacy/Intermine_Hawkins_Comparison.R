library(tidyverse)
library(ggplot2)
library(xlsx)
library(gridExtra)
library(ggpubr)
library(plyr)

#Need to fix mapping de chromosomes to proper factor (9 esta incorrecto porque cuando haces sorting, es #6 en romano. Cambiar a #)
values <- c("chrI",    "chrII",  "chrIII",  "chrIV",   "chrV",    "chrVI",   "chrVII",  "chrVIII", "chrIX",   "chrX",   
           "chrXI",   "chrXII",  "chrXIII", "chrXIV",  "chrXV",   "chrXVI")

path_1 = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/Query Output/ARS_info_2020_4_26.csv"
read.csv(path_1, header = TRUE) %>% mutate_if(is.factor, as.character) -> ARS_df
ARS_df$Chr <- as.factor(ARS_df$Chr)
ARS_df$Chr <- factor(ARS_df$Chr, levels = values)

path_2 = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/Hawkins_Data/Origin_Data_Hawkins_2013.xlsx"
read.xlsx(path_2, header = TRUE, sheetIndex = 1) %>% mutate_if(is.factor, as.character) -> Hawkins_df
Hawkins_df <- Hawkins_df[0:459,c(0:6)]
colnames(Hawkins_df) <- c("Chr", "Position", "T_half", "T_wdith", "Competence", "Efficiency")
Hawkins_df$Chr <- as.factor(Hawkins_df$Chr)
Hawkins_df$Chr <- mapvalues(Hawkins_df$Chr, from = levels(Hawkins_df$Chr), to = levels(ARS_df$Chr))


path_3 = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/Query Output/Chr_info_2020-04-06.csv"
read.csv(path_3, header = TRUE) -> CHR_df
CHR_df$Chr <- factor(CHR_df$Chr, levels = values)


ARS_array <- array(data = NA, length(ARS_df$Chr))
ARS_df['OriginDist'] <- ARS_array
for (index in 1:length(ARS_df$Chr)){
  
  if (index == length(ARS_df$Chr)){
    ARS_df$OriginDist[index] <- CHR_df$Length[which(CHR_df$Chr == ARS_df$Chr[index])]-ARS_df$Start[index]
  } 
  else{
    if (ARS_df$Chr[index] == ARS_df$Chr[index+1]){
      ARS_df$OriginDist[index] <- ARS_df$Start[index+1] - ARS_df$Start[index]
    } 
    else{
      ARS_df$OriginDist[index] <- CHR_df$Length[which(CHR_df$Chr == ARS_df$Chr[index])]-ARS_df$Start[index]
    }
  }
}

ARS_array <- array(data = NA, length(Hawkins_df$Chr))
Hawkins_df['OriginDist'] <- ARS_array
for (index in 1:length(Hawkins_df$Chr)){
  
  if (index == length(Hawkins_df$Chr)){
    Hawkins_df$OriginDist[index] <- CHR_df$Length[which(CHR_df$Chr == Hawkins_df$Chr[index])]-Hawkins_df$Position[index]
  } 
  else{
    if (Hawkins_df$Chr[index] == Hawkins_df$Chr[index+1]){
      Hawkins_df$OriginDist[index] <- Hawkins_df$Position[index+1] - Hawkins_df$Position[index]
    } 
    else{
      Hawkins_df$OriginDist[index] <- CHR_df$Length[which(CHR_df$Chr == Hawkins_df$Chr[index])]-Hawkins_df$Position[index]
    }
  }
}
#Distribution of Origin Distances by Chromosome for Intermine Dataframe and Hawkins 2013 Dataframe

# p1 <- ARS_df %>% group_by(Chr) %>% ggplot(aes(x = OriginDist, fill = Chr)) + geom_density(alpha = .4) + facet_wrap(Chr~., nrow = 2)
# p1
# ggsave("OriginDist_ARS_df_perChr_faceted.png")
# p2 <- Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = OriginDist, fill = Chr)) + geom_density(alpha = .4) + facet_wrap(Chr~., nrow = 2)
# p2
# ggsave("OriginDist_Hawkins_df_perChr_faceted.png")


#Numnber of Distances over 50kb in each chromosome for Both Dataframes
# num_over_50kb_ARS <- ARS_df %>% filter(OriginDist > 50000) %>% group_by(Chr) %>% nrow() %>% as.character()
# p3 <- ARS_df %>% filter(OriginDist > 50000) %>% group_by(Chr) %>% ggplot(aes(x = Chr, y = ..count..)) + geom_bar()

# num_over_50kb_Hawkins <- Hawkins_df %>% filter(OriginDist >50000) %>% group_by(Chr) %>% nrow() %>% as.character()
# p4 <- Hawkins_df %>% filter(OriginDist >50000) %>% group_by(Chr) %>% ggplot(aes(x = Chr, y =  ..count..)) + geom_bar()

# ggarrange(p3, p4, ncol = 2, labels = c(paste(num_over_50kb_ARS, '(# OriginDist over 50kb)', sep = ""), paste(num_over_50kb_Hawkins, '(# OriginDist over 50kb)', sep = "")))
# ggsave("Intermine_vs_Hawkins_OriginDist_comp.png")

#Determine OriginDist distribution normalizd to length of chromosome

Hawkins_df['OriginDist_to_Chr_Length'] <- ARS_array

for (index in 1:length(Hawkins_df$Chr)) {
  Hawkins_df$OriginDist_to_Chr_Length[index] <- Hawkins_df$OriginDist[index] / CHR_df$Length[which(CHR_df$Chr == Hawkins_df$Chr[index])]
}

p5 <- Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = OriginDist_to_Chr_Length, fill = Chr)) + geom_density(alpha = .4) + facet_wrap(Chr~., nrow = 2)
p5
#Determining closest intermine ARS to Hawkins origin position
ARS_array <- array(data = NA, length(Hawkins_df$Chr))
Hawkins_df['ClosestARS'] <- ARS_array
Hawkins_df['ClosestARS_Start'] <- ARS_array
Hawkins_df['ClosestARS_End'] <-ARS_array
Hawkins_df['Distance_to_closest_ARS'] <- ARS_array

for (index in 1:length(Hawkins_df$Chr)){
  
  df_holder <- ARS_df %>% filter(Chr == Hawkins_df$Chr[index])
  start_holder <- abs(df_holder$Start - Hawkins_df$Position[index])
  end_holder <- abs(df_holder$End - Hawkins_df$Position[index])
  Hawkins_df$Distance_to_closest_ARS[index] <- min(c(min(start_holder), min(end_holder)))
  
  if (which.min(start_holder) == which.min(end_holder)){
    Hawkins_df$ClosestARS[index] <- df_holder$Standard_Name[which.min(start_holder)]
    Hawkins_df$ClosestARS_Start[index] <- df_holder$Start[which.min(start_holder)]
    Hawkins_df$ClosestARS_End[index] <- df_holder$End[which.min(start_holder)]
  } else{
    
    if(which.min(c(start_holder[which.min(start_holder)], end_holder[which.min(end_holder)])) == 1){
      
      Hawkins_df$ClosestARS[index] <- df_holder$Standard_Name[which.min(start_holder)]
      Hawkins_df$ClosestARS_Start[index] <- df_holder$Start[which.min(start_holder)]
      Hawkins_df$ClosestARS_End[index] <- df_holder$End[which.min(start_holder)]
      
    } else{
      
      Hawkins_df$ClosestARS[index] <- df_holder$Standard_Name[which.min(end_holder)]
      Hawkins_df$ClosestARS_Start[index] <- df_holder$Start[which.min(end_holder)]
      Hawkins_df$ClosestARS_End[index] <- df_holder$End[which.min(end_holder)]
      }
    }
    
}
p6 <- Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = Distance_to_closest_ARS)) + geom_density(fill ="red") +theme_classic()
p6
p7 <- Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = Distance_to_closest_ARS, fill = Chr)) + geom_density(alpha = .4) + facet_wrap(Chr~., nrow = 2) 
p7
p8 <- Hawkins_df %>% group_by(Chr) %>% ggplot(aes(x = Distance_to_closest_ARS)) + geom_histogram(fill ="red") +theme_classic()
p8

mean(Hawkins_df$Distance_to_closest_ARS)
median(Hawkins_df$Distance_to_closest_ARS)
sd(Hawkins_df$Distance_to_closest_ARS)

setwd("C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/Hawkins_Data/")
write.xlsx(Hawkins_df, "Hawkins_df_with_closest_ARS.xlsx", row.names = FALSE)
#No funciono. 
# ARS_df <- ARS_df %>% group_by(Chr) %>% mutate(OriginDist = ifelse(Chr == Chr, diff(Start),  CHR_df$Length[which(CHR_df == Chr)]-Start))
# dist <- ARS_df %>% group_by(Chr) %>% filter(OriginDist < 0) %>% mutate(dist = CHR_df$Length[which(CHR_df == as.character(Chr))]-Start) %>% select(dist)
# ARS_df$OriginDist[which(ARS_df$OriginDist < 0)] <- dist$dist
# ARS_df$OriginDist[352] <- CHR_df$Length[which(CHR_df == as.character(ARS_df$Chr[352]))]-ARS_df$Start[352]

# Hawkins_df <- Hawkins_df %>% group_by(Chr) %>% mutate(OriginDist = ifelse(Chr == Chr, diff(Position),  CHR_df$Length[which(CHR_df == Chr)]-Position))
# dist <- Hawkins_df %>% filter(OriginDist < 0) %>% mutate(dist = CHR_df$Length[which(CHR_df == as.character(Chr))]-Position) %>% select(dist)
# Hawkins_df$OriginDist[which(Hawkins_df$OriginDist < 0)] <- dist$dist
# Hawkins_df$OriginDist[459] <- CHR_df$Length[which(CHR_df == as.character(Hawkins_df$Chr[459]))]-Hawkins_df$Position[459]
# 
# diff(ARS_df$Chr)
