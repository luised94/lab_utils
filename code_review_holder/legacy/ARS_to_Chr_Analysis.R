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


CHR_df <- CHR_df %>% mutate(Proportion_Length = (Length)/sum(Length))
total_ARS <- length(ARS_df$Chr)
total_ARS_Hawkins <- length(Hawkins_df$Chr)
#Determining proportions of ARS to proportion of total length for yeast chromosomes
CHR_df['numARS'] <- ARS_df %>% select(Chr) %>% group_by(Chr) %>% count() %>% pull(freq)
CHR_df <- CHR_df %>% mutate(ARS_to_Chr_Length = (numARS/Length)*1000)
CHR_df <- CHR_df %>% mutate(Proportion_of_ARS = (numARS/total_ARS))
CHR_df <- CHR_df %>% mutate(Ratio_ARS_to_Length = (Proportion_of_ARS/Proportion_Length))
CHR_df <- CHR_df %>% mutate(Kb_per_ARS = 1/Proportion_of_ARS)

CHR_df['numARS_Hawkins'] <- Hawkins_df %>% select(Chr) %>% group_by(Chr) %>% count() %>% pull(freq)
CHR_df <- CHR_df %>% mutate(Hawkins_ARS_to_Chr_Length = (numARS_Hawkins/Length)*1000)
CHR_df <- CHR_df %>% mutate(Proportion_Hawkins = (numARS_Hawkins/total_ARS_Hawkins))
CHR_df <- CHR_df %>% mutate(Ratio_Hawkins = (Proportion_Hawkins/Proportion_Length))
CHR_df <- CHR_df %>% mutate(Kb_per_ARS_Hawkins = 1/Proportion_Hawkins)

#Calculating differences between intermine database and Hawkins 2013 study 
CHR_df <- CHR_df %>% mutate(Diff = numARS - numARS_Hawkins)
CHR_df <- CHR_df %>% mutate(highlight = ifelse(Chr == 'chrXII', T, F))

proportion_of_length_plot <- CHR_df %>% ggplot(aes_string("Chr", colnames(CHR_df[colnames(CHR_df[3])]))) + geom_col() + theme_classic()

numARS_plot <- CHR_df %>% ggplot(aes_string("Chr", colnames(CHR_df[colnames(CHR_df[4])]))) + geom_col() + theme_classic()
p1 <- CHR_df %>% ggplot(aes_string("Chr", colnames(CHR_df[colnames(CHR_df[5])]))) + geom_col(aes(fill = highlight)) + theme_classic()
p2 <- CHR_df %>% ggplot(aes_string("Chr", colnames(CHR_df[colnames(CHR_df[6])]))) + geom_col(aes(fill = highlight)) + theme_classic()
p3 <- CHR_df %>% ggplot(aes_string("Chr", colnames(CHR_df[colnames(CHR_df[7])]))) + geom_col(aes(fill = highlight)) + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme_classic()
p4 <- CHR_df %>% ggplot(aes_string("Chr", colnames(CHR_df[colnames(CHR_df[8])]))) + geom_col(aes(fill = highlight)) + theme_classic()


ggarrange(p2,                                                 # First row with scatter plot
          ggarrange(p3, p4, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
          nrow = 2, 
          labels = "A"                                        # Labels of the scatter plot
) 
ggsave("Intermine_ARS_Chromosome_Length_Comparison.png")

p5 <- CHR_df %>% ggplot(aes_string("Chr", colnames(CHR_df[colnames(CHR_df[11])]))) + geom_col(aes(fill = highlight)) + theme_classic()
p6 <- CHR_df %>% ggplot(aes_string("Chr", colnames(CHR_df[colnames(CHR_df[12])]))) + geom_col(aes(fill = highlight)) + geom_hline(yintercept=1, linetype="dashed", color = "red") + theme_classic()
p7 <- CHR_df %>% ggplot(aes_string("Chr", colnames(CHR_df[colnames(CHR_df[13])]))) + geom_col(aes(fill = highlight)) + theme_classic()

ggarrange(p5,                                                 # First row with scatter plot
          ggarrange(p6, p7, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
          nrow = 2, 
          labels = "A"                                        # Labels of the scatter plot
) 
ggsave("Hawkins_ARS_Chromosome_Length_Comparison.png")

p8 <- CHR_df %>% ggplot(aes(Chr, Diff)) + geom_col(aes(fill = highlight)) + theme_classic()
p8
ggsave("Diff_#ARS_Intermine_Hawkins.png")
# indexes <- c(2,3)
# plot_vector <- list()
# for (index in indexes){
#   plot_holder <- CHR_df %>% ggplot(aes_string("Chr", colnames(CHR_df[colnames(CHR_df[index])]))) + geom_col() + theme_classic()
#   
#   }

write.xlsx(CHR_df, "ARS_Chromosome_Analysis.xlsx", row.names = FALSE)





# sum(CHR_df$numARS)
# sum(CHR_df$numARS_Hawkins)
# sum(CHR_df$Proportion_of_ARS)
# sum(CHR_df$Proportion_Hawkins)
# # ARS_df %>% ggplot(aes(Chr, y = ..count..)) + geom_bar() + theme_classic() + labs(title = "#ARS per Chromosome", x = "Chrom Num", y = "# of ARS")
# arsperChr %>% ggplot(aes(Chr, y = Diff)) + geom_col() + theme_classic()
# ARS_per_Kb <- CHR_df %>% ggplot(aes(Chr, numARS)) + geom_col() + theme_classic()
# ARS_per_Kb
# 
# ggarrange(bxp, dp, bp + rremove("x.text"), 
#           labels = c("A", "B", "C"),
#           ncol = 2, nrow = 2)

