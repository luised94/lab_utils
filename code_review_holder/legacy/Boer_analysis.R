library(readxl)
library(ggplot2)
library(tidyverse)
library(plyr)
path_to_Boer = "C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/Boer_Data/TF_Parameters_Boer_2019.xlsx"
Boer_df <- as.data.frame(read_excel(path_to_Boer, sheet = "positional params - basic",skip = 1)) 
Boer_columns <- colnames(Boer_df)
Boer_df %>% mutate_at(Boer_columns[c(3,6,7)], as.numeric) -> Boer_df


Boer_df %>% filter(activity < 0 & TFAlias == "ORC1")
max(Boer_df$activity, na.rm=TRUE); min(Boer_df$activity, na.rm = TRUE)
