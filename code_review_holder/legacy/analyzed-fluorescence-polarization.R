# Install the 'readxl' package if it's not already installed
# install.packages("readxl")

# Load the 'readxl' package
library(readxl)
library(tidyverse)
library(broom)
library(ggplot2)
library(data.table)
library(stats)
file_path <- "C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Experiments\\ATP hydrolysis by Orc1_Orc4\\Assays\\DNA binding\\2023_05_11 Fluorescence polarization\\Data\\2023_05_17 250 KGlut 20mM KCl WT ORC Annie vs Anusha DNA.txt"

# Set the file path of the Excel file

# Read the Excel file into a dataframe
df <- read.delim(file_path, header = TRUE, sep = "\t", skip = 1)
 
df <- df %>%
  .[,3:ncol(.)] %>%
  .[, colSums(is.na(.)) < nrow(.)] %>%
  .[rowSums(is.na(.)) < ncol(.),]

# as.character(as.numeric(rownames(df[duplicated(df),])):rownames(df[nrow(df),]))
df <- df[as.character(as.numeric(rownames(df[duplicated(df),])):rownames(df[nrow(df),])),] 

names(df) <- seq_along(1:ncol(df))
df <- df[-1,]
df <- df[,seq_len(ncol(df)) %% 2 == 0]

df_samples <- read_xlsx("C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Experiments\\ATP hydrolysis by Orc1_Orc4\\Assays\\DNA binding\\2023_05_11 Fluorescence polarization\\23_05_17 250 KGlut WTORC.xlsx", sheet = "Wells")
Concentrations <- as.numeric(unlist(unique(df_samples[,"Concentration"])))
DNA_length <- as.numeric(unlist(unique(df_samples[,"DNA"])))
Sample_labels <- c()
repeat_label <- seq(from=1, to=3, by=1)
for (i in 1:length(DNA_length)){
  experiment_label <- unlist(rep(DNA_length[i], each = 3))
  Sample_labels <- c(Sample_labels, paste(experiment_label, repeat_label, sep="."))
  # Sample_labels <- c(Sample_labels, unlist(rep(DNA_length[i], each = 3)))
}

df <- rbind(df, Concentrations)
df <- t(df)
df <- as.data.frame(df)

colnames(df) <- as.character(c(Sample_labels, "Concentrations"))


df_substracted <- mapply(function(exp_repeat){
  exp_repeat - as.numeric(tail(exp_repeat, 1))
  
}, exp_repeat = df[colnames(df) != "Concentrations"]) 

df_substracted<-as.data.frame(cbind(df_substracted, Concentrations))[1:dim(df_substracted)[1]-1,]

#https://datacornering.com/dplyr-error-in-select-unused-argument/
df_substracted %>%
  dplyr::select(contains("50")) %>%
  rowMeans()
  
apply(df_substracted, 1, sd)  
  dplyr::select(contains("50"))
  
apply(df_substracted, 1, stderr, 3)


DF %>%
  transmute(ID,
            Mean = rowMeans(across(C1:C3)))
# df_substracted %>%
#   ggplot(aes(x=.[,3], y = Concentrations)) +
#   geom_po

stderr <- function(x, n) {
    sd(x)/sqrt(n)
  }

#rownames(df) <- c(as.character(Sample_labels), "Concentration")

# for (i in 1:nrow(df)){
#   if (i >= 3) {
#     rownames(df)[i] <- paste(as.character(Sample_labels)[i], as.character(i %% 3), sep = ".")
#   } else {
#     rownames(df)[i] <- paste(as.character(Sample_labels)[i], as.character(i), sep = ".")
#   }
# }