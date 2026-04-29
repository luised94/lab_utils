library(ggplot2)
library(xlsx)
library(tidyverse)
#ATPase data must have same format as the Orc5 and 6 data.
create_experiment_label <- function(number) {
  label <- paste("Exp", number)
  return(label)
}
process_atpase_data <- function(file_path, sheet_index, Timepoints, sample_names, experiment_number, indices) {
  library(xlsx)
  library(tidyverse)
  
  atpase_data <- read.xlsx(file_path, header = FALSE, sheetIndex = sheet_index)
  atpase_data <- as.data.frame(atpase_data)
  colnames(atpase_data) <- c("Idx", "ADP")
  
  atpase_data <- mutate(atpase_data, Even_Odds = Idx %% 2)
  atpase_data$ATP <- NA
  
  for (i in 1:length(atpase_data$Idx)) {
    if (!atpase_data$Even_Odds[i]) {
      atpase_data$ATP[i] <- atpase_data$ADP[i - 1]
    }
  }
  
  atpase_data <- na.omit(atpase_data)
  atpase_data <- atpase_data[, grep("ADP|ATP", colnames(atpase_data))]
  
  Time_points <- Timepoints
  atpase_data$Timepoint <- NA
  atpase_data[1, 'Timepoint'] <- 90
  
  for (i in 2:length(atpase_data$ADP)) {
    if (i %% 4 == 2) {
      atpase_data$Timepoint[i] <- Time_points[1]
    } else if (i %% 4 == 3) {
      atpase_data$Timepoint[i] <- Time_points[2]
    } else if (i %% 4 == 0) {
      atpase_data$Timepoint[i] <- Time_points[3]
    } else if (i %% 4 == 1) {
      atpase_data$Timepoint[i] <- Time_points[4]
    }
  }
  
  atpase_data$Sample <- NA
  atpase_data[1, 'Sample'] <- "No_ORC"
  
  indices <- indices
  tracker <- 1
  for (i in indices) {
    if (indices[tracker] < indices[length(indices)]) {
      while (i < indices[tracker + 1]) {
        atpase_data$Sample[i] <- sample_names[tracker]
        i <- i + 1
      }
      tracker <- tracker + 1
    }
  }
  
  atpase_data$Sample <- as.factor(atpase_data$Sample)
  atpase_data$Experiment <- as.factor(create_experiment_label(experiment_number))
  
  atpase_data <- mutate(atpase_data, Percent = ADP / (ADP + ATP))
  atpase_data <- mutate(atpase_data, Percent_corr = Percent-Percent[1]) #Background subtraction
  return(atpase_data)
}

####First Orc6 EK experiment ####
#Data from ImageJ analysis is stored in sheets 6-11.
path_1 = "C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Experiments\\ATP hydrolysis by Orc1_Orc4\\Assays\\ATPase\\2021_07_29 ORC4R O5 and O6 EK\\2022_03_01 ORC4R O5,O6 EK.xlsx"
sheet_index = 6
Time_points = c(0, 15, 45, 90)
Sample_Names <- c("WT","4R", "Orc6_E304K_4R")
Experiment_number = 1
Indices = c(2, 6, 10, 14)
atpase_data <- process_atpase_data(path_1, 
                                   sheet_index = sheet_index, 
                                   Timepoints = Time_points, 
                                   sample_names = Sample_Names, 
                                   experiment_number = Experiment_number,
                                   indices = Indices)
all_experiments_df <- data.frame()
all_experiments_df <- rbind(all_experiments_df, atpase_data)
# lapply(1:3, function(x){
#   x
#   Sample_Names
# })
# atpase_data = read.xlsx(path_1, header = FALSE, sheetIndex = 6)

####Second Orc6 EK Experiment ####
sheet_index = 7
Experiment_number = 2

# atpase_data = read.xlsx(path_1, header = FALSE, sheetIndex = 7) 

atpase_data <- process_atpase_data(path_1, 
                                   sheet_index = sheet_index, 
                                   Timepoints = Time_points, 
                                   sample_names = Sample_Names, 
                                   experiment_number = Experiment_number,
                                   indices = Indices)
all_experiments_df <- rbind(all_experiments_df, atpase_data)
####Third Orc6 EK Experiment####
sheet_index = 8
Experiment_number = 3

# atpase_data = read.xlsx(path_1, header = FALSE, sheetIndex = 8) 

atpase_data <- process_atpase_data(path_1, 
                                   sheet_index = sheet_index, 
                                   Timepoints = Time_points, 
                                   sample_names = Sample_Names, 
                                   experiment_number = Experiment_number,
                                   indices = Indices)

all_experiments_df <- rbind(all_experiments_df, atpase_data)

####First Orc5 EK Experiment####
sheet_index = 9
Experiment_number = 4
Sample_Names <- c("WT","4R", "Orc5_E104K_4R")
# atpase_data = read.xlsx(path_1, header = FALSE, sheetIndex = sheet_index) 

atpase_data <- process_atpase_data(path_1, 
                                   sheet_index = sheet_index, 
                                   Timepoints = Time_points, 
                                   sample_names = Sample_Names, 
                                   experiment_number = Experiment_number,
                                   indices = Indices)
all_experiments_df <- rbind(all_experiments_df, atpase_data)
####Second Orc5 EK Experiment####
sheet_index = 10
Experiment_number = 5
atpase_data = read.xlsx(path_1, header = FALSE, sheetIndex = sheet_index) 

atpase_data <- process_atpase_data(path_1, 
                                   sheet_index = sheet_index, 
                                   Timepoints = Time_points, 
                                   sample_names = Sample_Names, 
                                   experiment_number = Experiment_number,
                                   indices = Indices)
all_experiments_df <- rbind(all_experiments_df, atpase_data)
####Third Orc5 EK Experiment####
sheet_index = 11
Experiment_number = 6
atpase_data = read.xlsx(path_1, header = FALSE, sheetIndex = sheet_index) 

atpase_data <- process_atpase_data(path_1, 
                                   sheet_index = sheet_index, 
                                   Timepoints = Time_points, 
                                   sample_names = Sample_Names, 
                                   experiment_number = Experiment_number,
                                   indices = Indices)
all_experiments_df <- rbind(all_experiments_df, atpase_data)
 
####Experiment Analysis####
number_of_samples = count(all_experiments_df, Sample, Timepoint)
# all_experiments_df <- all_experiments_df %>% mutate(ADP_per_ORC = ADP_pmols/pmol_ORC)
all_experiments_df['Sample_Size'] <-array(data = NA, length(all_experiments_df$ADP))
for (i in 1:length(all_experiments_df$ADP)){
  all_experiments_df$Sample_Size[i] <- number_of_samples[number_of_samples$Sample == all_experiments_df$Sample[i], 3][1] 
}
all_experiments_df2 <- all_experiments_df
experiments456 <- lapply(c("Exp 4", "Exp 5", "Exp 6"), function(x){
  grepl(x, all_experiments_df$Experiment)
})
experiments456 <- Reduce(`|`, lapply(experiments456, unlist))

filtered_df <- all_experiments_df %>% filter(!(Sample == "WT" & experiments456) & !(Sample == "Orc5_E104K_4R" & Experiment == "Exp 4"))
# 
# 
summary_df <- filtered_df %>% group_by(Timepoint, Sample) %>% summarise(avg = mean(Percent_corr), std = sd(Percent_corr))
# summary_df <- all_experiments_df %>% group_by(Timepoint, Sample) %>% summarise(avg = mean(ADP_per_ORC), std = sd(ADP_per_ORC))
# summary_df['SEM'] <-array(data = NA, length(summary_df$avg))
# for (i in 1:length(summary_df$SEM)){
  # summary_df$SEM[i] <- summary_df$std[i]/sqrt(all_experiments_df[all_experiments_df$Sample == summary_df$Sample[i], 10][1])
# }
# 
# summary_df %>% ggplot(aes(x = Timepoint, y = avg, color = Sample)) + geom_line() + geom_point()
summary_df %>% ggplot(aes(x = Timepoint, y = avg, color = Sample)) + geom_line() + geom_point() + geom_errorbar(aes(x=Timepoint, ymin=avg-std, ymax=avg+std), width = .25) +
  ggtitle("ORC ATPase Timecourse ") +
  xlab("Time (mins)") + ylab("Percent") +
  theme_classic()
# 
# summary_df
# 
# summary_df %>% filter(Sample != "No_ORC") %>% ggplot(aes(x = Timepoint, y = avg, color = Sample)) + geom_line() + geom_point() + geom_errorbar(aes(x=Timepoint, ymin=avg-SEM, ymax=avg+SEM), width = .25) +
#   ggtitle("ORC ATPase Timecourse (n=3-5) ") + theme(plot.title = element_text(hjust = 0.5)) +
#   xlab("Time (mins)") + ylab("(pmol ADP)/(pmol ORC)") 
# 
# 
# #Doing after initial lines since I had an error and wanted to make sure this and preivous code work by themselves 
# all_experiments_df <- all_experiments_df %>% mutate(ADP_per_ORC_adj = ADP_per_ORC/25)
# summary_adj_df <- all_experiments_df %>% group_by(Timepoint, Sample) %>% summarise(avg = mean(ADP_per_ORC_adj), std = sd(ADP_per_ORC_adj))
# summary_adj_df['SEM_adj'] <-array(data = NA, length(summary_adj_df$avg))
# for (i in 1:length(summary_adj_df$SEM_adj)){
#   summary_adj_df$SEM_adj[i] <- summary_adj_df$std[i]/sqrt(all_experiments_df[all_experiments_df$Sample == summary_adj_df$Sample[i], 10][1])
# }
# 
# summary_adj_df %>% filter(Sample != "No_ORC") %>% ggplot(aes(x = Timepoint, y = avg, color = Sample)) + geom_line() + geom_point() + geom_errorbar(aes(x=Timepoint, ymin=avg-SEM_adj, ymax=avg+SEM_adj), width = .25) +
#   ggtitle("ORC ATPase Timecourse (n=3-5) ") + theme(plot.title = element_text(hjust = 0.5)) +
#   xlab("Time (mins)") + ylab("(pmol ADP)/(pmol ORC)") 
# 
# 
# 

####Code for Troubleshooting####
# colnames(atpase_data) <- c("Idx", "ADP")
# 
# atpase_data <- mutate(atpase_data,  Even_Odds = Idx %% 2)
# atpase_data['ATP'] <- array(data = NA, length(atpase_data$Idx))
# 
# #For all rows, if row is odd, add ADP number of next row to the ATP row.
# for (i in 1:length(atpase_data$Idx)){
#   if (!atpase_data[,3][i]){
#     atpase_data$ATP[i] <- atpase_data$ADP[i-1]
#   }
# }
# 
# atpase_data <- na.omit(atpase_data)
# atpase_data <- atpase_data[,grep("ADP|ATP", colnames(atpase_data))]
# 
# 
# atpase_data['Timepoint'] <- array(data = NA, length(atpase_data$Idx))
# atpase_data[1, 'Timepoint'] <- 90
# for (i in 2:length(atpase_data$ADP)){
#   if(i %% 4 == 2){
#     atpase_data$Timepoint[i] <- Time_points[1]
#   }
#   else if(i %% 4 == 3){
#     atpase_data$Timepoint[i] <- Time_points[2]
#   }
#   else if(i %% 4 == 0){
#     atpase_data$Timepoint[i] <- Time_points[3]
#   }
#   else if(i %% 4 == 1){
#     atpase_data$Timepoint[i] <- Time_points[4]
#   }
# }
# 
# 
# 
# atpase_data['Sample'] <- array(data = NA, length(atpase_data$Idx))
# atpase_data[1, 'Sample'] <- "No_ORC"
# indices <- c(2,6,10,14)
# tracker = 1
# for (i in indices){
#   if (indices[tracker] < indices[length(indices)]){
#     #print(tracker)
#     while(i < indices[tracker+1]){
#       atpase_data$Sample[i] <- Sample_Names[tracker]
#       i <- i + 1
#     }
#     tracker <- tracker + 1
#   }
# }
# 
# atpase_data$Sample <- as.factor(atpase_data$Sample)
# atpase_data['Experiment'] <- as.factor('Exp 1')
# 
# atpase_data <- mutate(atpase_data, Percent = ADP/(ADP+ATP)) #Calculate ADP out of total intensity in sample
# atpase_data <- mutate(atpase_data, Percent_corr = Percent-Percent[1]) #Background subtraction
# # atpase_data <- mutate(atpase_data, ADP_pmols = Percent_corr*pmols_rad_ATP/active_frac_ATP)
# all_experiments_df <- data.frame()
# all_experiments_df <- rbind(all_experiments_df, atpase_data)
# all_experiments_df <- atpase_data
# atpase_data %>% ggplot(aes(x = Timepoint, y = Percent, color = Sample)) + geom_line() + geom_point()
# uCi_used = 5
# days_past_peak = 2
# Ci_fraction = exp(-log(2)*days_past_peak/14.26)
# Ci_per_mmol = 3000*Ci_fraction
# pmols_rad_ATP = (1/Ci_per_mmol)*(uCi_used*(10^-6))*(10^9)
# pmols_buff_ATP = (100*25)
# active_frac_ATP = pmols_rad_ATP/(pmols_buff_ATP+pmols_rad_ATP)
# pmol_ORC = 4.88# atpase_data <- as.data.frame(atpase_data)


