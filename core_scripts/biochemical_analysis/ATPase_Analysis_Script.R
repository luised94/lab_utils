library(xlsx)
library(tidyverse)
library(ggplot2)
#make sure you install 64-bit java
#rm(list = ls())
#Copy and pasted code to analyze all experiments. 
#### Paths to files ####
#Use variables with Luis when in lab and liusm when at home. 
path_1 = "C:/Users/Luis/Desktop/2020_09_03 ATPase Analysis of 4R supps/2020_07_13 Orc3 P481L 4R/2020_07_13 Orc3 P481L 4R ATPase.xlsx"
path_2 = "C:/Users/Luis/Desktop/2020_09_03 ATPase Analysis of 4R supps/2020_07_27 Orc3 P481L 4R #2/2020_07_27 Orc3 P481L 4R ATPase #2.xlsx"
path_3 = "C:/Users/Luis/Desktop/2020_09_03 ATPase Analysis of 4R supps/2020_08_13 Orc3 PL 4R #3 and Orc4 PS 4R/2020_08_13 Orc3 PL 4R #3 and Orc4 PS 4R #1.xlsx"
path_4 = "C:/Users/Luis/Desktop/2020_09_03 ATPase Analysis of 4R supps/2020_08_24 Orc4 PS 4R #2 and #3/2020_08_24 Orc4 PS 4R #2 and #3.xlsx"
path_5 = "C:/Users/Luis/Desktop/2020_09_03 ATPase Analysis of 4R supps/2020_09_21 Orc1 EK 4R/2020_09_21 Orc1 EK 4R.xlsx"

# path_1 = "C:/Users/liusm/Desktop/2020_09_03 ATPase Analysis of 4R supps/2020_07_13 Orc3 P481L 4R/2020_07_13 Orc3 P481L 4R ATPase.xlsx"
# path_2 = "C:/Users/liusm/Desktop/2020_09_03 ATPase Analysis of 4R supps/2020_07_27 Orc3 P481L 4R #2/2020_07_27 Orc3 P481L 4R ATPase #2.xlsx"
# path_3 = "C:/Users/liusm/Desktop/2020_09_03 ATPase Analysis of 4R supps/2020_08_13 Orc3 PL 4R #3 and Orc4 PS 4R/2020_08_13 Orc3 PL 4R #3 and Orc4 PS 4R #1.xlsx"
# path_4 = "C:/Users/liusm/Desktop/2020_09_03 ATPase Analysis of 4R supps/2020_08_24 Orc4 PS 4R #2 and #3/2020_08_24 Orc4 PS 4R #2 and #3.xlsx"
# path_5 = "C:/Users/liusm/Desktop/2020_09_03 ATPase Analysis of 4R supps/2020_09_21 Orc1 EK 4R/2020_09_21 Orc1 EK 4R.xlsx"

###### First Experiment ######
#Calculations for radioactivity 
uCi_used = 5
days_past_peak = 2
Ci_fraction = exp(-log(2)*days_past_peak/14.26)
Ci_per_mmol = 3000*Ci_fraction
pmols_rad_ATP = (1/Ci_per_mmol)*(uCi_used*(10^-6))*(10^9)
pmols_buff_ATP = (100*25)
active_frac_ATP = pmols_rad_ATP/(pmols_buff_ATP+pmols_rad_ATP)
pmol_ORC = 4.88

atpase_data = read.xlsx(path_1, header = FALSE, sheetIndex = 4) 
atpase_data <- as.data.frame(atpase_data)
colnames(atpase_data) <- c("Idx", "ADP")

atpase_data <- mutate(atpase_data,  Even_Odds = Idx %% 2)
atpase_data['ATP'] <-array(data = NA, length(atpase_data$Idx)) 

#For all rows, if row is odd, add ADP number of next row to the ATP row. 
for (i in 1:length(atpase_data$Idx)){
  if (atpase_data[,3][i]){
    atpase_data$ATP[i] <- atpase_data$ADP[i+1]
  }
}

atpase_data <- na.omit(atpase_data)
atpase_data <- atpase_data[,c(2,4)]

Time_points = c(0, 15, 45, 90)
atpase_data['Timepoint'] <-array(data = NA, length(atpase_data$Idx)) 
for (i in 1:length(atpase_data$ADP)){
  if(i %% 4 == 1){
    atpase_data$Timepoint[i] <- Time_points[1]
  }
  else if(i %% 4 == 2){
    atpase_data$Timepoint[i] <- Time_points[2]
  }
  else if(i %% 4 == 3){
    atpase_data$Timepoint[i] <- Time_points[3]
  }
  else if(i %% 4 == 0){
    atpase_data$Timepoint[i] <- Time_points[4]
  }
}

Sample_Names <- c("No_ORC", "WT", "WT/DNA","Orc3_P481L_4R")
atpase_data['Sample'] <- array(data = NA, length(atpase_data$Idx)) 
indices <- c(1,5,9,13,17)
tracker = 1
for (i in indices){
  if (indices[tracker] < 17){
    #print(tracker)
    while(i < indices[tracker+1]){
      atpase_data$Sample[i] <- Sample_Names[tracker]
      i <- i + 1
    }
    tracker <- tracker + 1 
  }
}
atpase_data$Sample <- as.factor(atpase_data$Sample)
atpase_data['Experiment'] <- as.factor('Exp 1')

atpase_data <- mutate(atpase_data, Percent = ADP/(ADP+ATP)) #Calculate ADP out of total intensity in sample
atpase_data <- mutate(atpase_data, Percent_corr = Percent-Percent[4]) #Background subtraction 
atpase_data <- mutate(atpase_data, ADP_pmols = Percent_corr*pmols_rad_ATP/active_frac_ATP)

all_experiments_df <- atpase_data
atpase_data %>% ggplot(aes(x = Timepoint, y = Percent, color = Sample)) + geom_line() + geom_point()
atpase_data %>% ggplot(aes(x = Timepoint, y = ADP_pmols, color = Sample)) + geom_line() + geom_point()
####### Second Experiment 2020_07_27 ######
uCi_used = 5; days_past_peak = 10;Ci_fraction = exp(-log(2)*days_past_peak/14.26)
Ci_per_mmol = 3000*Ci_fraction; pmols_rad_ATP = (1/Ci_per_mmol)*(uCi_used*(10^-6))*(10^9)
pmols_buff_ATP = (100*25); active_frac_ATP = pmols_rad_ATP/(pmols_buff_ATP+pmols_rad_ATP)
pmol_ORC = 4.88

atpase_data = read.xlsx(path_2, header = FALSE, sheetIndex = 4)
atpase_data <- as.data.frame(atpase_data)
colnames(atpase_data) <- c("Idx", "ADP")

atpase_data <- mutate(atpase_data,  Even_Odds = Idx %% 2)
atpase_data['ATP'] <-array(data = NA, length(atpase_data$Idx)) 

#For all rows, if row is odd, add ADP number of next row to the ATP row. 
for (i in 1:length(atpase_data$Idx)){
  if (atpase_data[,3][i]){
    atpase_data$ATP[i] <- atpase_data$ADP[i+1]
  }
}

atpase_data <- na.omit(atpase_data)
atpase_data <- atpase_data[,c(2,4)]

Time_points = c(0, 15, 45, 90)
atpase_data['Timepoint'] <-array(data = NA, length(atpase_data$Idx)) 
for (i in 1:length(atpase_data$ADP)){
  if(i %% 4 == 1){
    atpase_data$Timepoint[i] <- Time_points[1]
  }
  else if(i %% 4 == 2){
    atpase_data$Timepoint[i] <- Time_points[2]
  }
  else if(i %% 4 == 3){
    atpase_data$Timepoint[i] <- Time_points[3]
  }
  else if(i %% 4 == 0){
    atpase_data$Timepoint[i] <- Time_points[4]
  }
}
atpase_data$Timepoint[17] <- 90
Sample_Names <- c("WT", "WT/DNA", "4R","Orc3_P481L_4R")
atpase_data['Sample'] <- array(data = NA, length(atpase_data$Idx)) 
indices <- c(1,5,9,13,17)
tracker = 1
for (i in indices){
  if (indices[tracker] < 17){
    #print(tracker)
    while(i < indices[tracker+1]){
      atpase_data$Sample[i] <- Sample_Names[tracker]
      i <- i + 1
    }
    tracker <- tracker + 1 
  }
}
atpase_data$Sample[17] <- "No_ORC"
atpase_data$Sample <- as.factor(atpase_data$Sample)
atpase_data['Experiment'] <- as.factor('Exp 2')

atpase_data <- mutate(atpase_data, Percent = ADP/(ADP+ATP)) #Calculate ADP out of total intensity in sample
atpase_data <- mutate(atpase_data, Percent_corr = Percent-Percent[17]) #Background subtraction 
atpase_data <- mutate(atpase_data, ADP_pmols = Percent_corr*pmols_rad_ATP/active_frac_ATP)
 
all_experiments_df <- bind_rows(all_experiments_df, atpase_data)

atpase_data %>% ggplot(aes(x = Timepoint, y = Percent, color = Sample)) + geom_line() + geom_point()
atpase_data %>% ggplot(aes(x = Timepoint, y = ADP_pmols, color = Sample)) + geom_line() + geom_point()

##### Third Experiment 2020_08_13 Orc3 PL 4R #3 #####
uCi_used = 5; days_past_peak = -8;Ci_fraction = exp(-log(2)*days_past_peak/14.26)
Ci_per_mmol = 3000*Ci_fraction; pmols_rad_ATP = (1/Ci_per_mmol)*(uCi_used*(10^-6))*(10^9)
pmols_buff_ATP = (100*25); active_frac_ATP = pmols_rad_ATP/(pmols_buff_ATP+pmols_rad_ATP)
pmol_ORC = 4.88

atpase_data = read.xlsx(path_3, header = FALSE, sheetIndex = 4)
atpase_data <- as.data.frame(atpase_data)
colnames(atpase_data) <- c("Idx", "ADP")

atpase_data <- mutate(atpase_data,  Even_Odds = Idx %% 2)
atpase_data['ATP'] <-array(data = NA, length(atpase_data$Idx)) 

#For all rows, if row is odd, add ADP number of next row to the ATP row. 
for (i in 1:length(atpase_data$Idx)){
  if (atpase_data[,3][i]){
    atpase_data$ATP[i] <- atpase_data$ADP[i+1]
  }
}

atpase_data <- na.omit(atpase_data)
atpase_data <- atpase_data[,c(2,4)]


Time_points = c(0, 15, 45, 90)
atpase_data['Timepoint'] <-array(data = NA, length(atpase_data$Idx)) 
for (i in 1:length(atpase_data$ADP)){
  if(i %% 4 == 1){
    atpase_data$Timepoint[i] <- Time_points[1]
  }
  else if(i %% 4 == 2){
    atpase_data$Timepoint[i] <- Time_points[2]
  }
  else if(i %% 4 == 3){
    atpase_data$Timepoint[i] <- Time_points[3]
  }
  else if(i %% 4 == 0){
    atpase_data$Timepoint[i] <- Time_points[4]
  }
}
atpase_data$Timepoint[17] <- 90
Sample_Names <- c("WT", "WT/DNA", "4R","Orc3_P481L_4R")
atpase_data['Sample'] <- array(data = NA, length(atpase_data$Idx)) 
indices <- c(1,5,9,13,17)
tracker = 1
for (i in indices){
  if (indices[tracker] < 17){
    #print(tracker)
    while(i < indices[tracker+1]){
      atpase_data$Sample[i] <- Sample_Names[tracker]
      i <- i + 1
    }
    tracker <- tracker + 1 
  }
}
atpase_data$Sample[17] <- "No_ORC"
atpase_data$Sample <- as.factor(atpase_data$Sample)
atpase_data['Experiment'] <- as.factor('Exp 3')

atpase_data <- mutate(atpase_data, Percent = ADP/(ADP+ATP)) #Calculate ADP out of total intensity in sample
atpase_data <- mutate(atpase_data, Percent_corr = Percent-Percent[17]) #Background subtraction 
atpase_data <- mutate(atpase_data, ADP_pmols = Percent_corr*pmols_rad_ATP/active_frac_ATP)

all_experiments_df <- bind_rows(all_experiments_df, atpase_data)
atpase_data %>% ggplot(aes(x = Timepoint, y = Percent, color = Sample)) + geom_line() + geom_point()
atpase_data %>% ggplot(aes(x = Timepoint, y = ADP_pmols, color = Sample)) + geom_line() + geom_point()
s
##### Third Experiment 2020_08_13 Orc4 PS 4R #1 #####
uCi_used = 5; days_past_peak = -8;Ci_fraction = exp(-log(2)*days_past_peak/14.26)
Ci_per_mmol = 3000*Ci_fraction; pmols_rad_ATP = (1/Ci_per_mmol)*(uCi_used*(10^-6))*(10^9)
pmols_buff_ATP = (100*25); active_frac_ATP = pmols_rad_ATP/(pmols_buff_ATP+pmols_rad_ATP)
pmol_ORC = 4.88

atpase_data = read.xlsx(path_3, header = FALSE, sheetIndex = 5)
atpase_data <- as.data.frame(atpase_data)
colnames(atpase_data) <- c("Idx", "ADP")

atpase_data <- mutate(atpase_data,  Even_Odds = Idx %% 2)
atpase_data['ATP'] <-array(data = NA, length(atpase_data$Idx)) 

#For all rows, if row is odd, add ADP number of next row to the ATP row. 
for (i in 1:length(atpase_data$Idx)){
  if (atpase_data[,3][i]){
    atpase_data$ATP[i] <- atpase_data$ADP[i+1]
  }
}

atpase_data <- na.omit(atpase_data)
atpase_data <- atpase_data[,c(2,4)]

Time_points = c(0, 15, 45, 90)
atpase_data['Timepoint'] <-array(data = NA, length(atpase_data$Idx)) 
for (i in 1:length(atpase_data$ADP)){
  if(i %% 4 == 1){
    atpase_data$Timepoint[i] <- Time_points[1]
  }
  else if(i %% 4 == 2){
    atpase_data$Timepoint[i] <- Time_points[2]
  }
  else if(i %% 4 == 3){
    atpase_data$Timepoint[i] <- Time_points[3]
  }
  else if(i %% 4 == 0){
    atpase_data$Timepoint[i] <- Time_points[4]
  }
}
atpase_data$Timepoint[17] <- 90
Sample_Names <- c("WT", "WT/DNA", "4R","Orc4_P225S_4R")
atpase_data['Sample'] <- array(data = NA, length(atpase_data$Idx)) 
indices <- c(1,5,9,13,17)
tracker = 1
for (i in indices){
  if (indices[tracker] < 17){
    #print(tracker)
    while(i < indices[tracker+1]){
      atpase_data$Sample[i] <- Sample_Names[tracker]
      i <- i + 1
    }
    tracker <- tracker + 1 
  }
}
atpase_data$Sample[17] <- "No_ORC"
atpase_data$Sample <- as.factor(atpase_data$Sample)
atpase_data['Experiment'] <- as.factor('Exp 3.1')

atpase_data <- mutate(atpase_data, Percent = ADP/(ADP+ATP)) #Calculate ADP out of total intensity in sample
atpase_data <- mutate(atpase_data, Percent_corr = Percent-Percent[17]) #Background subtraction 
atpase_data <- mutate(atpase_data, ADP_pmols = Percent_corr*pmols_rad_ATP/active_frac_ATP)

all_experiments_df <- bind_rows(all_experiments_df, atpase_data)
atpase_data %>% ggplot(aes(x = Timepoint, y = Percent, color = Sample)) + geom_line() + geom_point()
atpase_data %>% ggplot(aes(x = Timepoint, y = ADP_pmols, color = Sample)) + geom_line() + geom_point()

##### Fourth Experiment 2020_08_24 Orc4 PS 4R #2 and #3 #####
uCi_used = 5; days_past_peak = 2;Ci_fraction = exp(-log(2)*days_past_peak/14.26)
Ci_per_mmol = 3000*Ci_fraction; pmols_rad_ATP = (1/Ci_per_mmol)*(uCi_used*(10^-6))*(10^9)
pmols_buff_ATP = (100*25); active_frac_ATP = pmols_rad_ATP/(pmols_buff_ATP+pmols_rad_ATP)
pmol_ORC = 4.88

atpase_data = read.xlsx(path_4, header = FALSE, sheetIndex = 4)
atpase_data <- as.data.frame(atpase_data)
colnames(atpase_data) <- c("Idx", "ADP")

atpase_data <- mutate(atpase_data,  Even_Odds = Idx %% 2)
atpase_data['ATP'] <-array(data = NA, length(atpase_data$Idx)) 

#For all rows, if row is odd, add ADP number of next row to the ATP row. 
for (i in 1:length(atpase_data$Idx)){
  if (atpase_data[,3][i]){
    atpase_data$ATP[i] <- atpase_data$ADP[i+1]
  }
}

atpase_data <- na.omit(atpase_data)
atpase_data <- atpase_data[,c(2,4)]

Time_points = c(0, 15, 45, 90)
atpase_data['Timepoint'] <-array(data = NA, length(atpase_data$Idx)) 
for (i in 1:length(atpase_data$ADP)){
  if(i %% 4 == 1){
    atpase_data$Timepoint[i] <- Time_points[1]
  }
  else if(i %% 4 == 2){
    atpase_data$Timepoint[i] <- Time_points[2]
  }
  else if(i %% 4 == 3){
    atpase_data$Timepoint[i] <- Time_points[3]
  }
  else if(i %% 4 == 0){
    atpase_data$Timepoint[i] <- Time_points[4]
  }
}
atpase_data$Timepoint[17] <- 90
Sample_Names <- c("WT", "4R", "Orc4_P225S_4R","Orc4_P225S_4R")
atpase_data['Sample'] <- array(data = NA, length(atpase_data$Idx)) 
indices <- c(1,5,9,13,17)
tracker = 1
for (i in indices){
  if (indices[tracker] < 17){
    #print(tracker)
    while(i < indices[tracker+1]){
      atpase_data$Sample[i] <- Sample_Names[tracker]
      i <- i + 1
    }
    tracker <- tracker + 1 
  }
}
atpase_data$Sample[17] <- "No_ORC"
atpase_data$Sample <- as.factor(atpase_data$Sample)
atpase_data['Experiment'] <- as.factor('Exp 4')

atpase_data <- mutate(atpase_data, Percent = ADP/(ADP+ATP)) #Calculate ADP out of total intensity in sample
atpase_data <- mutate(atpase_data, Percent_corr = Percent-Percent[17]) #Background subtraction 
atpase_data <- mutate(atpase_data, ADP_pmols = Percent_corr*pmols_rad_ATP/active_frac_ATP)

all_experiments_df <- bind_rows(all_experiments_df, atpase_data)
atpase_data %>% ggplot(aes(x = Timepoint, y = Percent, color = Sample)) + geom_line() + geom_point()
atpase_data %>% ggplot(aes(x = Timepoint, y = ADP_pmols, color = Sample)) + geom_line() + geom_point()

##### Fifth Experiment 2020_09_21 Orc1 EK 4R #1#####
uCi_used = 5; days_past_peak = 2;Ci_fraction = exp(-log(2)*days_past_peak/14.26)
Ci_per_mmol = 3000*Ci_fraction; pmols_rad_ATP = (1/Ci_per_mmol)*(uCi_used*(10^-6))*(10^9)
pmols_buff_ATP = (100*25); active_frac_ATP = pmols_rad_ATP/(pmols_buff_ATP+pmols_rad_ATP)
pmol_ORC = 4.88

atpase_data = read.xlsx(path_5, header = FALSE, sheetIndex = 4)
atpase_data <- as.data.frame(atpase_data)
colnames(atpase_data) <- c("Idx","ADP")

atpase_data <- mutate(atpase_data,  Even_Odds = Idx %% 2)
atpase_data['ATP'] <-array(data = NA, length(atpase_data$Idx)) 

#For all rows, if row is odd, add ADP number of next row to the ATP row. 
for (i in 1:length(atpase_data$Idx)){
  if (atpase_data[,3][i]){
    atpase_data$ATP[i] <- atpase_data$ADP[i+1]
  }
}

atpase_data <- na.omit(atpase_data)
atpase_data <- atpase_data[,c(2,4)]

Time_points = c(0, 15, 45, 90)
atpase_data['Timepoint'] <-array(data = NA, length(atpase_data$Idx)) 
for (i in 1:length(atpase_data$ADP)){
  if(i %% 4 == 1){
    atpase_data$Timepoint[i] <- Time_points[1]
  }
  else if(i %% 4 == 2){
    atpase_data$Timepoint[i] <- Time_points[2]
  }
  else if(i %% 4 == 3){
    atpase_data$Timepoint[i] <- Time_points[3]
  }
  else if(i %% 4 == 0){
    atpase_data$Timepoint[i] <- Time_points[4]
  }
}
atpase_data$Timepoint[17] <- 90
Sample_Names <- c("WT", "WT/DNA", "4R","Orc1_E495K_4R")
atpase_data['Sample'] <- array(data = NA, length(atpase_data$Idx)) 
indices <- c(1,5,9,13,17)
tracker = 1
for (i in indices){
  if (indices[tracker] < 17){
    #print(tracker)
    while(i < indices[tracker+1]){
      atpase_data$Sample[i] <- Sample_Names[tracker]
      i <- i + 1
    }
    tracker <- tracker + 1 
  }
}
atpase_data$Sample[17] <- "No_ORC"
atpase_data$Sample <- as.factor(atpase_data$Sample)
atpase_data['Experiment'] <- as.factor('Exp 5')

atpase_data <- mutate(atpase_data, Percent = ADP/(ADP+ATP)) #Calculate ADP out of total intensity in sample
atpase_data <- mutate(atpase_data, Percent_corr = Percent-Percent[17]) #Background subtraction 
atpase_data <- mutate(atpase_data, ADP_pmols = Percent_corr*pmols_rad_ATP/active_frac_ATP)

all_experiments_df <- bind_rows(all_experiments_df, atpase_data)
atpase_data %>% ggplot(aes(x = Timepoint, y = Percent, color = Sample)) + geom_line() + geom_point()
atpase_data %>% ggplot(aes(x = Timepoint, y = ADP_pmols, color = Sample)) + geom_line() + geom_point()

##### Fifth Experiment 2020_09_21 Orc1 EK 4R #2#####
uCi_used = 5; days_past_peak = 2;Ci_fraction = exp(-log(2)*days_past_peak/14.26)
Ci_per_mmol = 3000*Ci_fraction; pmols_rad_ATP = (1/Ci_per_mmol)*(uCi_used*(10^-6))*(10^9)
pmols_buff_ATP = (100*25); active_frac_ATP = pmols_rad_ATP/(pmols_buff_ATP+pmols_rad_ATP)
pmol_ORC = 4.88

atpase_data = read.xlsx(path_5, header = FALSE, sheetIndex = 5)
atpase_data <- as.data.frame(atpase_data)
colnames(atpase_data) <- c("Idx","ADP")

atpase_data <- mutate(atpase_data,  Even_Odds = Idx %% 2)
atpase_data['ATP'] <-array(data = NA, length(atpase_data$Idx)) 

#For all rows, if row is odd, add ADP number of next row to the ATP row. 
for (i in 1:length(atpase_data$Idx)){
  if (atpase_data[,3][i]){
    atpase_data$ATP[i] <- atpase_data$ADP[i+1]
  }
}

atpase_data <- na.omit(atpase_data)
atpase_data <- atpase_data[,c(2,4)]

Time_points = c(0, 15, 45, 90)
atpase_data['Timepoint'] <-array(data = NA, length(atpase_data$Idx)) 
for (i in 1:length(atpase_data$ADP)){
  if(i %% 4 == 1){
    atpase_data$Timepoint[i] <- Time_points[1]
  }
  else if(i %% 4 == 2){
    atpase_data$Timepoint[i] <- Time_points[2]
  }
  else if(i %% 4 == 3){
    atpase_data$Timepoint[i] <- Time_points[3]
  }
  else if(i %% 4 == 0){
    atpase_data$Timepoint[i] <- Time_points[4]
  }
}
atpase_data$Timepoint[17] <- 90
Sample_Names <- c("WT", "WT/DNA", "4R","Orc1_E495K_4R")
atpase_data['Sample'] <- array(data = NA, length(atpase_data$Idx)) 
indices <- c(1,5,9,13,17)
tracker = 1
for (i in indices){
  if (indices[tracker] < 17){
    #print(tracker)
    while(i < indices[tracker+1]){
      atpase_data$Sample[i] <- Sample_Names[tracker]
      i <- i + 1
    }
    tracker <- tracker + 1 
  }
}
atpase_data$Sample[17] <- "No_ORC"
atpase_data$Sample <- as.factor(atpase_data$Sample)
atpase_data['Experiment'] <- as.factor('Exp 5.1')

atpase_data <- mutate(atpase_data, Percent = ADP/(ADP+ATP)) #Calculate ADP out of total intensity in sample
atpase_data <- mutate(atpase_data, Percent_corr = Percent-Percent[17]) #Background subtraction 
atpase_data <- mutate(atpase_data, ADP_pmols = Percent_corr*pmols_rad_ATP/active_frac_ATP)

all_experiments_df <- bind_rows(all_experiments_df, atpase_data)
atpase_data %>% ggplot(aes(x = Timepoint, y = Percent, color = Sample)) + geom_line() + geom_point()
atpase_data %>% ggplot(aes(x = Timepoint, y = ADP_pmols, color = Sample)) + geom_line() + geom_point()

##### Analysis #####
pmol_ORC = 5
number_of_samples = count(all_experiments_df, Sample, Timepoint)
all_experiments_df <- all_experiments_df %>% mutate(ADP_per_ORC = ADP_pmols/pmol_ORC)
all_experiments_df['Sample_Size'] <-array(data = NA, length(all_experiments_df$ADP))
for (i in 1:length(all_experiments_df$ADP)){
  all_experiments_df$Sample_Size[i] <- number_of_samples[number_of_samples$Sample == all_experiments_df$Sample[i], 3][1] 
}

summary_df <- all_experiments_df %>% group_by(Timepoint, Sample) %>% summarise(avg = mean(ADP_pmols), std = sd(ADP_pmols))
summary_df <- all_experiments_df %>% group_by(Timepoint, Sample) %>% summarise(avg = mean(ADP_per_ORC), std = sd(ADP_per_ORC))
summary_df['SEM'] <-array(data = NA, length(summary_df$avg))
for (i in 1:length(summary_df$SEM)){
 summary_df$SEM[i] <- summary_df$std[i]/sqrt(all_experiments_df[all_experiments_df$Sample == summary_df$Sample[i], 10][1])
}

summary_df %>% ggplot(aes(x = Timepoint, y = avg, color = Sample)) + geom_line() + geom_point()
summary_df %>% ggplot(aes(x = Timepoint, y = avg, color = Sample)) + geom_line() + geom_point() + geom_errorbar(aes(x=Timepoint, ymin=avg-std, ymax=avg+std), width = .25) +
  ggtitle("ORC ATPase Timecourse ") + 
  xlab("Time (mins)") + ylab("(pmol ADP)/(pmol ORC)") +
  theme_classic()

summary_df %>% filter(Sample != "No_ORC") %>% ggplot(aes(x = Timepoint, y = avg, color = Sample)) + geom_line() + geom_point() + geom_errorbar(aes(x=Timepoint, ymin=avg-SEM, ymax=avg+SEM), width = .25) +
  ggtitle("ORC ATPase Timecourse (n=3-5) ") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time (mins)") + ylab("(pmol ADP)/(pmol ORC)") 
 

#Doing after initial lines since I had an error and wanted to make sure this and preivous code work by themselves 
all_experiments_df <- all_experiments_df %>% mutate(ADP_per_ORC_adj = ADP_per_ORC/25)
summary_adj_df <- all_experiments_df %>% group_by(Timepoint, Sample) %>% summarise(avg = mean(ADP_per_ORC_adj), std = sd(ADP_per_ORC_adj))
summary_adj_df['SEM_adj'] <-array(data = NA, length(summary_adj_df$avg))
for (i in 1:length(summary_adj_df$SEM_adj)){
  summary_adj_df$SEM_adj[i] <- summary_adj_df$std[i]/sqrt(all_experiments_df[all_experiments_df$Sample == summary_adj_df$Sample[i], 10][1])
}

summary_adj_df %>% filter(Sample != "No_ORC") %>% ggplot(aes(x = Timepoint, y = avg, color = Sample)) + geom_line() + geom_point() + geom_errorbar(aes(x=Timepoint, ymin=avg-SEM_adj, ymax=avg+SEM_adj), width = .25) +
  ggtitle("ORC ATPase Timecourse (n=3-5) ") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time (mins)") + ylab("(pmol ADP)/(pmol ORC)") 
 