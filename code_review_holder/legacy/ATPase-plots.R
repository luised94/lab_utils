scripts_path <- "./script/"
# Function to calculate summarize dataframe subset 
subset_and_summarize_df <- function(data_df, combination, column_string, sample_string) {
  to_subset <- data_df[,column_string] %in% combination & data_df[,'Sample'] == sample_string
  subset_df <- data_df[to_subset, ]
  summary_df <- aggregate(Percent_corr ~ Timepoint + Sample, data = subset_df, 
                          FUN = function(x) c(avg = mean(x), std = sd(x)))
  return(summary_df)
}

# Function to calculate the weighted sum, Negative weight for std, positive weight for avg
calculate_weighted_sum <- function(df, weights = c(1,-1)) {
  #Accessing is determined by output of aggregate function.
  #Extract the avg and std determined by aggregate function. 
  weighted_sum <- sum(as.matrix(df[,3]) %*% weights) #sum(rowSums(as.matrix(df[, 3]) * weights)) 
  return(weighted_sum)
}
#Function to simplify filtering and determining weighted sums of combinations. 
run_experiment_combinations <- function(data_df, experiment_col, experiment_pattern, sample_string, weights = c(1,-1)) {
  experiment_numbers <- unique(data_df[grep(experiment_pattern, data_df$Sample),]$Experiment)
  
  exp_num_comb <- combn(experiment_numbers, 3, simplify = FALSE)
  
  avg_sd_list <- lapply(exp_num_comb, subset_and_summarize_df, data_df = data_df, column_string = experiment_col, sample_string = sample_string)
  
  weighted_sums <- unlist(lapply(avg_sd_list, calculate_weighted_sum, weights = weights))
  
  return(weighted_sums)
}

source(paste(scripts_path, "ATPase-data-process-5,6EK.R", sep = ""))
source(paste(scripts_path, "ATPase-data-process-1,3,4.R", sep = ""))

#Extract columns present in all_experiments_df2
all_experiments_df <- rbind(all_experiments_df1[,names(all_experiments_df2)], all_experiments_df2)
# unique(all_experiments_df$Sample) #Use to make sure all samples were there. 
#Remove other variables cause they arent required. 
rm(list = ls()[!grepl("^all_experiments_df$", ls())])

#Summarize all data and plot to verify. 
summary_df <- all_experiments_df %>% group_by(Timepoint, Sample) %>% summarise(avg = mean(Percent_corr), std = sd(Percent_corr))
summary_df %>% ggplot(aes(x = Timepoint, y = avg, color = Sample)) + geom_line() + geom_point() + geom_errorbar(aes(x=Timepoint, ymin=avg-std, ymax=avg+std), width = .25) +
  ggtitle("ORC ATPase Timecourse ") +
  xlab("Time (mins)") + ylab("Percent") +
  theme_classic()

#Determine all experiments that are 4, 5 and 6.
experiments456 <- lapply(c("Exp 4", "Exp 5", "Exp 6"), function(x){
  grepl(x, all_experiments_df$Experiment)
})
experiments456 <- Reduce(`|`, lapply(experiments456, unlist))
# levels(summary_df$Sample)

#Filter certain experiments since I dont want to plot them.
filtered_df <- all_experiments_df %>% filter(!(Sample == "WT" & experiments456) & !(Sample == "WT/DNA") & !(Sample == "No_ORC") & !(Sample == "Orc5_E104K_4R" & Experiment == "Exp 4"))
summary_df <- filtered_df %>% group_by(Timepoint, Sample) %>% summarise(avg = mean(Percent_corr), std = sd(Percent_corr))
summary_df %>% ggplot(aes(x = Timepoint, y = avg, color = Sample)) + geom_line() + geom_point() + geom_errorbar(aes(x=Timepoint, ymin=avg-std, ymax=avg+std), width = .25) +
    ggtitle("ORC ATPase Timecourse ") +
    xlab("Time (mins)") + ylab("Percent hydrolysis") +
    theme_classic() + scale_colour_brewer(palette="Spectral")

####Determining which combinations for WT to plot####
#Use already filtered to generate combinations to look for. 
experiment_numbers <- unique(filtered_df[grep("^WT$", filtered_df$Sample),]$Experiment)

exp_num_comb <- combn(experiment_numbers, 3, simplify = FALSE)

#Test function and then run for all combinations. 
#subset_and_summarize_df(all_experiments_df, exp_num_comb[[2]], 'Experiment', 'WT')
avg_sd_list <- lapply(exp_num_comb, subset_and_summarize_df, data_df = all_experiments_df, column_string = 'Experiment', sample_string = 'WT')

weighted_sums <- unlist(lapply(avg_sd_list, calculate_weighted_sum))
#run_experiment_combinations(filtered_df, "Experiment", "^WT$", 'WT')

# Find the index of the dataframe with the optimal weighted sum
index_optimal_weighted_sum <- which.max(weighted_sums)
index_worse_weighted_sum <- which.min(weighted_sums)
# Get the corresponding dataframe with the optimal combination for a quick check
# optimal_combination_df <- avg_sd_list[[index_optimal_weighted_sum]]
# worse_combination_df <- avg_sd_list[[index_worse_weighted_sum]]

# For optimal weighted sum
optimal_summary <- all_experiments_df %>%
  filter(Sample == "WT" & Experiment %in% unlist(exp_num_comb[index_optimal_weighted_sum])) %>%
  group_by(Timepoint, Sample) %>%
  summarise(avg = mean(Percent_corr), std = sd(Percent_corr))

# For worse weighted sum
worse_summary <- all_experiments_df %>%
  filter(Sample == "WT" & Experiment %in% unlist(exp_num_comb[index_worse_weighted_sum])) %>%
  group_by(Timepoint, Sample) %>%
  summarise(avg = mean(Percent_corr), std = sd(Percent_corr))

# Combine the dataframes and create a column for the type (optimal or worse)
combined_summary <- rbind(transform(optimal_summary, Type = "Optimal"),
                          transform(worse_summary, Type = "Worse"))

# Plot using ggplot2
ggplot(combined_summary, aes(x = Timepoint, y = avg, group = Type, color = Type)) +
  geom_line() +
  geom_errorbar(aes(ymin = avg - std, ymax = avg + std), width = 0.2) +
  labs(title = "Average and Standard Deviation Comparison",
       x = "Time (s)",
       y = "Percent") +
  theme_classic() 

#Choose experiments for WT
wt_experiment_combination <- unlist(exp_num_comb[index_optimal_weighted_sum])
wt_experiments <- lapply(wt_experiment_combination, function(x){
  grepl(x, all_experiments_df$Experiment)
})
wt_experiments <- Reduce(`|`, lapply(wt_experiments, unlist))
####Repeat analysis for 4R samples####
experiment_numbers <- unique(filtered_df[grep("^4R$", filtered_df$Sample),]$Experiment)
exp_num_comb <- combn(experiment_numbers, 3, simplify = FALSE)
weighted_sums <- run_experiment_combinations(filtered_df, "Experiment", "^4R$", '4R', weights = c(-.2,-1))

# Find the index of the dataframe with the optimal weighted sum
index_optimal_weighted_sum <- which.max(weighted_sums)
index_worse_weighted_sum <- which.min(weighted_sums)
# Get the corresponding dataframe with the optimal combination for a quick check
# optimal_combination_df <- avg_sd_list[[index_optimal_weighted_sum]]
# worse_combination_df <- avg_sd_list[[index_worse_weighted_sum]]

# For optimal weighted sum
optimal_summary <- all_experiments_df %>%
  filter(Sample == "4R" & Experiment %in% unlist(exp_num_comb[index_optimal_weighted_sum])) %>%
  group_by(Timepoint, Sample) %>%
  summarise(avg = mean(Percent_corr), std = sd(Percent_corr))

# For worse weighted sum
worse_summary <- all_experiments_df %>%
  filter(Sample == "4R" & Experiment %in% unlist(exp_num_comb[index_worse_weighted_sum])) %>%
  group_by(Timepoint, Sample) %>%
  summarise(avg = mean(Percent_corr), std = sd(Percent_corr))

# Combine the dataframes and create a column for the type (optimal or worse)
combined_summary <- rbind(transform(optimal_summary, Type = "Optimal"),
                          transform(worse_summary, Type = "Worse"))

# Plot using ggplot2
ggplot(combined_summary, aes(x = Timepoint, y = avg, group = Type, color = Type)) +
  geom_line() +
  geom_errorbar(aes(ymin = avg - std, ymax = avg + std), width = 0.2) +
  labs(title = "Average and Standard Deviation Comparison",
       x = "Time (s)",
       y = "Percent") +
  theme_classic()

#Choose experiments for 4R
fourr_experiment_combination <- unlist(exp_num_comb[index_optimal_weighted_sum])
fourr_experiments <- lapply(fourr_experiment_combination, function(x){
  grepl(x, all_experiments_df$Experiment)
})
fourr_experiments <- Reduce(`|`, lapply(fourr_experiments, unlist))


####Generate plots for thesis####
plotting_df <- all_experiments_df %>% filter(
  (Sample == "WT" & wt_experiments) | (Sample == "4R" & fourr_experiments)) 

plotting_df <- all_experiments_df %>% filter(!(Sample == "WT/DNA") &
                                !(Sample == "No_ORC") &
                                !(Sample == "Orc5_E104K_4R" & Experiment == "Exp 4") &
                                  !(Sample == "WT") & !(Sample == "4R"))
    
plotting_df <- all_experiments_df %>% filter(
  (Sample == "WT" & wt_experiments) |
    (Sample == "4R" & fourr_experiments) | 
    (!(Sample == "WT/DNA") &
       !(Sample == "No_ORC") &
       !(Sample == "Orc5_E104K_4R" & Experiment == "Exp 4") &
       !(Sample == "WT") & !(Sample == "4R")))

summary_plot_df <- plotting_df %>% group_by(Timepoint, Sample) %>%
  summarise(avg = mean(Percent_corr), std = sd(Percent_corr))

# Plot using ggplot2
summary_plot_df %>% ggplot(aes(x = Timepoint, y = avg, color = Sample)) + geom_line() + geom_point() + geom_errorbar(aes(x=Timepoint, ymin=avg-std, ymax=avg+std), width = .25) +
  ggtitle("ORC ATPase Timecourse ") +
  xlab("Time (mins)") + ylab("Percent hydrolysis") +
  theme_classic() + scale_colour_brewer(palette="Dark2")

ggsave("./output/ATPase_plot_classic_dark2.svg")
# summary_df <- all_experiments_df %>% group_by(Timepoint, Sample) %>% summarise(avg = mean(ADP_per_ORC), std = sd(ADP_per_ORC))
# summary_df['SEM'] <-array(data = NA, length(summary_df$avg))
# for (i in 1:length(summary_df$SEM)){
# summary_df$SEM[i] <- summary_df$std[i]/sqrt(all_experiments_df[all_experiments_df$Sample == summary_df$Sample[i], 10][1])
# }
# 
# summary_df %>% ggplot(aes(x = Timepoint, y = avg, color = Sample)) + geom_line() + geom_point()
# weighted_sum <- 0
# for (i in 1:nrow(df)){
#   print(weights[1] * df[i,3][1] + weights[2] * df[i,3][2])
#   weighted_sum <- weighted_sum + weights[1] * df[i,3][1] + weights[2] * df[i,3][2]
# }
# avg_sd_df[[1]][1,3][1]
# df <- avg_sd_df[[1]]
# as.matrix(df[,3]) %*% weights
# weights <- c(1, -1)  # Negative weight for std, positive weight for avg
# Assuming df[, 3] is a 2 by 4 matrix, and weights is a vector with two elements
# weighted_sums <- sum(rowSums(as.matrix(df[, 3]) * weights))
# print(weighted_sums)
# to_subset <- all_experiments_df[,'Experiment']  %in% exp_num_comb[[2]] & all_experiments_df[,'Sample'] =='WT'
# subset_df <- all_experiments_df[to_subset, ]
# summary_df <- aggregate(Percent_corr ~ Timepoint + Sample, data = subset_df, 
#                     FUN = function(x) c(avg = mean(x), std = sd(x)))
# 

