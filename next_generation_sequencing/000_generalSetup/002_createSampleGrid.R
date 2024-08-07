
rm(list = ls())
category_names <- c("rescue_allele", "mcm_tag", "auxin_treatment", "cell_cycle", "antibody")

strain_source <- c("lemr", "oa")
rescue_allele <- c("none", "wt")
mcm_tag <- c("none", "2", "4", "7")
auxin_treatment <- c("no", "yes")
cell_cycle <- c("G1", "M")
antibody <- c("Input", "ProtG", "ALFA", "HM1108", "74", "CHA", "11HA")

variables_for_sample_grid <- ls()[!grepl("category_names", ls())]
sample_combinations <- expand.grid(lapply(variables_for_sample_grid, get))
names(sample_combinations) <- variables_for_sample_grid
sample_combinations <- sample_combinations[,c("strain_source", "rescue_allele", "mcm_tag", "auxin_treatment", "cell_cycle", "antibody")]

sample_combinations$full_name <- apply(sample_combinations, 1, paste, collapse = "_")
sample_combinations$short_name <- apply(sample_combinations[,!grepl("full_name", colnames(sample_combinations))], 1, function(row) paste0(substr(row, 1, 1), collapse = ""))

is_input <- sample_combinations[, "rescue_allele"] == "none" &
            sample_combinations[, "mcm_tag"] == "none" &
            sample_combinations[, "cell_cycle"] == "M" &
            sample_combinations[, "antibody"] == "Input" &
            ((sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "auxin_treatment"] == "no") | (sample_combinations[, "strain_source"] == "lemr"))

print(sum(is_input))
print(sample_combinations[is_input,])


