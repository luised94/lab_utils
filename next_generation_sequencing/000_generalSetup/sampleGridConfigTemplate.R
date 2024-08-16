#categories <- list(
#    category_names = c("rescue_allele", "mcm_tag", "auxin_treatment", "cell_cycle", "antibody")
#    strain_source = c("lemr", "oa")
#    rescue_allele = c("none", "wt")
#    mcm_tag = c("none", "2", "4", "7")
#    mcm_tag = c("none", "2", "7")
#    auxin_treatment = c("no", "yes")
#    cell_cycle = c("G1", "M")
#    antibody = c("Input", "ProtG", "ALFA", "HM1108", "74", "CHA", "11HA")
#)
#
#subsetting_rules <- list(
#    list(column = "treatment", condition = "DMSO", time_condition = "!= 0"),
#    list(column = "treatment", condition = "DMSO", time_condition = "!= 0"),
#)
#

category_names <- c("rescue_allele", "mcm_tag", "auxin_treatment", "cell_cycle", "antibody")
strain_source <- c("lemr", "oa")
rescue_allele <- c("none", "wt")
mcm_tag <- c("none", "2", "4", "7")
mcm_tag <- c("none", "2", "7")
auxin_treatment <- c("no", "yes")
cell_cycle <- c("G1", "M")
antibody <- c("Input", "ProtG", "ALFA", "HM1108", "74", "CHA", "11HA")

variables_to_exclude <- "category_names|experiment_dir|args|dropbox_dir|variables_to_exclude"
print(ls()[!grepl(variables_to_exclude, ls())])
variables_for_sample_grid <- ls()[!grepl(variables_to_exclude, ls())]
sample_combinations <- expand.grid(lapply(variables_for_sample_grid, get))
names(sample_combinations) <- variables_for_sample_grid
sample_combinations <- sample_combinations[,c("strain_source", "rescue_allele", "mcm_tag", "auxin_treatment", "cell_cycle", "antibody")]

is_input <- sample_combinations[, "rescue_allele"] == "none" &
            sample_combinations[, "mcm_tag"] == "none" &
            sample_combinations[, "cell_cycle"] == "M" &
            sample_combinations[, "antibody"] == "Input" &
            ((sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "auxin_treatment"] == "no") | (sample_combinations[, "strain_source"] == "lemr"))

print(sum(is_input))
#print(sample_combinations[is_input,])


is_protg <- sample_combinations[, "rescue_allele"] == "wt" &
            sample_combinations[, "mcm_tag"] == "none" &
            sample_combinations[, "cell_cycle"] == "M" &
            sample_combinations[, "antibody"] == "ProtG" &
            sample_combinations[, "strain_source"] == "oa" &
            sample_combinations[, "auxin_treatment"] == "no"

print(sum(is_protg))
#print(sample_combinations[is_input,])
#print(sample_combinations[is_input,])


is_alfa <- sample_combinations[, "rescue_allele"] == "none" &
            sample_combinations[, "mcm_tag"] == "none" &
            sample_combinations[, "cell_cycle"] == "M" &
            sample_combinations[, "antibody"] == "ALFA" &
            ((sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "auxin_treatment"] == "no") | (sample_combinations[, "strain_source"] == "lemr"))

print(sum(is_alfa))
#print(sample_combinations[is_alfa,])


is_1108 <- sample_combinations[, "rescue_allele"] == "none" &
            sample_combinations[, "mcm_tag"] == "none" &
            sample_combinations[, "cell_cycle"] == "M" &
            sample_combinations[, "antibody"] == "HM1108" &
            ((sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "auxin_treatment"] == "no") | (sample_combinations[, "strain_source"] == "lemr"))


print(sum(is_1108))
#print(sample_combinations[is_1108,])


is_174 <- sample_combinations[, "antibody"] == "74" &
          sample_combinations[, "auxin_treatment"] == "no" &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "rescue_allele"] == "none") &
          !(sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "rescue_allele"] == "wt") &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "mcm_tag"] == "7") &
          !(sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "mcm_tag"] == "2")
print(sum(is_174))
#print(sample_combinations[is_174,])

is_cha <- sample_combinations[, "antibody"] == "CHA" &
          sample_combinations[, "auxin_treatment"] == "no" &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "rescue_allele"] == "none") &
          !(sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "rescue_allele"] == "wt") &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "mcm_tag"] == "7") &
          !(sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "mcm_tag"] == "2")
print(sum(is_cha))
#print(sample_combinations[is_cha,])


is_11HA <- sample_combinations[, "antibody"] == "11HA" &
          sample_combinations[, "auxin_treatment"] == "no" &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "rescue_allele"] == "none") &
          !(sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "rescue_allele"] == "wt") &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "mcm_tag"] == "7") &
          !(sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "mcm_tag"] == "2") &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "mcm_tag"] == "none" & sample_combinations[, "rescue_allele"] == "wt" & sample_combinations[, "cell_cycle"] == "M")
print(sum(is_11HA))
#print(sample_combinations[is_11HA,])

