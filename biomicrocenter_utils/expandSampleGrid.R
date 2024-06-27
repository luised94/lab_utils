rescue_allele <- c("none","WT","4R")
suppressor_allele <- c("none","1EK","4PS","T2")
cell_cycle <- c("G1","G2")
treatment <- c("No","Yes")
antibody <- c("ALFA","V5", "174", "185", "ProtG")

none_samples <- expand.grid(rescue_allele[1],suppressor_allele[1],cell_cycle,treatment,antibody)
 wt_samples <- expand.grid(rescue_allele[2],suppressor_allele[1],cell_cycle,treatment[2],antibody)
 fourr_samples <- expand.grid(rescue_allele[3],suppressor_allele,cell_cycle,treatment[2],antibody)
 one185wt_samples <- expand.grid(rescue_allele[2],suppressor_allele[1],cell_cycle,treatment[2],"185")
 one185none_samples <- expand.grid(rescue_allele[1],suppressor_allele[1],cell_cycle, treatment,"185")
one185fourr_samples <- expand.grid(rescue_allele[3],suppressor_allele[c(1,3)],cell_cycle,treatment[2],"185")

all_samples <- rbind(none_samples,wt_samples,fourr_samples, one185wt_samples,one185fourr_samples, one185none_samples)
number_of_antibody_samples <- nrow(all_samples)
number_of_input <- nrow(unique(all_samples[,1:2]))
number_of_protg <- 1
#print(unique(all_samples[,1:2]))
#print(number_of_input)
#print(sum(number_of_antibody_samples,number_of_input,number_of_protg))

input_samples <- expand.grid( rescue_allele, suppressor_allele, cell_cycle[1], treatment[2], "INP")
none_samples <- input_samples[input_samples[ , 1 ] == "none" & input_samples[ , 2 ] == "none", ]
wt_samples <- input_samples[input_samples[ , 1 ] == "WT" & input_samples[ , 2] == "none", ]
fourr_samples <- input_samples[input_samples[ , 1] == "4R" , ]
#print(input_samples)
#print(input_samples[, 1 ] == "none") 
#print(none_samples)
#print(wt_samples)
#print(fourr_samples)

input_samples <- rbind(none_samples, wt_samples, fourr_samples)
print(input_samples)

all_samples <- rbind(all_samples, input_samples)
print(nrow(all_samples))

protg_sample <- c("none", "none", "G1", "No", "ProtG") 
levels(all_samples[ , 5]) <- c(levels(all_samples[ , 5]), "ProtG")
all_samples <- rbind(all_samples, protg_sample)
print(nrow(all_samples))

all_samples[ , 6] <- apply(all_samples, 1, function(row) paste(row, collapse = "_"))

#write.table(all_samples, "samples.tsv", sep="\t", row.names=FALSE)
