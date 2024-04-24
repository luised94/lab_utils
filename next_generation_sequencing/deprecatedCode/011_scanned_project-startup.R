#Initial project setup ----
Sys.time()

#Create directories and subdirectories to store files that will be generated during the analysis.
#This should be run in the directory that will be used for analysis

#Add this to all subdirectories
files_suffix <- "-files"

#Create folders that will have data.
#Can modify this for other types of projects/analysis
folder_data <- "./data"

ngs_data_folders <- c("fastq", "genome", "peak", "alignment","fastqc","processed-fastq", "bw")
data_dirs <- c()

if (!dir.exists(folder_data)){
  dir.create(folder_data)
  
  #Create folders for data that will be generated. 
  lapply(ngs_data_folders, function(x){
    dir.create(paste(paste(folder_data, x, sep = "/"), files_suffix, sep = ""))
  })
  
  #Store folder names for later use.  
  data_dirs <- c(data_dirs, paste(paste(folder_data, ngs_data_folders, sep = "/"), files_suffix, sep = ""))
} else {
  print(paste(folder_data, "folder exists"))
}


#Create folder that holds code used in the project 
folder_code <- "./scripts"

code_folders <- c("python", "R", "command-line")
code_dirs <- c()
#All projects will have code, so no need to check the name
if (!dir.exists(folder_code)){
  dir.create(folder_code)
  
  #Create folders for data that will be generated. 
  lapply(code_folders, function(x){
    dir.create(paste(paste(folder_code, x, sep = "/"), files_suffix, sep = ""))
  })
}

#Create folder to output project reports or figures.
folder_figures <- "./reports"
figure_folders <- c("plots", "reports")
figure_dirs <- c()
#All projects will have code, so no need to check the name
if (!dir.exists(folder_figures)){
  dir.create(folder_figures)
  
  #Create folders for data that will be generated. 
  lapply(figure_folders, function(x){
    dir.create(paste(paste(folder_figures, x, sep = "/"), files_suffix, sep = ""))
  })
}

log_folder <- "./logs"
if (!dir.exists(log_folder)) {
  dir.create(log_folder)
}

print("Folders created")



