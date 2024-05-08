

df_sample_info <- read.delim("./documentation/sample_info_table.csv", header = TRUE, sep = ",")
sample_ids <- df_sample_info$sample_ID
base_dir <- "./qualityControl"

qualityControl_dir <- list.dirs("./qualityControl", recursive = FALSE)
contains_sample_id <- grepl(sample_ids[1], list.dirs("./qualityControl", recursive = FALSE))
subset_qualityControl_dir <- qualityControl_dir[contains_sample_id][1]
fastqc_file_path <- list.files(subset_qualityControl_dir, pattern = "fastqc_data", full.names = TRUE)	
summary_file_path <- list.files(subset_qualityControl_dir, pattern = "summary", full.names = TRUE)	

data_blocks <- list()
con <- file(fastqc_file_path, open = "r")
counter <- 0
while (TRUE) {
  message <- sprintf("Loop counter: %s", counter)
  print(message)
  line <- readLines(con, n = 1, warn = FALSE)
#  line_content <- sprintf("Content of line is: %s", line)
#  print(line_content)

  if (length(line) == 0) { 
	print("Line is length 0, Exiting")
	break
  }
  counter <- counter + 1
  
  if (grepl(">>", line)) {
    block_name <- sub(">> ", "", line)
    data_chunk <- ""
    message <- sprintf("Block name is %s", block_name)
    print(message)
    while (TRUE) {
      line <- readLines(con, n = 1)
      if (grepl(">>END_MODULE", line)) { 
	message <- "Reached end of module"
        print(message)
	break 
	}
      data_chunk <- paste(data_chunk, line, sep = "\n")
      print(data_chunk)
    }
    data_blocks[[block_name]] <- data_chunk
  }
}

close(con)
head(data_blocks)
#  message <- sprintf("Loop counter: %s", counter)
#  print(message)
#  line <- readLines(con, n = 1, warn = FALSE)
#  line_content <- sprintf("Content of line is: %s", line)
#  print(line_content)
#  if (length(line) == 0) { 
#	print("Line is length 0, Exiting")
#	break
#  }
#  if (grepl(">>", line)) {
#    block_name <- sub(">> ", "", line)
#    data_chunk <- ""
#    message <- sprintf("Block name is %s", block_name)
#    print(message)
#    while (TRUE) {
#      line <- readLines(con, n = 1)
#      if (grepl(">>END_MODULE", line)) { 
#	message <- "Reached end of module"
#        print(message)
#	break 
#	}
#      data_chunk <- paste(data_chunk, line, sep = "\n")
#      print(data_chunk)
#    }
#    data_blocks[[block_name]] <- data_chunk
#  }
#  counter <- counter + 1
parse_summary <- function(file_path) {
#TODO after working out logic of reading summary file, convert into function, add more descriptive name and Roxygen2 style documentation
}
