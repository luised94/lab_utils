library(xlsx)
library(tidyverse)

FILENAME <- "Analysis.xlsx"
OUTPUT_DIRECTORY <- "~/data/loading_analysis"
EXPERIMENT_DIRECTORY <- "Lab/Experiments/Loading/2022_12_18 Loading Assays Repeats for publication"
# required columns: Intensity	Lane	ORC 	Suppressor	kGlut	Input
# @NOTE: An additional space in the excel sheet cell with cause column name to have dot at the end.
#   > "ORC " would be read in as "ORC."
REQUIRED_COLUMNS <- c(
  "NA.", "Intensity", "Lane",
  "ORC", "Suppressor",
  "kGlut", "Input"
)

DROPBOX_PATH <- Sys.getenv("DROPBOX_PATH")
FILE_PATH <- file.path(DROPBOX_PATH, EXPERIMENT_DIRECTORY, FILENAME)
if (nchar(DROPBOX_PATH) == 0) {
  stop("DROPBOX_PATH not defined: ")
}


if (!dir.exists(OUTPUT_DIRECTORY)) {
  dir.create(OUTPUT_DIRECTORY, showWarnings = FALSE, recursive = TRUE)
}

if (!file.exists(FILE_PATH)) {
  stop("FILE_PATH does not exist: ", FILE_PATH)
}

message("Paths for input and output set...")
sheet_indices <- c(2, 3, 4)
df_lst <- vector(mode = "list", length = length(sheet_indices))
loading_df <- data.frame(

)
df_count <- 0

message("Reading in loading assay sheets...")
for (sheet_idx in sheet_indices){
  df_count <- df_count + 1
  df_lst[[df_count]] <- read.xlsx(FILE_PATH, header = TRUE, sheetIndex = sheet_idx)
  colnames(df_lst[[df_count]])
  #all(unlist(lapply(REQUIRED_COLUMNS, function(x) {x %in% colnames(df_lst[[df_count]])})))

}

stop("Breakpoint...")

df1 <- read.xlsx(FILE_PATH, header = TRUE, sheetIndex = 2)
df2 <- read.xlsx(FILE_PATH, header = TRUE, sheetIndex =  3) 
df3 <- read.xlsx(FILE_PATH, header = TRUE, sheetIndex =  4)

df1 <- df1 %>%
  mutate(Intensity = Intensity - df1$Intensity[2],
         Experiment = "Exp_1") %>% 
  slice(-2) %>% 
  mutate(Intensity = (Intensity/Intensity[Input == "yes"])*.5) %>%
  filter(Input == "no")

df2 <- df2 %>%
  mutate(Intensity = Intensity - df2$Intensity[2],  
         Experiment = "Exp_2") %>%
  slice(-2) %>%
  mutate(Intensity = (Intensity/Intensity[Input == "yes"])*.5) %>%
  filter(Input == "no")

df3 <- df3 %>%
  mutate(Intensity = Intensity - df3$Intensity[2],
         Experiment = "Exp_3") %>%
  slice(-2) %>% 
  mutate(Intensity = (Intensity/Intensity[Input == "yes"])*.5) %>%
  filter(Input == "no")  

df <- bind_rows(df1, df2, df3)

df <- df %>%
  mutate(Label = case_when(
    ORC. == "WT" & Suppressor == "None" ~ "WT",
    ORC. == "RA" & Suppressor == "None" ~ "ORC4R", 
    ORC. == "RA" & Suppressor == "4PS" ~ "+4sofr",
    ORC. == "RA" & Suppressor == "1EK" ~ "+1sofr",
    ORC. == "RA" & Suppressor == "3PL" ~ "+3sofr"
  ))

# Calculate mean and sd for each kGlut 
df_summary <- df %>% 
  filter(!(Label %in% c("+1sofr","+3sofr"))) %>%
  group_by(kGlut, Label) %>%
  summarise(pmol_MCM = mean(Intensity), 
            sd = sd(Intensity))

# Plot
ggplot(df_summary, aes(x = kGlut, y = pmol_MCM, color = Label)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = pmol_MCM-sd, ymax = pmol_MCM+sd), width = .1) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "kGlut", 
       y = "MCM (pmol)",
       title = "Intensity vs kGlut",
       color = "Protein") +
  theme_classic()
