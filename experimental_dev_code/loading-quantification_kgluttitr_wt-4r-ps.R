library(xlsx)
library(tidyverse)

FILENAME <- "Analysis.xlsx"
OUTPUT_DIRECTORY <- "~/data/loading_analysis"
EXPERIMENT_DIRECTORY <- "Lab/Experiments/Loading/2022_12_18 Loading Assays Repeats for publication"
sheet_indices <- c(2, 3, 4)
EXPECTED_NUMBER_OF_ROWS <- 14
EXPECTED_NUMBER_OF_COLS <- 7
COLUMN_TO_REMOVE <- "NA." # Column has row numbers from imagej export/copy-paste
# Define custom orderings
factor_order <- list(
  "Suppressor" = c("None", "1EK", "3PL", "4PS"),
  "kGlut" = c("250", "300", "350"), # In inverse order (as factor, alphabetical would be "250", ...)
  "Label" = c("WT", "ORC4R",  "+1sofr",  "+3sofr", "+4sofr") # Adjust this as needed
)


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
df_lst <- vector(mode = "list", length = length(sheet_indices))
temp_df <- data.frame(matrix(
  NA,
  nrow = EXPECTED_NUMBER_OF_ROWS,
  ncol = EXPECTED_NUMBER_OF_COLS
))
loading_df <- data.frame(matrix(
  NA,
  nrow = EXPECTED_NUMBER_OF_ROWS * length(sheet_indices),
  ncol = EXPECTED_NUMBER_OF_COLS
))

message("Reading in loading assay sheets...")

df_count <- 0
for (sheet_idx in sheet_indices){
  df_count <- df_count + 1
  temp_df <- read.xlsx(FILE_PATH, header = TRUE, sheetIndex = sheet_idx)

  stopifnot(
    "Number of rows in temp_df does not match EXPECTED_NUMBER_OF_ROWS." =
      nrow(temp_df) == EXPECTED_NUMBER_OF_ROWS,
    "Number of columns in temp_df does not match EXPECTED_NUMBER_OF_ROWS." =
      ncol(temp_df) == EXPECTED_NUMBER_OF_COLS,
    "df_lst[[df_count]] colnames not identical to REQUIRED_COLUMNS." =
      identical(REQUIRED_COLUMNS, colnames(temp_df))
  )

  temp_df <- temp_df %>%
    mutate(
      Intensity = Intensity - Intensity[2],
      Experiment = paste0("Exp_", df_count)
    ) %>%
    slice(-2) %>%
    mutate(Intensity = (Intensity / Intensity[Input == "yes"]) * 0.5) %>%
    filter(Input == "no")

  df_lst[[df_count]] <- temp_df


} # end read-data for loop

loading_df <- do.call(rbind, df_lst)
loading_df <- loading_df[, names(loading_df) != COLUMN_TO_REMOVE]

message("loading_df preparation complete...")

loading_df <- loading_df %>%
  mutate(Label = case_when(
    ORC == "WT" & Suppressor == "None" ~ "WT",
    ORC == "RA" & Suppressor == "None" ~ "ORC4R",
    ORC == "RA" & Suppressor == "4PS" ~ "+4sofr",
    ORC == "RA" & Suppressor == "1EK" ~ "+1sofr",
    ORC == "RA" & Suppressor == "3PL" ~ "+3sofr"
  ))

for (col_name in names(factor_order)) {

  loading_df[[col_name]] <- factor(
    loading_df[[col_name]],
    levels = factor_order[[col_name]],
    ordered = TRUE
  )
}

#mutate(
#  Suppressor = factor(Suppressor, levels = suppressor_levels),
#  kGlut = factor(kGlut, levels = kGlut_levels),
#  Label = factor(Label, levels = label_levels)
#)


#COLUMNS_TO_IGNORE <- c(
#  "Intensity", "Lane", "Input"
#)
#lapply(colnames(loading_df), function(column_name){
#  if (!column_name %in% COLUMNS_TO_IGNORE){
#    unique(loading_df[, column_name])
#
#  }
#})

#sapply(loading_df, function(x) if(!is.numeric(x)) unique(x))

# Calculate mean and sd for each kGlut
df_summary <- loading_df %>% 
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

pdf(
  file = "~/output.pdf",
  width = 10, height = 8,
  bg = "white",
  compress = TRUE,
  colormodel = "srgb", useDingbats = FALSE
)

ggplot(df_summary, aes(x = kGlut, y = pmol_MCM, fill = Label)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(aes(ymin = pmol_MCM - sd, ymax = pmol_MCM + sd),
                position = position_dodge(width = 0.9), width = 0.2) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "kGlut", y = "MCM (pmol)", title = "Intensity vs kGlut", fill = "Protein") +
  theme_classic()
dev.off()

pdf(
  file = "~/output_1.pdf",
  width = 10, height = 8,
  bg = "white",
  compress = TRUE,
  colormodel = "srgb", useDingbats = FALSE
)
ggplot(df_summary, aes(x = Label, y = pmol_MCM, fill = Label)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(aes(ymin = pmol_MCM - sd, ymax = pmol_MCM + sd),
                position = position_dodge(width = 0.9), width = 0.2) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Sample", y = "MCM (pmol)", title = "Intensity by Sample", fill = "Protein") +
  theme_classic()
dev.off()



ggplot(df_summary, aes(x = Label, y = pmol_MCM, fill = kGlut)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = pmol_MCM - sd, ymax = pmol_MCM + sd),
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Sample", y = "MCM (pmol)", title = "All Samples by Sample Type and kGlut", fill = "kGlut") +
  theme_classic()
