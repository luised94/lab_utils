library(xlsx)
library(tidyverse)
file_path = "C:\\Users\\liusm\\Dropbox (MIT)\\Lab\\Experiments\\ATP hydrolysis by Orc1_Orc4\\Assays\\Loading\\2022_12_18 Loading Assays Repeats for publication\\Analysis.xlsx" 
#file_path = "C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Experiments\\ATP hydrolysis by Orc1_Orc4\\Assays\\Loading\\2022_12_18 Loading Assays Repeats for publication\\Analysis.xlsx" 
# process_sheet <- function(df) {
#   df %>%
#     mutate(Intensity = Intensity - df$Intensity[2]) %>%
#     slice(-2) %>%
#     mutate(Intensity = (Intensity/Intensity[Input == "yes"]) * 0.5) %>% 
#     filter(Input == "no")
# }
# df1 <- read.xlsx(file_path, sheet = 1) %>% process_sheet() 
# df2 <- read.xlsx(file_path, sheet = 2) %>% process_sheet()
# df3 <- read.xlsx(file_path, sheet = 3) %>% process_sheet()
df1 <- read.xlsx(file_path, header = TRUE, sheetIndex = 1)
df2 <- read.xlsx(file_path, header = TRUE, sheetIndex =  2) 
df3 <- read.xlsx(file_path, header = TRUE, sheetIndex =  3)

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

#Make bigger to save higher res file or modify the parameters in the function.
ggsave("loading-kgluttitr-wt,4r,4rps.svg", plot = last_plot(), device = "svg")

#Post processing in illustrator