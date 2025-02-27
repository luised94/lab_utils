library(flowCore)
library(ggcyto)
library(ggridges)
library(tidyr)
#Read in file, default is to linearize
#Add sample(100:500,50) to which.lines argument to extract a sample of the FCS.
fcs_file <- "C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Experiments\\ATP hydrolysis by Orc1_Orc4\\Assays\\Yeast Genetics\\2022_12_14 CHIP-seq_AID2_ORC-MCM-supps\\Exp_20230522_AID2ALFA_G2intoG15Ph2_timecourse\\01-Tube-E5.fcs"
fcs_dir <- "C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Experiments\\ATP hydrolysis by Orc1_Orc4\\Assays\\Yeast Genetics\\2022_12_14 CHIP-seq_AID2_ORC-MCM-supps\\Exp_20230712_bar1,cdc152,wt_controls\\"

fitc_values <- lapply(list.files(fcs_dir, pattern = ".fcs", full.names = TRUE), function(fcs_file) {
  flow_exp <- read.FCS(fcs_file,
                       transformation = "linearize",
                       alter.names = TRUE)
  fcs_label <- gsub(".fcs", "", basename(fcs_file))  # Extract label from file name
  data.frame(FITC_Value = na.omit(exprs(flow_exp)[,"FL1.A"]), FCS_File = fcs_label)
})

sample_df <- do.call(rbind, fitc_values)

sample_df$FCS_File <- gsub("-", "_", sample_df$FCS_File) 

sample_df <- sample_df[!grepl("sytox", sample_df$FCS_File), ]

data_names <- c("Plate","Strain","Cell_Cycle_Stage","Time", "Well_Num")
sample_df <- separate(data = sample_df, col = FCS_File, into = data_names, sep = "_")

wrongly_labeled <- is.na(sample_df$Well_Num)

sample_df$Well_Num[is.na(sample_df$Well_Num)] <- sample_df$Time[is.na(sample_df$Well_Num)]
sample_df$Time[wrongly_labeled] <- NA


ggplot(sample_df, aes(x = FITC_Value, y = Time)) + 
  geom_density_ridges() + 
  coord_cartesian(xlim =  c(0, 200000)) + 
  xlim(0, 100000) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = expand_scale(mult = c(0.01, .7))) + 
  facet_wrap(~Strain) +
  theme_ridges()


ggplot(plot_data, aes(x = Concentrations, y = Value, color = Experiment)) +
  geom_point() +
  geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD), width = 0.1) +
  geom_line(data = fit_results, aes(group = Experiment, y = predict)) +
  labs(x = "Concentrations", y = "Value") +
  scale_color_discrete(name = "Experiment") +
  theme_minimal()

# 
# library(ggplot2)
# library(ggridges)
# 
# data <- data.frame(x = 1:5, y = rep(1, 5), height = c(0, 1, 3, 4, 2))
# ggplot(data, aes(x, y, height = height)) + geom_ridgeline()
# 
# 
# d <- data.frame(
#   x = rep(1:5, 3),
#   y = c(rep(0, 5), rep(1, 5), rep(2, 5)),
#   height = c(0, 1, 3, 4, 0, 1, 2, 3, 5, 4, 0, 5, 4, 4, 1)
# )
# 
# ggplot(d, aes(x, y, height = height, group = y)) + 
#   geom_ridgeline(fill = "lightblue")
# 
# ggplot(d, aes(x, y, height = height, group = y)) + 
#   geom_density_ridges(stat = "identity", scale = 1)
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species)) + geom_density_ridges()
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species)) + geom_density_ridges2()
# 
# # modified dataset that represents species as a number
# iris_num <- transform(iris, Species_num = as.numeric(Species))
# 
# # does not work, causes error
# # ggplot(iris_num, aes(x = Sepal.Length, y = Species)) + geom_density_ridges()
# 
# # works 
# ggplot(iris_num, aes(x = Sepal.Length, y = Species_num, group = Species_num)) + 
#   geom_density_ridges()
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species)) + 
#   geom_density_ridges(rel_min_height = 0.01)
# 
# # scale = 0.9, not quite touching
# ggplot(iris, aes(x = Sepal.Length, y = Species)) + geom_density_ridges(scale = 0.9)
# 
# # scale = 1, exactly touching
# ggplot(iris, aes(x = Sepal.Length, y = Species)) + geom_density_ridges(scale = 1)
# 
# # scale = 5, substantial overlap
# ggplot(iris, aes(x = Sepal.Length, y = Species)) + geom_density_ridges(scale = 5)
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species)) + 
#   geom_density_ridges(scale = 1) + facet_wrap(~Species)
# 
# d <- data.frame(
#   x = rep(1:5, 3) + c(rep(0, 5), rep(0.3, 5), rep(0.6, 5)),
#   y = c(rep(0, 5), rep(1, 5), rep(3, 5)),
#   height = c(0, 1, 3, 4, 0, 1, 2, 3, 5, 4, 0, 5, 4, 4, 1))
# 
# ggplot(d, aes(x, y, height = height, group = y, fill = factor(x+y))) +
#   geom_ridgeline_gradient() +
#   scale_fill_viridis_d(direction = -1, guide = "none")
# 
# ggplot(lincoln_weather, aes(x = `Mean Temperature [F]`, y = Month, fill = stat(x))) +
#   geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
#   scale_fill_viridis_c(name = "Temp. [F]", option = "C") +
#   labs(title = 'Temperatures in Lincoln NE in 2016')
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species)) +
#   stat_density_ridges(quantile_lines = TRUE)
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species)) +
#   stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7)
# 
# ggplot(iris, aes(x=Sepal.Length, y=Species, fill = factor(stat(quantile)))) +
#   stat_density_ridges(
#     geom = "density_ridges_gradient", calc_ecdf = TRUE,
#     quantiles = 4, quantile_lines = TRUE
#   ) +
#   scale_fill_viridis_d(name = "Quartiles")
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species, fill = factor(stat(quantile)))) +
#   stat_density_ridges(
#     geom = "density_ridges_gradient",
#     calc_ecdf = TRUE,
#     quantiles = c(0.025, 0.975)
#   ) +
#   scale_fill_manual(
#     name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
#     labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
#   )
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
#   scale_fill_viridis_c(name = "Tail probability", direction = -1)
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species)) +
#   geom_density_ridges(jittered_points = TRUE)
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species)) +
#   geom_density_ridges(
#     jittered_points = TRUE, position = "raincloud",
#     alpha = 0.7, scale = 0.9
#   )
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species)) +
#   geom_density_ridges(
#     jittered_points = TRUE,
#     position = position_points_jitter(width = 0.05, height = 0),
#     point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
#   )
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species, fill = Species)) +
#   geom_density_ridges(
#     aes(point_color = Species, point_fill = Species, point_shape = Species),
#     alpha = .2, point_alpha = 1, jittered_points = TRUE
#   ) +
#   scale_point_color_hue(l = 40) +
#   scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23))
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species, fill = Species)) +
#   geom_density_ridges(
#     aes(point_shape = Species, point_fill = Species, point_size = Petal.Length), 
#     alpha = .2, point_alpha = 1, jittered_points = TRUE
#   ) +
#   scale_point_color_hue(l = 40) + scale_point_size_continuous(range = c(0.5, 4)) +
#   scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23))
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species)) +
#   geom_density_ridges(
#     jittered_points = TRUE, quantile_lines = TRUE, scale = 0.9, alpha = 0.7,
#     vline_size = 1, vline_color = "red",
#     point_size = 0.4, point_alpha = 1,
#     position = position_raincloud(adjust_vlines = TRUE)
#   )
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species, height = stat(density))) + 
#   geom_density_ridges(stat = "density")
# 
