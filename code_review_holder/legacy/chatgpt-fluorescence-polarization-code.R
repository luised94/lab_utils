# Install required packages if not already installed
install.packages(c("dplyr", "broom", "ggplot2"))

# Load required packages
library(dplyr)
library(broom)
library(ggplot2)

# Create a dataframe with the actual data
data <- data.frame(
  Exp1_Rep1 = c(319.9013, 226.0579, 216.8701, 214.6696, 215.8831, 216.8711, 215.6911, 183.9771),
  Exp1_Rep2 = c(315.187, 223.1486, 216.6287, 223.8176, 217.4015, 204.0279, 219.9958, 182.8383),
  Exp1_Rep3 = c(313.66662, 231.20887, 222.79615, 226.10572, 99.02211, 203.59574, 212.47239, 189.44633),
  Exp2_Rep1 = c(375.203, 176.0811, 278.0506, 276.1944, 284.5559, 267.1035, 269.9042, 247.5806),
  Exp2_Rep2 = c(280.9988, 286.316, 274.2547, 279.1131, 280.5905, 258.6079, 267.6636, 244.8262),
  Exp2_Rep3 = c(158.3907, 288.1005, 275.8863, 278.5713, 277.6004, 267.401, 264.9308, 242.9766),
  Concentrations = c(128, 64, 32, 16, 8, 4, 2, 1)
)

# Define the Hill equation function
hill_equation <- function(x, Vmax, Kd, n) {
  Vmax * (x^n) / (Kd^n + x^n)
}

# Fit Hill equation to each set of experiments
fit_results <- data %>%
  summarise_all(list(~ list(nls(
    . ~ hill_equation(Concentrations, Vmax, Kd, n),
    start = list(Vmax = 300, Kd = 50, n = 2)
  )))) %>%
  tidyr::unnest(everything()) %>%
  broom::tidy()

# Create a new dataframe with the average and standard deviation for each experiment
plot_data <- data %>%
  summarise_all(list(
    Avg = mean,
    SD = sd
  )) %>%
  tidyr::gather(Experiment, Value, -Concentrations)

# Plot the data
ggplot(plot_data, aes(x = Concentrations, y = Value, color = Experiment)) +
  geom_point() +
  geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD), width = 0.1) +
  geom_line(data = fit_results, aes(group = Experiment, y = predict)) +
  labs(x = "Concentrations", y = "Value") +
  scale_color_discrete(name = "Experiment") +
  theme_minimal()
