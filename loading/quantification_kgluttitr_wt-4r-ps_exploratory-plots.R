# Date created: 2026-04-01
# Exploratory plots extracted from quantification_kgluttitr_wt-4r-ps.R
# These are alternative visualizations of the multi-salt loading data.
# The primary plot (faceted_by_kglut_plot) remains in the main script.
# Usage: source("quantification_kgluttitr_wt-4r-ps_exploratory-plots.R")
# Requires: quantification_kgluttitr_wt-4r-ps.R (sources it for data frames)

source("quantification_kgluttitr_wt-4r-ps.R")

# ==============================================================================
# Derived data for exploratory plots
# ==============================================================================

# Fold change relative to WT with propagated error
df_normalized <- df_summary %>%
  group_by(kGlut) %>%
  mutate(
    wt_val = Rel_to_Input_mean[Label == "WT"][1],
    wt_sd = Rel_to_Input_sd[Label == "WT"][1],
    fold_change = Rel_to_Input_mean / wt_val,
    fold_change_sd = fold_change * sqrt((Rel_to_Input_sd/Rel_to_Input_mean)^2 + (wt_sd/wt_val)^2)
  ) %>%
  ungroup()

df_summary_300 <- df_summary %>% filter(kGlut == "300")

# ==============================================================================
# Exploratory plots
# ==============================================================================

intensity_vs_kglut_plot <- ggplot(df_summary, aes(x = kGlut, y = Rel_to_Input_mean, color = Label, fill = Label)) +
  geom_line(linewidth = 1.2, aes(group = Label)) +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_errorbar(aes(ymin = Rel_to_Input_mean - Rel_to_Input_sd, ymax = Rel_to_Input_mean + Rel_to_Input_sd), 
                width = 0.15, linewidth = 0.6) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "kGlut Concentration (mM)", 
       y = "MCM (pmol)", 
       title = "Intensity vs kGlut",
       color = "Protein",
       fill = "Protein") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
  )

all_samples_dodged_plot <- ggplot(df_summary, aes(x = Label, y = Rel_to_Input_mean, fill = Label, group = kGlut)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, 
           color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Rel_to_Input_mean - Rel_to_Input_sd, ymax = Rel_to_Input_mean + Rel_to_Input_sd),
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.6) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Sample Type", 
       y = "MCM (pmol)", 
       title = "MCM Levels by Sample Type and kGlut Concentration", 
       fill = "Sample") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(face = "bold")
  )

kglut_300_only_plot <- ggplot(df_summary_300, aes(x = Label, y = Rel_to_Input_mean, fill = Label)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Rel_to_Input_mean - Rel_to_Input_sd, ymax = Rel_to_Input_mean + Rel_to_Input_sd), 
                width = 0.25, linewidth = 0.6) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Sample Type", 
       y = "MCM (pmol)", 
       title = "300 mM kGlut Only", 
       fill = "Sample") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(face = "bold")
  )

normalized_to_wt_plot <- ggplot(df_normalized, aes(x = kGlut, y = fold_change, color = Label, fill = Label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_line(linewidth = 1.2, aes(group = Label)) +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_errorbar(aes(ymin = fold_change - fold_change_sd, 
                    ymax = fold_change + fold_change_sd), 
                width = 0.15, linewidth = 0.6) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "kGlut Concentration (mM)", 
       y = "Fold Change (relative to WT)", 
       title = "Normalized to Wildtype",
       color = "Sample",
       fill = "Sample") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
  )

faceted_by_label_plot <- ggplot(df_summary, aes(x = kGlut, y = Rel_to_Input_mean, fill = kGlut)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Rel_to_Input_mean - Rel_to_Input_sd, ymax = Rel_to_Input_mean + Rel_to_Input_sd), 
                width = 0.25, linewidth = 0.6) +
  facet_wrap(~Label, nrow = 1) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "kGlut Concentration (mM)", 
       y = "MCM (pmol)", 
       title = "kGlut Dose-Response by Sample Type",
       fill = "kGlut") +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  )

# ==============================================================================
# Save exploratory plots
# ==============================================================================
exploratory_plot_names <- c(
    "intensity_vs_kglut_plot",
    "all_samples_dodged_plot",
    "kglut_300_only_plot",
    "normalized_to_wt_plot",
    "faceted_by_label_plot"
)

for (plot_name in exploratory_plot_names) {
    plot_filepath <- file.path(OUTPUT_DIRECTORY, paste0(plot_name, ".pdf"))
    if (!file.exists(plot_filepath) || OVERWRITE_PLOTS) {
        ggsave(plot_filepath, get(plot_name), width = 8, height = 5)
        cat("Saved plot:", basename(plot_filepath), "\n")
    } else {
        cat("Skipped plot (already exists):", basename(plot_filepath), "\n")
    }
}

message("Exploratory plots script complete.")
