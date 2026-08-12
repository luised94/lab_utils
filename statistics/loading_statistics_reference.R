# Namespace convention: tidyverse verbs used unqualified via library(tidyverse).
# All other packages namespaced (e.g., nlme::lme).
# Date created: 2026-08-11
# Purpose: a SELF-CONTAINED, ARCHIVAL reference for the lab's house-style
# inferential statistics on blocked replicate data. It depends on NOTHING in the
# gel densitometry repo and reads NO external file: it SIMULATES a small blocked
# dataset (fixed seed, reference pinned to 100) and runs the full machinery on it,
# so the code and its expected shape survive on their own.
#
# WHY THIS IS SEPARATE AND NOT WIRED INTO THE PIPELINE: for the current loading
# and radioactivity experiments the reported quantity is the descriptive summary
# (mean, standard deviation, standard error) plus percent- and fold-difference
# matrices, produced by aggregate_loading_replicates.R. Inferential testing is
# appropriate for other cases (more data, or a different design). This script is
# kept as a future reference for that machinery, in the same house style as the
# deposited analysis scripts, without adding statistics to the pipeline output.
#
# Usage:
#   Rscript loading_statistics_reference.R [output_directory]
# With no argument it prints everything to the console and writes nothing. With
# an output_directory it also writes the result CSVs there.
#
# House-style contract reproduced here (identical to the deposited scripts):
#   - The REPORTED test is a paired t-test blocked on replicate.
#   - A paired Wilcoxon signed-rank test is a DECORATIVE companion only (at n=3
#     its two-sided p floors near 0.25 and can never reach significance).
#   - Holm is applied WITHIN one explicitly named family (the WT-contrasts within
#     a single condition); non-WT contrasts are reported with raw p OUTSIDE it.
#   - EXACT p-values, never stars.
#   - A mixed model (random intercept per replicate) is a PRE-COMMITTED
#     consistency DIAGNOSTIC only; it is wrapped so non-convergence flags rather
#     than crashes, and it never licenses swapping in a friendlier p-value.
#   - n=3 CAVEAT: normality is untestable at n=3; the parametric paired test is
#     justified by the paired/blocked design, not by a Shapiro pass. The
#     reference is pinned to 100 with zero within-group variance, so
#     variance-homogeneity tests are undefined and per-group variances are
#     reported directly.

library(tidyverse)

# ==============================================================================
# Configuration
# ==============================================================================
RANDOM_SEED <- 20260811L
REPLICATE_COUNT <- 3L
# The reference condition is pinned to 100 in every replicate (zero variance),
# exactly as the deposited scripts normalize WT = 100 per replicate. The other
# conditions are drawn around plausible per-replicate means; these numbers are
# SYNTHETIC and exist only to exercise the machinery.
REFERENCE_CONDITION_LABEL <- "WT"
EFFECT_SIZE_BASELINE_LABEL <- "ORC4R"   # the mutant baseline effect sizes measure against
SIMULATED_CONDITION_MEANS <- c(
    "WT" = 100, "ORC4R" = 35, "+1sofa" = 78, "+3sofa" = 22, "+4sofa" = 61
)
SIMULATED_WITHIN_CONDITION_SD <- 8      # spread of non-reference conditions across replicates

output_directory_argument <- commandArgs(trailingOnly = TRUE)
output_directory <- if (length(output_directory_argument) >= 1) output_directory_argument[1] else NA_character_
if (!is.na(output_directory) && !dir.exists(output_directory)) {
    dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
}

# ==============================================================================
# Simulate a blocked dataset: every condition present in every replicate, so the
# blocks are complete and every paired contrast uses the full n. The reference is
# set to exactly 100 in each replicate (zero within-group variance by construction).
# ==============================================================================
set.seed(RANDOM_SEED)
condition_labels <- names(SIMULATED_CONDITION_MEANS)

simulated_rows_label <- character(0)
simulated_rows_replicate <- integer(0)
simulated_rows_value <- numeric(0)
for (condition_label in condition_labels) {
    for (replicate_number in seq_len(REPLICATE_COUNT)) {
        simulated_value <- if (condition_label == REFERENCE_CONDITION_LABEL) {
            100
        } else {
            rnorm(1, mean = SIMULATED_CONDITION_MEANS[[condition_label]],
                  sd = SIMULATED_WITHIN_CONDITION_SD)
        }
        simulated_rows_label <- c(simulated_rows_label, condition_label)
        simulated_rows_replicate <- c(simulated_rows_replicate, replicate_number)
        simulated_rows_value <- c(simulated_rows_value, simulated_value)
    }
}
blocked_data <- data.frame(
    label = factor(simulated_rows_label, levels = condition_labels),
    replicate = simulated_rows_replicate,
    percent_reference = simulated_rows_value,
    stringsAsFactors = FALSE
)
message("Simulated blocked data: ", nrow(blocked_data), " rows (",
        length(condition_labels), " conditions x ", REPLICATE_COUNT, " replicates).")
print(blocked_data)

# Structural guards, so the machinery below can assume complete blocks.
label_replicate_counts <- blocked_data %>% count(label)
stopifnot(
    "Each condition must appear once per replicate (complete blocks)." =
        all(label_replicate_counts$n == REPLICATE_COUNT),
    "Reference condition must be pinned to 100 in every replicate." =
        all(blocked_data$percent_reference[blocked_data$label == REFERENCE_CONDITION_LABEL] == 100)
)

# ==============================================================================
# Assumption diagnostics (transparency only; see n=3 caveat in the header)
# ==============================================================================
diagnostics_per_label <- blocked_data %>%
    group_by(label) %>%
    summarise(
        n_replicates = n(),
        mean_percent_reference = mean(percent_reference),
        sd_percent_reference = sd(percent_reference),
        variance_percent_reference = var(percent_reference),
        .groups = "drop"
    )
within_label_residuals <- blocked_data %>%
    group_by(label) %>%
    mutate(within_label_residual = percent_reference - mean(percent_reference)) %>%
    ungroup()
pooled_residual_shapiro <- shapiro.test(within_label_residuals$within_label_residual)
diagnostics_normality <- data.frame(
    diagnostic = "shapiro_wilk_pooled_within_label_residuals",
    shapiro_w = unname(pooled_residual_shapiro$statistic),
    shapiro_p_value = pooled_residual_shapiro$p.value,
    n_residuals = nrow(within_label_residuals),
    caveat = "n=3 per group: normality untestable; reported for transparency only.",
    stringsAsFactors = FALSE
)
message("Assumption diagnostics computed (per-label variances + pooled-residual Shapiro).")

# ==============================================================================
# WT-contrasts: paired t-test blocked on replicate is the reported result.
# Decorative paired Wilcoxon companion. Holm within this one family.
# ==============================================================================
reference_values_by_replicate <- blocked_data %>%
    filter(label == REFERENCE_CONDITION_LABEL) %>%
    arrange(replicate) %>%
    pull(percent_reference)

wt_contrast_labels <- setdiff(condition_labels, REFERENCE_CONDITION_LABEL)

wt_contrast_label <- character(0)
wt_contrast_n_pairs <- integer(0)
wt_contrast_estimate_reference_minus_group <- numeric(0)
wt_contrast_ci_lower <- numeric(0)
wt_contrast_ci_upper <- numeric(0)
wt_contrast_t_statistic <- numeric(0)
wt_contrast_df <- numeric(0)
wt_contrast_raw_p_value <- numeric(0)
wt_contrast_wilcoxon_p_value <- numeric(0)
for (comparison_label in wt_contrast_labels) {
    group_values_by_replicate <- blocked_data %>%
        filter(label == comparison_label) %>%
        arrange(replicate) %>%
        pull(percent_reference)
    stopifnot(
        "Comparison group lacks a full set of paired values (incomplete block)." =
            length(group_values_by_replicate) == REPLICATE_COUNT
    )
    paired_t <- t.test(reference_values_by_replicate, group_values_by_replicate, paired = TRUE)
    paired_wilcoxon <- suppressWarnings(
        wilcox.test(reference_values_by_replicate, group_values_by_replicate, paired = TRUE)
    )
    wt_contrast_label <- c(wt_contrast_label, comparison_label)
    wt_contrast_n_pairs <- c(wt_contrast_n_pairs, length(group_values_by_replicate))
    wt_contrast_estimate_reference_minus_group <- c(wt_contrast_estimate_reference_minus_group, unname(paired_t$estimate))
    wt_contrast_ci_lower <- c(wt_contrast_ci_lower, paired_t$conf.int[1])
    wt_contrast_ci_upper <- c(wt_contrast_ci_upper, paired_t$conf.int[2])
    wt_contrast_t_statistic <- c(wt_contrast_t_statistic, unname(paired_t$statistic))
    wt_contrast_df <- c(wt_contrast_df, unname(paired_t$parameter))
    wt_contrast_raw_p_value <- c(wt_contrast_raw_p_value, paired_t$p.value)
    wt_contrast_wilcoxon_p_value <- c(wt_contrast_wilcoxon_p_value, paired_wilcoxon$p.value)
}
wt_contrast_holm_p_value <- p.adjust(wt_contrast_raw_p_value, method = "holm")
wt_contrast_results <- data.frame(
    comparison = paste0(REFERENCE_CONDITION_LABEL, "_vs_", wt_contrast_label),
    group_label = wt_contrast_label,
    test = "paired t-test, block = replicate",
    holm_family = paste0(length(wt_contrast_labels), " WT-contrasts (single condition)"),
    n_pairs = wt_contrast_n_pairs,
    estimate_reference_minus_group = wt_contrast_estimate_reference_minus_group,
    ci_lower = wt_contrast_ci_lower,
    ci_upper = wt_contrast_ci_upper,
    t_statistic = wt_contrast_t_statistic,
    df = wt_contrast_df,
    raw_p_value = wt_contrast_raw_p_value,
    holm_adjusted_p_value = wt_contrast_holm_p_value,
    wilcoxon_p_value_decorative = wt_contrast_wilcoxon_p_value,
    stringsAsFactors = FALSE
)
message("WT-contrast paired t-tests + decorative Wilcoxon + Holm computed.")
print(wt_contrast_results)

# ==============================================================================
# Effect sizes: (group - baseline) paired estimate with 95% CI. Magnitude and
# direction only; at n=3 the CIs are wide (df = 2). Reported OUTSIDE the Holm
# family (this is not a WT-contrast).
# ==============================================================================
baseline_values_by_replicate <- blocked_data %>%
    filter(label == EFFECT_SIZE_BASELINE_LABEL) %>%
    arrange(replicate) %>%
    pull(percent_reference)
effect_size_labels <- setdiff(condition_labels, c(REFERENCE_CONDITION_LABEL, EFFECT_SIZE_BASELINE_LABEL))
effect_size_label <- character(0)
effect_size_estimate_group_minus_baseline <- numeric(0)
effect_size_ci_lower <- numeric(0)
effect_size_ci_upper <- numeric(0)
for (group_label in effect_size_labels) {
    group_values_by_replicate <- blocked_data %>%
        filter(label == group_label) %>%
        arrange(replicate) %>%
        pull(percent_reference)
    effect_t <- t.test(group_values_by_replicate, baseline_values_by_replicate, paired = TRUE)
    effect_size_label <- c(effect_size_label, group_label)
    effect_size_estimate_group_minus_baseline <- c(effect_size_estimate_group_minus_baseline, unname(effect_t$estimate))
    effect_size_ci_lower <- c(effect_size_ci_lower, effect_t$conf.int[1])
    effect_size_ci_upper <- c(effect_size_ci_upper, effect_t$conf.int[2])
}
effect_size_results <- data.frame(
    group_label = effect_size_label,
    quantity = paste0("ordinary paired effect size (group - ", EFFECT_SIZE_BASELINE_LABEL, ")"),
    estimate_group_minus_baseline = effect_size_estimate_group_minus_baseline,
    ci_lower = effect_size_ci_lower,
    ci_upper = effect_size_ci_upper,
    stringsAsFactors = FALSE
)
message("Effect sizes vs baseline computed.")
print(effect_size_results)

# ==============================================================================
# Pre-committed mixed-model consistency diagnostic (DIAGNOSTIC ONLY).
# random intercept per replicate; the reference is the factor base level, so each
# label coefficient is (group - reference). Wrapped so non-convergence FLAGS the
# diagnostic rather than crashing; a non-converging diagnostic is itself a flag.
# ==============================================================================
mixed_model_fit <- tryCatch(
    nlme::lme(percent_reference ~ label, random = ~ 1 | replicate, data = blocked_data),
    error = function(e) e
)
if (inherits(mixed_model_fit, "error")) {
    mixed_model_results <- data.frame(
        term = NA_character_, estimate_group_minus_reference = NA_real_,
        std_error = NA_real_, t_value = NA_real_, model_p_value = NA_real_,
        converged = FALSE,
        note = paste0("lme did not converge: ", conditionMessage(mixed_model_fit),
                      " -- the paired t-test stands as reported."),
        stringsAsFactors = FALSE
    )
    message("Mixed-model diagnostic DID NOT CONVERGE; flagged. Paired t-test stands.")
} else {
    mixed_model_fixed_effects <- summary(mixed_model_fit)$tTable
    label_coefficient_rows <- rownames(mixed_model_fixed_effects) != "(Intercept)"
    mixed_model_results <- data.frame(
        term = rownames(mixed_model_fixed_effects)[label_coefficient_rows],
        estimate_group_minus_reference = mixed_model_fixed_effects[label_coefficient_rows, "Value"],
        std_error = mixed_model_fixed_effects[label_coefficient_rows, "Std.Error"],
        t_value = mixed_model_fixed_effects[label_coefficient_rows, "t-value"],
        model_p_value = mixed_model_fixed_effects[label_coefficient_rows, "p-value"],
        converged = TRUE,
        note = "DIAGNOSTIC ONLY; reported result is the paired t-test.",
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    message("Mixed-model diagnostic converged; flagged as consistency check.")
}
print(mixed_model_results)

# ==============================================================================
# Optional CSV output
# ==============================================================================
if (!is.na(output_directory)) {
    write.csv(blocked_data, file.path(output_directory, "statistics_reference_simulated_data.csv"), row.names = FALSE)
    write.csv(diagnostics_per_label, file.path(output_directory, "statistics_reference_diagnostics_per_label.csv"), row.names = FALSE)
    write.csv(diagnostics_normality, file.path(output_directory, "statistics_reference_diagnostics_normality.csv"), row.names = FALSE)
    write.csv(wt_contrast_results, file.path(output_directory, "statistics_reference_wt_contrasts.csv"), row.names = FALSE)
    write.csv(effect_size_results, file.path(output_directory, "statistics_reference_effect_sizes.csv"), row.names = FALSE)
    write.csv(mixed_model_results, file.path(output_directory, "statistics_reference_mixed_model_diagnostic.csv"), row.names = FALSE)
    message("Wrote result CSVs to: ", output_directory)
}
message("Statistics reference script complete.")
