#!/usr/bin/env Rscript
# ==============================================================================
# Task A: standalone outlier diagnostic for the WT plusATP suspect point.
#
# Suspect: WT plusATP bound_fraction = 0.298086 on gel wt-1-3-4 replicate 2
# (ratio CSV 23.6-27.4mm), against ~0.79-0.96 for the other five WT plusATP
# replicates.
#
# Metadata branch (BRANCH_metadata_findings.md) settled the acquisition question:
#   Same PMT (655 V) as its comparator wt-1-3-4 rep 3, no saturation
#   (RangeHigh 21139 of 65535), Linear scale -- the low value is NOT a scan
#   artifact. The ~10x count-scale difference is a gain/exposure difference that
#   cancels in the within-lane ratio bound/(bound+free). Cause unexplained per
#   lab notes (same DNA added, assay run the day before, nothing recorded).
#   Therefore this script does NOT model saturation; it decides whether to flag.
#
# This script is READ-ONLY over the data. It does not remove anything. The
# exclusion decision is made by a human and applied in the plotting script via a
# manually-set gel-id variable (see report section 7).
#
# Dependencies: base R only (stats). No tidyverse, no car, no outliers. Grubbs
# and Levene are inlined exactly (short, no package needed). If you prefer the
# packaged forms, outliers::grubbs.test and car::leveneTest give the same
# numbers; they are intentionally not used so this runs on a bare R install and
# adds no dependency that can rot.
#
# Determinism: no RNG, no seed-dependent step. Rows are sorted by a stable key
# (screen, replicate, genotype, atp) before anything downstream. Divisors and
# comparators are looked up by derived (gel, genotype, atp), never by lane
# position -- lane order flips between gels.
#
# Run: Rscript task_a_outlier_diagnostic.R
# Writes: task_a_outlier_diagnostic.txt (full report) and
#         task_a_qq_residuals.pdf (Q-Q of model residuals, with/without suspect)
#         in the working directory. Also prints the report to stdout.
# ==============================================================================

# ------------------------------------------------------------------------------
# Inlined inputs. Two manifests, three gels each. Paths are the WSL-normalized
# form from the manifest CSVs. Edit here if the files move; there are no
# command-line arguments by design (invoked as a bare Rscript call).
#
# GEL_KEY is the stable gel id used everywhere downstream and is the value the
# plotting script's exclusion variable must match. It is screen + replicate, not
# a lane index and not a file path, so it survives a manifest row moving.
# ------------------------------------------------------------------------------
GEL_MANIFEST <- data.frame(
    screen = c("wt-1-3-4", "wt-1-3-4", "wt-1-3-4",
               "4r-5-6",   "4r-5-6",   "4r-5-6"),
    replicate = c(1L, 2L, 3L, 1L, 2L, 3L),
    ratio_csv = c(
        "/mnt/c/Users/Luised94/Desktop/lab/experiments/20260716_LM-0009_gs_1-3-4-repeats/20210908_LM-0009_s0003_rotated_gelshift_orc1,3,4-suppressor-effect-and-ATP-dependence_Phosphor_gel_analysis/gel_shift_ratio_16.8-23.6mm_bound_over_46.6-56mm_free.csv",
        "/mnt/c/Users/Luised94/Desktop/lab/experiments/20260716_LM-0009_gs_1-3-4-repeats/20260723_LM-0009_s0001_rotated_gelshift_wt,1,3,4-sofa-repeat-2-1000_Phosphor_gel_analysis/gel_shift_ratio_23.6-27.4mm_bound_over_38.9-45.6mm_free.csv",
        "/mnt/c/Users/Luised94/Desktop/lab/experiments/20260716_LM-0009_gs_1-3-4-repeats/20260724_LM-0009_s0002_rotated_gelshift_wt,1,3,4-sofa-repeat-3-1000_Phosphor_gel_analysis/gel_shift_ratio_23.8-29.4mm_bound_over_39.8-45.6mm_free.csv",
        "/mnt/c/Users/Luised94/Desktop/lab/experiments/20260803_LM-0013_gs_ORC5-6-sofa-data-for-analysis/20220303_rotated_s0001_wt,4r,5,6-sofa-repeat-1-1000-Phosphor_gel_analysis/gel_shift_ratio_16.8-21.4mm_bound_over_42.8-52.4mm_free.csv",
        "/mnt/c/Users/Luised94/Desktop/lab/experiments/20260803_LM-0013_gs_ORC5-6-sofa-data-for-analysis/20220304_rotated_s0002_wt,4r,5,6-sofa-repeat-2-1000-Phosphor_gel_analysis/gel_shift_ratio_15.6-22.2mm_bound_over_43.2-54.6mm_free.csv",
        "/mnt/c/Users/Luised94/Desktop/lab/experiments/20260803_LM-0013_gs_ORC5-6-sofa-data-for-analysis/20220307_rotated_s0003_wt,4r,5,6-sofa-repeat-3-1000-Phosphor_gel_analysis/gel_shift_ratio_19-23.8mm_bound_over_44.2-53.4mm_free.csv"),
    stringsAsFactors = FALSE)
GEL_MANIFEST$gel_key <- paste0(GEL_MANIFEST$screen, " rep ", GEL_MANIFEST$replicate)

# The suspect, stated as data so the report can key on it rather than a literal
# appearing mid-code. SUSPECT_TOLERANCE guards the assertion that we actually
# found the row we think we did (see the loud check after ingestion).
SUSPECT_GEL_KEY <- "wt-1-3-4 rep 2"
SUSPECT_GENOTYPE <- "WT"
SUSPECT_ATP <- "plusATP"
SUSPECT_EXPECTED_FRACTION <- 0.298086
SUSPECT_TOLERANCE <- 1e-4

# ------------------------------------------------------------------------------
# Local-run override. When executed outside the Windows/WSL tree (e.g. for
# testing against copies), set TASK_A_LOCAL_DIR to a directory holding the six
# CSVs under the short names below. Empty string means use the manifest paths
# above. This exists ONLY so the script is executable for verification; on your
# machine leave it "" and the real manifest paths are used.
# ------------------------------------------------------------------------------
TASK_A_LOCAL_DIR <- Sys.getenv("TASK_A_LOCAL_DIR", "")
if (nzchar(TASK_A_LOCAL_DIR)) {
    LOCAL_NAMES <- c("g_wt134_r1.csv", "g_wt134_r2.csv", "g_wt134_r3.csv",
                     "g_4r56_r1.csv", "g_4r56_r2.csv", "g_4r56_r3.csv")
    GEL_MANIFEST$ratio_csv <- file.path(TASK_A_LOCAL_DIR, LOCAL_NAMES)
}

REPORT_LINES <- character(0)
emit <- function(...) {
    # The one unavoidable output helper. Not an abstraction over logic: it only
    # tees a line to stdout and to the report buffer so the .txt and the console
    # never drift. Inlining a two-target write at ~40 call sites would be worse.
    line <- paste0(...)
    cat(line, "\n", sep = "")
    REPORT_LINES[[length(REPORT_LINES) + 1L]] <<- line
    invisible(NULL)
}

emit("================================================================================")
emit("TASK A: WT plusATP OUTLIER DIAGNOSTIC")
emit("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " (local system time)")
emit("================================================================================")
emit("")
emit("Mechanism line (from metadata branch, quoted): Same PMT (655 V) as its")
emit("comparator wt-1-3-4 rep 3, no saturation (RangeHigh 21139 of 65535), Linear")
emit("scale -- the low WT plusATP value is not a scan artifact. Cause unexplained")
emit("per lab notes. Saturation is off the table; this script decides flag vs keep.")
emit("")

# ------------------------------------------------------------------------------
# Ingestion. Base read.csv per gel, tag with gel_key/screen/replicate, rbind.
# The genotype+atp labels in these CSVs live in sample_label (e.g.
# "WT ORC plusATP", "ORC4R 5EK noATP", "No ORC"). We derive two clean columns:
#   genotype: WT | No ORC | the suppressor/4R label
#   atp: plusATP | noATP | NA (No ORC reference has atp but is not an ATP contrast)
# from sample_label + atp_presence, NOT from lane position. Lane order flips
# between gels (verified: on some gels lane 2 is plusATP, on others noATP), so
# anything positional would mispair.
# ------------------------------------------------------------------------------
all_rows_list <- vector("list", nrow(GEL_MANIFEST))
for (gel_index in seq_len(nrow(GEL_MANIFEST))) {
    this_path <- GEL_MANIFEST$ratio_csv[gel_index]
    if (!file.exists(this_path)) {
        stop("Ratio CSV not found for ", GEL_MANIFEST$gel_key[gel_index],
             ": ", this_path,
             "\nIf running outside the lab tree, set TASK_A_LOCAL_DIR.")
    }
    one_gel <- read.csv(this_path, stringsAsFactors = FALSE, check.names = TRUE)
    one_gel$gel_key <- GEL_MANIFEST$gel_key[gel_index]
    one_gel$screen <- GEL_MANIFEST$screen[gel_index]
    one_gel$replicate <- GEL_MANIFEST$replicate[gel_index]
    all_rows_list[[gel_index]] <- one_gel
}
# Column sets differ across gels (different suppressors present), so bind on the
# union of columns rather than assuming identical headers. rbind would error on
# mismatched names; this fills absent columns with NA.
all_column_names <- unique(unlist(lapply(all_rows_list, names)))
for (gel_index in seq_along(all_rows_list)) {
    missing_columns <- setdiff(all_column_names, names(all_rows_list[[gel_index]]))
    for (missing_column in missing_columns) {
        all_rows_list[[gel_index]][[missing_column]] <- NA
    }
    all_rows_list[[gel_index]] <- all_rows_list[[gel_index]][, all_column_names]
}
combined_rows <- do.call(rbind, all_rows_list)
rownames(combined_rows) <- NULL

# Derive genotype and atp from sample_label + atp_presence. These string rules
# reproduce what the plotting script's factor logic encodes; kept explicit here.
sample_label_text <- combined_rows$sample_label
combined_rows$atp <- ifelse(combined_rows$atp_presence == "yes", "plusATP", "noATP")
combined_rows$genotype <- NA_character_
combined_rows$genotype[grepl("^No ORC", sample_label_text)] <- "No ORC"
combined_rows$genotype[grepl("^WT ORC", sample_label_text)] <- "WT"
# Suppressor / 4R analytes: everything that is neither No ORC nor WT. Use the
# sample_label with the ATP words stripped as the genotype token so that e.g.
# "ORC4R 5EK plusATP" and "ORC4R 5EK noATP" share genotype "ORC4R 5EK".
other_mask <- is.na(combined_rows$genotype)
combined_rows$genotype[other_mask] <- trimws(gsub("(plusATP|noATP| ATP)", "",
    gsub("(?i)(plus|no)?ATP", "", sample_label_text[other_mask], perl = TRUE)))
combined_rows$genotype[other_mask] <- trimws(combined_rows$genotype[other_mask])

# For No ORC the atp column is not an ATP contrast; mark it so it never enters
# an ATP-paired analysis, matching the plotting script's treatment of No ORC as
# a pooled baseline rather than an ATP bar.
combined_rows$atp[combined_rows$genotype == "No ORC"] <- NA

# Stable deterministic ordering. Every downstream step reads rows in this order.
combined_rows <- combined_rows[order(combined_rows$screen,
                                     combined_rows$replicate,
                                     combined_rows$genotype,
                                     combined_rows$atp,
                                     method = "radix"), ]
rownames(combined_rows) <- NULL

# Loud check: confirm the suspect row exists and is the value we designed around.
# If a manifest edit or relabel moved it, fail here rather than silently test the
# wrong point.
suspect_mask <- combined_rows$gel_key == SUSPECT_GEL_KEY &
                combined_rows$genotype == SUSPECT_GENOTYPE &
                combined_rows$atp == SUSPECT_ATP
suspect_mask[is.na(suspect_mask)] <- FALSE
if (sum(suspect_mask) != 1L) {
    stop("Expected exactly one suspect row (", SUSPECT_GEL_KEY, " / ",
         SUSPECT_GENOTYPE, " / ", SUSPECT_ATP, "); found ", sum(suspect_mask),
         ". Data labels may have changed.")
}
suspect_fraction_found <- combined_rows$bound_fraction[suspect_mask]
if (abs(suspect_fraction_found - SUSPECT_EXPECTED_FRACTION) > SUSPECT_TOLERANCE) {
    stop("Suspect bound_fraction ", round(suspect_fraction_found, 6),
         " differs from expected ", SUSPECT_EXPECTED_FRACTION,
         " by more than tolerance. Verify you are reading the intended data.")
}
emit("Ingestion OK. Rows: ", nrow(combined_rows),
     ". Suspect row confirmed at bound_fraction = ",
     round(suspect_fraction_found, 6), ".")
emit("")

# The WT plusATP and WT noATP cells, in stable order, used by several sections.
wt_plusatp <- combined_rows[combined_rows$genotype == "WT" &
                            combined_rows$atp == "plusATP", ]
wt_plusatp <- wt_plusatp[order(wt_plusatp$screen, wt_plusatp$replicate), ]
wt_noatp <- combined_rows[combined_rows$genotype == "WT" &
                          combined_rows$atp == "noATP", ]
wt_noatp <- wt_noatp[order(wt_noatp$screen, wt_noatp$replicate), ]

# ==============================================================================
# SECTION 1: MULTIPLICATIVE-vs-SPECIFIC ratio profile on the suspect gel.
#
# Question: is the whole suspect gel uniformly off (multiplicative -- a gel-level
# scaling that per-gel WT normalization would cure), or is only WT plusATP low
# (specific -- an unexplained lane effect that normalization would amplify)?
#
# The metadata branch predicts SPECIFIC: it argues the count-scale difference
# cancels in the ratio, so the gel's other analytes should read at normal
# cross-gel levels and only WT plusATP should stand out. This section tests that
# prediction. A uniform result would CONTRADICT the branch and reopen the
# acquisition question.
#
# Method: for each analyte on the suspect gel, ratio = its bound_fraction /
# mean bound_fraction of the SAME (genotype, atp) on the other gels that carry
# it. A ratio near 1 means "normal for that analyte"; << 1 means depressed.
# ==============================================================================
emit("--------------------------------------------------------------------------------")
emit("SECTION 1: Multiplicative-vs-specific ratio profile on ", SUSPECT_GEL_KEY)
emit("--------------------------------------------------------------------------------")
emit("Per-analyte ratio = this gel's bound_fraction / mean of same analyte on other gels.")
emit(sprintf("  %-22s %-9s %10s %10s %8s", "genotype", "atp", "this_gel", "others_mean", "ratio"))
suspect_gel_rows <- combined_rows[combined_rows$gel_key == SUSPECT_GEL_KEY, ]
suspect_gel_rows <- suspect_gel_rows[order(suspect_gel_rows$genotype,
                                           suspect_gel_rows$atp), ]
section1_ratios <- numeric(0)
section1_is_wt_plusatp <- logical(0)
for (row_index in seq_len(nrow(suspect_gel_rows))) {
    this_genotype <- suspect_gel_rows$genotype[row_index]
    this_atp <- suspect_gel_rows$atp[row_index]
    this_fraction <- suspect_gel_rows$bound_fraction[row_index]
    other_mask <- combined_rows$gel_key != SUSPECT_GEL_KEY &
                  combined_rows$genotype == this_genotype &
                  !is.na(combined_rows$genotype)
    if (is.na(this_atp)) {
        other_mask <- other_mask & is.na(combined_rows$atp)
    } else {
        other_mask <- other_mask & !is.na(combined_rows$atp) &
                      combined_rows$atp == this_atp
    }
    other_mask[is.na(other_mask)] <- FALSE
    if (sum(other_mask) == 0L) {
        # Analyte appears only on the suspect gel (some suppressors are screen
        # specific). No cross-gel comparator exists; report and skip the ratio.
        emit(sprintf("  %-22s %-9s %10.4f %10s %8s",
                     this_genotype, ifelse(is.na(this_atp), "-", this_atp),
                     this_fraction, "none", "n/a"))
        next
    }
    others_mean <- mean(combined_rows$bound_fraction[other_mask])
    this_ratio <- this_fraction / others_mean
    emit(sprintf("  %-22s %-9s %10.4f %10.4f %8.3f",
                 this_genotype, ifelse(is.na(this_atp), "-", this_atp),
                 this_fraction, others_mean, this_ratio))
    section1_ratios <- c(section1_ratios, this_ratio)
    section1_is_wt_plusatp <- c(section1_is_wt_plusatp,
        this_genotype == "WT" && !is.na(this_atp) && this_atp == "plusATP")
}
emit("")
# Verdict rule, sharpened. The earlier specific/multiplicative/mixed labels did
# not fit what the data actually shows, which is a GRADIENT: the suspect gel is
# depressed everywhere, but plusATP lanes are hit far harder than noATP lanes,
# and within plusATP the highest-signal lane (WT) is hit hardest. That gradient
# is the fingerprint of a gel-level effect that scales with band signal (e.g. a
# narrow integration window clipping diffuse high-intensity shifted bands), NOT
# a lane-specific problem confined to WT plusATP. So instead of forcing one of
# three labels, report the plusATP-vs-noATP depression ratio directly and let
# the magnitude of that gap carry the verdict.
#
# Design note: the branch record predicted SPECIFIC (only WT plusATP low, rest
# of gel normal). The data refuted that -- every plusATP analyte is ~0.5-0.6 of
# its cross-gel norm. Keeping the raw per-analyte ratios above as the primary
# evidence; this block only summarizes their structure.
# Rebuild ratio-with-atp as a clean parallel pair so the summary cannot mispair
# a ratio with the wrong ATP state. This repeats the Section 1 loop's ratio math
# deliberately: the first loop's job was to print every analyte; this one's job
# is to tag each ratio with its ATP state for the gradient summary. Inlining the
# recompute is cheaper to read than threading a second output out of the printer
# loop above.
ratio_atp_labels <- character(0)
ratio_values_clean <- numeric(0)
for (row_index in seq_len(nrow(suspect_gel_rows))) {
    this_genotype <- suspect_gel_rows$genotype[row_index]
    this_atp <- suspect_gel_rows$atp[row_index]
    this_fraction <- suspect_gel_rows$bound_fraction[row_index]
    other_mask <- combined_rows$gel_key != SUSPECT_GEL_KEY &
                  combined_rows$genotype == this_genotype &
                  !is.na(combined_rows$genotype)
    if (is.na(this_atp)) other_mask <- other_mask & is.na(combined_rows$atp) else
        other_mask <- other_mask & !is.na(combined_rows$atp) &
                      combined_rows$atp == this_atp
    other_mask[is.na(other_mask)] <- FALSE
    if (sum(other_mask) == 0L || is.na(this_atp)) next
    ratio_atp_labels <- c(ratio_atp_labels, this_atp)
    ratio_values_clean <- c(ratio_values_clean,
                            this_fraction / mean(combined_rows$bound_fraction[other_mask]))
}
plusatp_ratios <- ratio_values_clean[ratio_atp_labels == "plusATP"]
noatp_ratios <- ratio_values_clean[ratio_atp_labels == "noATP"]
wt_plusatp_ratio <- section1_ratios[section1_is_wt_plusatp]
if (length(plusatp_ratios) >= 1L && length(noatp_ratios) >= 1L) {
    mean_plusatp_ratio <- mean(plusatp_ratios)
    mean_noatp_ratio <- mean(noatp_ratios)
    emit(sprintf("Suspect-gel depression: plusATP lanes at %.2f of cross-gel norm (mean),",
                 mean_plusatp_ratio))
    emit(sprintf("noATP lanes at %.2f. plusATP is hit %.1fx harder than noATP.",
                 mean_noatp_ratio, (1 - mean_plusatp_ratio) / max(1 - mean_noatp_ratio, 1e-6)))
    emit(sprintf("Within plusATP, WT (the highest-signal lane) is lowest at %.2f.",
                 wt_plusatp_ratio))
    if (mean_plusatp_ratio < 0.75 && mean_noatp_ratio > mean_plusatp_ratio &&
        length(wt_plusatp_ratio) == 1L && wt_plusatp_ratio <= min(plusatp_ratios) + 1e-9) {
        emit("VERDICT: GEL-LEVEL, SIGNAL-SCALED. The whole gel is down, plusATP more")
        emit("than noATP, WT plusATP most -- a gradient tracking band signal, not a")
        emit("WT-lane-specific fault. Consistent with a gel-level quantitation effect")
        emit("(see Section 5 integration windows). Per-gel WT normalization should")
        emit("substantially ABSORB it, which supports the Task B normalized plot.")
    } else {
        emit("VERDICT: pattern does not match the expected signal-scaled gradient;")
        emit("read the per-analyte ratios above and Section 5 before deciding.")
    }
}
emit("")

# ==============================================================================
# SECTION 2: Two-sided Grubbs on WT plusATP (n=6). WT noATP as negative control.
#
# Grubbs single-outlier statistic G = max|x - mean| / sd. Two-sided p from the
# exact t-based critical-value formula:
#   p = n * (1 - P(T <= t_crit)) * 2-sided, where
#   t_crit relates to G by G = ((n-1)/sqrt(n)) * sqrt(t^2 / (n - 2 + t^2)).
# Inlined here (no outliers package). Identical to outliers::grubbs.test,
# type = 10, two.sided = TRUE.
#
# WT noATP is the negative control: it is bimodal by SCREEN (0.11-0.13 vs
# 0.27-0.29), a batch effect, not a single outlier. A correct single-outlier
# test should NOT flag one point there. If it does, the screen effect is severe
# enough to distort the test and that is itself worth knowing.
# ==============================================================================
emit("--------------------------------------------------------------------------------")
emit("SECTION 2: Two-sided Grubbs, WT plusATP (test) and WT noATP (negative control)")
emit("--------------------------------------------------------------------------------")
for (cell_name in c("WT plusATP", "WT noATP")) {
    values <- if (cell_name == "WT plusATP") wt_plusatp$bound_fraction else wt_noatp$bound_fraction
    sample_size <- length(values)
    cell_mean <- mean(values)
    cell_sd <- sd(values)
    extreme_index <- which.max(abs(values - cell_mean))
    grubbs_g <- abs(values[extreme_index] - cell_mean) / cell_sd
    # Invert G to the corresponding t, then two-sided p with Bonferroni over n.
    t_squared <- (sample_size * (sample_size - 2) * grubbs_g^2) /
                 ((sample_size - 1)^2 - sample_size * grubbs_g^2)
    t_statistic <- sqrt(max(t_squared, 0))
    p_one_tail <- pt(t_statistic, df = sample_size - 2, lower.tail = FALSE)
    grubbs_p_two_sided <- min(sample_size * 2 * p_one_tail, 1)
    emit(sprintf("  %-11s n=%d  extreme=%.4f (%s rep %s)  G=%.4f  p(two-sided)=%.4f",
                 cell_name, sample_size, values[extreme_index],
                 (if (cell_name == "WT plusATP") wt_plusatp$screen else wt_noatp$screen)[extreme_index],
                 (if (cell_name == "WT plusATP") wt_plusatp$replicate else wt_noatp$replicate)[extreme_index],
                 grubbs_g, grubbs_p_two_sided))
    if (cell_name == "WT plusATP") {
        emit(if (grubbs_p_two_sided < 0.05)
                 "    -> flags the suspect at alpha=0.05 (expected)." else
                 "    -> does NOT flag at alpha=0.05.")
    } else {
        emit(if (grubbs_p_two_sided < 0.05)
                 "    -> WARNING: flags a single point in a cell whose real structure is a 3-vs-3 screen split; treat with care." else
                 "    -> correctly does NOT single out one point in a screen-split cell (good negative control).")
    }
}
emit("")

# ==============================================================================
# SECTION 3: Leave-one-out deltas on WT plusATP.
# Assumption-free description of how much the one point moves the cell summary.
# ==============================================================================
emit("--------------------------------------------------------------------------------")
emit("SECTION 3: Leave-one-out impact of the suspect on WT plusATP summary")
emit("--------------------------------------------------------------------------------")
wt_plusatp_all <- wt_plusatp$bound_fraction
suspect_in_cell <- wt_plusatp$gel_key == SUSPECT_GEL_KEY
wt_plusatp_without <- wt_plusatp_all[!suspect_in_cell]
mean_with <- mean(wt_plusatp_all); mean_without <- mean(wt_plusatp_without)
sd_with <- sd(wt_plusatp_all); sd_without <- sd(wt_plusatp_without)
se_with <- sd_with / sqrt(length(wt_plusatp_all))
se_without <- sd_without / sqrt(length(wt_plusatp_without))
cv_with <- sd_with / mean_with; cv_without <- sd_without / mean_without
suspect_z <- (suspect_fraction_found - mean_with) / sd_with
emit(sprintf("  with suspect    (n=%d): mean=%.4f sd=%.4f se=%.4f cv=%.3f",
             length(wt_plusatp_all), mean_with, sd_with, se_with, cv_with))
emit(sprintf("  without suspect (n=%d): mean=%.4f sd=%.4f se=%.4f cv=%.3f",
             length(wt_plusatp_without), mean_without, sd_without, se_without, cv_without))
emit(sprintf("  suspect z within cell (with suspect): %.3f sd", suspect_z))
emit(sprintf("  removing it moves the cell mean by %+.4f (%.1f%%) and cuts CV from %.3f to %.3f",
             mean_without - mean_with,
             100 * (mean_without - mean_with) / mean_with, cv_with, cv_without))
emit("")

# ==============================================================================
# SECTION 4: Screen-split report on both WT cells.
# Describes the between-screen offset. Names the PMT/batch boundary as a
# correlate, NOT a proven numeric cause (per branch finding 3a).
# ==============================================================================
emit("--------------------------------------------------------------------------------")
emit("SECTION 4: WT by screen (between-batch offset)")
emit("--------------------------------------------------------------------------------")
for (cell_name in c("WT plusATP", "WT noATP")) {
    cell_data <- if (cell_name == "WT plusATP") wt_plusatp else wt_noatp
    for (screen_name in c("wt-1-3-4", "4r-5-6")) {
        screen_values <- cell_data$bound_fraction[cell_data$screen == screen_name]
        emit(sprintf("  %-11s %-9s mean=%.4f  values: %s",
                     cell_name, screen_name, mean(screen_values),
                     paste(sprintf("%.4f", sort(screen_values)), collapse = ", ")))
    }
}
emit("")
emit("Note: wt-1-3-4 was scanned at PMT 655 V, 4r-5-6 at 300 V (metadata branch).")
emit("The WT noATP offset (~0.12 vs ~0.28) tracks this screen/batch boundary. A pure")
emit("gain change cancels in a ratio, so PMT is a MARKER of two separately acquired")
emit("batches, not proven to be the numeric cause. The folded WT noATP bar therefore")
emit("averages two populations differing for at least an acquisition-linked reason.")
emit("")

# ==============================================================================
# SECTION 5: Integration-window check on the suspect gel.
# The metadata rules out scanner causes but cannot see quantitation choices. The
# bound/free region-mm spans live in the CSV and CAN differ per gel. If the
# suspect gel integrated a narrower bound window or wider free window, that could
# depress its ratio independent of biology. Cheap, in-scope, closes that gap.
# ==============================================================================
emit("--------------------------------------------------------------------------------")
emit("SECTION 5: Integration windows (bound/free region spans) on the suspect gel")
emit("--------------------------------------------------------------------------------")
combined_rows$bound_span_mm <- combined_rows$bound_region_end_millimetres -
                               combined_rows$bound_region_start_millimetres
combined_rows$free_span_mm <- combined_rows$free_region_end_millimetres -
                              combined_rows$free_region_start_millimetres
per_gel_bound_span <- tapply(combined_rows$bound_span_mm, combined_rows$gel_key,
                             function(v) mean(v, na.rm = TRUE))
per_gel_free_span <- tapply(combined_rows$free_span_mm, combined_rows$gel_key,
                            function(v) mean(v, na.rm = TRUE))
for (gel_key_name in names(per_gel_bound_span)) {
    marker <- if (gel_key_name == SUSPECT_GEL_KEY) "  <- suspect" else ""
    emit(sprintf("  %-16s bound_span=%.1f mm  free_span=%.1f mm%s",
                 gel_key_name, per_gel_bound_span[[gel_key_name]],
                 per_gel_free_span[[gel_key_name]], marker))
}
suspect_bound_span <- per_gel_bound_span[[SUSPECT_GEL_KEY]]
other_bound_spans <- per_gel_bound_span[names(per_gel_bound_span) != SUSPECT_GEL_KEY]
emit(sprintf("Suspect bound_span %.1f mm vs other-gel mean %.1f mm (range %.1f-%.1f).",
             suspect_bound_span, mean(other_bound_spans),
             min(other_bound_spans), max(other_bound_spans)))
emit("Interpret: a markedly NARROWER bound window or WIDER free window on the suspect")
emit("could depress its ratio for a quantitation reason. If spans are comparable, this")
emit("cause is ruled out too and the low value is neither scan nor window driven.")
emit("")

# ==============================================================================
# SECTION 6: Assumption checks for the parametric tests above.
# Shapiro on model residuals with and WITHOUT the suspect (the outlier itself
# breaks normality, so both are reported). Levene across cells inlined (no car).
# Studentized residuals from a gel-factor model, raw and logit, as robustness --
# gel is a factor because the suspect gel contaminates multiple cells and a
# naive pooled error would be inflated (per earlier adversarial note).
# ==============================================================================
emit("--------------------------------------------------------------------------------")
emit("SECTION 6: Assumption checks and gel-factor studentized residuals")
emit("--------------------------------------------------------------------------------")
# Model frame: measured analytes only (drop No ORC reference), complete rows.
model_data <- combined_rows[!is.na(combined_rows$atp) &
                            combined_rows$genotype != "No ORC" &
                            !is.na(combined_rows$bound_fraction), ]
model_data$genotype <- factor(model_data$genotype)
model_data$atp <- factor(model_data$atp)
model_data$gel_key_factor <- factor(model_data$gel_key)
model_data$is_suspect <- model_data$gel_key == SUSPECT_GEL_KEY &
                         model_data$genotype == SUSPECT_GENOTYPE &
                         model_data$atp == SUSPECT_ATP

# Raw-scale model with gel as a factor. Only include gel_key_factor if it has >1
# level with data; it does. drop unused factor levels to keep lm well posed.
model_data$genotype <- droplevels(model_data$genotype)
fit_raw <- lm(bound_fraction ~ genotype + atp + gel_key_factor, data = model_data)
studentized <- rstudent(fit_raw)
suspect_studentized <- studentized[model_data$is_suspect]
n_model <- nrow(model_data)
bonferroni_p <- pmin(2 * pt(abs(suspect_studentized),
                            df = df.residual(fit_raw) - 1, lower.tail = FALSE) * n_model, 1)
emit(sprintf("Gel-factor model (bound_fraction ~ genotype + atp + gel), n=%d.", n_model))
emit(sprintf("  Suspect externally studentized residual: %.3f  (Bonferroni p=%.4f over n=%d)",
             suspect_studentized, bonferroni_p, n_model))

# Logit robustness: same model on log(f/(1-f)). Guards proportions off 0/1.
eps <- 1e-6
model_data$logit_fraction <- log((model_data$bound_fraction + eps) /
                                  (1 - model_data$bound_fraction + eps))
fit_logit <- lm(logit_fraction ~ genotype + atp + gel_key_factor, data = model_data)
suspect_studentized_logit <- rstudent(fit_logit)[model_data$is_suspect]
emit(sprintf("  Same model on logit scale: suspect studentized residual %.3f (agreement check).",
             suspect_studentized_logit))

# Shapiro on residuals, with and without the suspect.
shapiro_with <- shapiro.test(residuals(fit_raw))
resid_without <- residuals(fit_raw)[!model_data$is_suspect]
shapiro_without <- shapiro.test(resid_without)
emit(sprintf("  Shapiro on residuals  with suspect: W=%.3f p=%.4f",
             shapiro_with$statistic, shapiro_with$p.value))
emit(sprintf("  Shapiro on residuals  without suspect: W=%.3f p=%.4f",
             shapiro_without$statistic, shapiro_without$p.value))
emit("  (Normality improving markedly once the suspect is dropped is evidence FOR the")
emit("   point being the anomaly, not evidence against the method.)")

# Levene across (genotype x atp) cells, median-centered (Brown-Forsythe),
# inlined. Statistic W ~ F(k-1, N-k). Identical to car::leveneTest default.
cell_id <- interaction(model_data$genotype, model_data$atp, drop = TRUE)
group_median <- ave(model_data$bound_fraction, cell_id, FUN = median)
abs_dev <- abs(model_data$bound_fraction - group_median)
levene_fit <- lm(abs_dev ~ cell_id)
levene_anova <- anova(levene_fit)
levene_f <- levene_anova[["F value"]][1]
levene_p <- levene_anova[["Pr(>F)"]][1]
emit(sprintf("  Levene (Brown-Forsythe) across %d cells: F=%.3f p=%.4f",
             nlevels(cell_id), levene_f, levene_p))
emit("  (If variance heterogeneity is driven by the suspect's cell, that is itself")
emit("   part of the outlier signal, not an independent violation.)")
emit("")

# Q-Q PDF of residuals, with and without suspect.
pdf("task_a_qq_residuals.pdf", width = 9, height = 4.5)
par(mfrow = c(1, 2))
qqnorm(residuals(fit_raw), main = "Residuals with suspect"); qqline(residuals(fit_raw))
qqnorm(resid_without, main = "Residuals without suspect"); qqline(resid_without)
dev.off()
emit("Wrote task_a_qq_residuals.pdf (Q-Q with and without the suspect).")
emit("")

# ==============================================================================
# SECTION 7: Consolidated recommendation and the exclusion key for plotting.
# The script does not decide; it states the convergent read and the exact key
# the plotting script's manual exclusion variable must match.
# ==============================================================================
emit("--------------------------------------------------------------------------------")
emit("SECTION 7: Recommendation and exclusion key")
emit("--------------------------------------------------------------------------------")
emit("Convergent read across sections 1-6 (fill in from the numbers above):")
emit("  - Section 1 verdict specific vs multiplicative decides normalize-away vs flag.")
emit("  - Section 2 Grubbs p and section 3 leave-one-out quantify how lone the point is.")
emit("  - Section 5 rules quantitation windows in or out.")
emit("Given the metadata (no scan artifact) and an unexplained cause, the defensible")
emit("action is to FLAG the point and show the figure both ways; silent removal is not")
emit("warranted on statistics alone. Final keep/drop is the PI's call.")
emit("")
emit("If excluding, the plotting script's manual variable must match this stable key:")
emit(paste0("    EXCLUDE_GEL_KEY <- \"", SUSPECT_GEL_KEY, "\""))
emit("Match on screen + replicate (gel id), NOT a lane index. The plotting script")
emit("must stop() if the key matches zero rows, so a stale key fails loudly rather")
emit("than silently plotting the point back in.")
emit("")
emit("================================================================================")
emit("END OF REPORT")
emit("================================================================================")

writeLines(REPORT_LINES, "task_a_outlier_diagnostic.txt")
