#!/usr/bin/env bash
#
# Per-session driver for the gel densitometry pipeline. Manually run, by you, from
# images/. This is a TEMPLATE: copy it to an active, gitignored working copy for the
# session (session.sh) and edit the pasted paths there, so the template stays clean:
#
#     cp session_template.sh session.sh    # session.sh is gitignored
#     nvim session.sh                      # paste your gel paths into GEL_PASTED_PATHS
#     bash session.sh                       # or source it to keep the vars in your shell
#
# It does NOT replace the scripts' own validation. It reduces the two real frictions
# of driving the CLI from WSL while the imaging/clicking happens in Windows:
#   1. Windows paths (from File Explorer "Copy as path") -> WSL paths, once, in bulk.
#   2. Looping the per-gel stages over a whole session's worth of gels, proceeding on
#      the gels whose inputs are present and reporting (not silently skipping) the rest.
# Schema validation stays in the scripts; this only checks that the inputs EXIST.
#
# Conventions match the repo's Python style where they apply to bash: full
# descriptive names, no abbreviations, comments state why, ASCII only. Every path is
# used double-quoted so spaces, commas, and square brackets in real gel filenames are
# inert; do not remove the quotes.

set -u  # unset variable is a bug; catch it. NOT set -e: a missing input on ONE gel
        # must not abort the session's ready gels (see the readiness loop below).

# -----------------------------------------------------------------------------
# 0. Paste your gel paths here, one per line, between the quotes.
# -----------------------------------------------------------------------------
# In Windows File Explorer, right-click the gel image (or its _gel_analysis folder)
# and "Copy as path". Paste here. "Copy as path" wraps the value in double quotes and
# uses backslashes; both are handled below, so paste verbatim. Either the .tif or the
# _gel_analysis folder works; both resolve to the same analysis directory.
#
# Example (delete these, paste yours):
#   "C:\Users\luis\gels\20260818 rotated, LM-0008 [s0012].tif"
#   "C:\Users\luis\gels\20260819 rotated, LM-0008 [s0013].tif"
GEL_PASTED_PATHS=(
    ""
)

# The picker / plot region, shared across every gel in the session. Kept here so one
# edit changes the whole session and every gel is extracted over the same window.
# (The aggregator currently requires the same selector across gels; see HOWTO.md for
# the offset caveat.) Millimetres.
REGION_START_MILLIMETRES=31.3
REGION_END_MILLIMETRES=46.1
REGION_NET_BASELINE="none"

GEL_ANALYSIS_SUFFIX="_gel_analysis"
LANE_PROFILE_FILENAME="manual_lane_profiles.csv"
LANE_ROI_ZIP_FILENAME="lane_rois.zip"
SAMPLE_SHEET_FILENAME="sample_sheet.csv"

# -----------------------------------------------------------------------------
# 1. Resolve each pasted Windows path to a WSL analysis directory.
# -----------------------------------------------------------------------------
# wslpath does the Windows->WSL conversion (C:\... -> /mnt/c/...). It is fed the raw
# pasted value with its outer double-quotes stripped; backslashes inside are fine for
# wslpath. Whether the paste was the image or the folder, we resolve to the folder:
# the export macro names it <image_stem>_gel_analysis next to the image.
RESOLVED_ANALYSIS_DIRECTORIES=()
for GEL_PASTED_PATH in "${GEL_PASTED_PATHS[@]}"; do
    if [ -z "$GEL_PASTED_PATH" ]; then
        continue  # skip the empty template placeholder line
    fi

    # Strip the surrounding double-quotes that "Copy as path" adds. Only the outer
    # pair; a path cannot legally contain a double-quote, so this is safe.
    GEL_PASTED_PATH_UNQUOTED="${GEL_PASTED_PATH%\"}"
    GEL_PASTED_PATH_UNQUOTED="${GEL_PASTED_PATH_UNQUOTED#\"}"

    # Convert to a WSL path. If it already looks like a POSIX path (starts with /),
    # wslpath is a no-op but still safe, so we do not special-case it.
    GEL_WSL_PATH="$(wslpath -u "$GEL_PASTED_PATH_UNQUOTED")"

    # Resolve to the analysis directory. Two paste cases, one result:
    #   - a directory already ending in _gel_analysis: use as is.
    #   - anything else (the image file): the sibling <stem>_gel_analysis.
    if [ -d "$GEL_WSL_PATH" ] && [ "${GEL_WSL_PATH%$GEL_ANALYSIS_SUFFIX}" != "$GEL_WSL_PATH" ]; then
        RESOLVED_ANALYSIS_DIRECTORY="$GEL_WSL_PATH"
    else
        PARENT_DIRECTORY="$(dirname "$GEL_WSL_PATH")"
        IMAGE_STEM="$(basename "${GEL_WSL_PATH%.*}")"
        RESOLVED_ANALYSIS_DIRECTORY="$PARENT_DIRECTORY/${IMAGE_STEM}${GEL_ANALYSIS_SUFFIX}"
    fi
    RESOLVED_ANALYSIS_DIRECTORIES+=("$RESOLVED_ANALYSIS_DIRECTORY")
done

if [ "${#RESOLVED_ANALYSIS_DIRECTORIES[@]}" -eq 0 ]; then
    echo "[session] no paths in GEL_PASTED_PATHS; paste at least one and re-run." >&2
    return 2 2>/dev/null || exit 2
fi

# -----------------------------------------------------------------------------
# 2. Readiness report. For each gel, are the three manual inputs present?
# -----------------------------------------------------------------------------
# The inputs a human must have produced before the CLI can run: the ImageJ export
# (profile CSV + ROI zip) and the authored sample sheet. Existence only; the scripts
# check schema. Gels missing any input are reported and then skipped by the loops in
# step 3, so a half-prepared gel never aborts the ready ones.
READY_ANALYSIS_DIRECTORIES=()
printf '%s\n' "[session] readiness (profile / roizip / sheet):"
for RESOLVED_ANALYSIS_DIRECTORY in "${RESOLVED_ANALYSIS_DIRECTORIES[@]}"; do
    HAS_PROFILE="no"; HAS_ROI_ZIP="no"; HAS_SAMPLE_SHEET="no"
    [ -f "$RESOLVED_ANALYSIS_DIRECTORY/$LANE_PROFILE_FILENAME" ] && HAS_PROFILE="yes"
    [ -f "$RESOLVED_ANALYSIS_DIRECTORY/$LANE_ROI_ZIP_FILENAME" ] && HAS_ROI_ZIP="yes"
    [ -f "$RESOLVED_ANALYSIS_DIRECTORY/$SAMPLE_SHEET_FILENAME" ] && HAS_SAMPLE_SHEET="yes"

    printf '  [%s/%s/%s] %s\n' \
        "$HAS_PROFILE" "$HAS_ROI_ZIP" "$HAS_SAMPLE_SHEET" \
        "$(basename "$RESOLVED_ANALYSIS_DIRECTORY")"

    if [ "$HAS_PROFILE" = "yes" ] && [ "$HAS_ROI_ZIP" = "yes" ] && [ "$HAS_SAMPLE_SHEET" = "yes" ]; then
        READY_ANALYSIS_DIRECTORIES+=("$RESOLVED_ANALYSIS_DIRECTORY")
    fi
done

if [ "${#READY_ANALYSIS_DIRECTORIES[@]}" -eq 0 ]; then
    echo "[session] no gel has all three inputs; nothing to run. Finish the export/sheet, re-run." >&2
    return 1 2>/dev/null || exit 1
fi
printf '%s\n' "[session] ${#READY_ANALYSIS_DIRECTORIES[@]} of ${#RESOLVED_ANALYSIS_DIRECTORIES[@]} gel(s) ready."

# -----------------------------------------------------------------------------
# 3. Per-gel stages, looped over the ready gels. Stop the session if a stage FAILS
#    on a gel (a real error, not a missing input): a bad analyze must not be plotted.
# -----------------------------------------------------------------------------
# These are the live pipeline stages in order, per gel. The cross-gel stages
# (validate_manifest, aggregate_repeats) run once, after, and are left commented
# because they need a manifest you point at your session's gels (see HOWTO.md).
for READY_ANALYSIS_DIRECTORY in "${READY_ANALYSIS_DIRECTORIES[@]}"; do
    echo "[session] === $(basename "$READY_ANALYSIS_DIRECTORY") ==="

    uv run validate_sample_sheet.py "$READY_ANALYSIS_DIRECTORY" || {
        echo "[session] validate_sample_sheet failed; fix the sheet, re-run this gel." >&2
        continue  # this gel is not usable; move to the next rather than plotting garbage
    }

    uv run analyze_gel.py "$READY_ANALYSIS_DIRECTORY" || {
        echo "[session] analyze_gel failed on this gel; skipping its downstream stages." >&2
        continue
    }

    uv run extract_lane_values.py "$READY_ANALYSIS_DIRECTORY" \
        --region "$REGION_START_MILLIMETRES" "$REGION_END_MILLIMETRES" \
        --net-baseline "$REGION_NET_BASELINE" || {
        echo "[session] extract failed on this gel; skipping plot." >&2
        continue
    }

    EXTRACT_CSV="$READY_ANALYSIS_DIRECTORY/extract_region_${REGION_START_MILLIMETRES}-${REGION_END_MILLIMETRES}mm_${REGION_NET_BASELINE}.csv"
    uv run plot_single_experiment.py "$READY_ANALYSIS_DIRECTORY" \
        --extract-csv "$EXTRACT_CSV" || {
        echo "[session] plot failed on this gel." >&2
        continue
    }
done

# -----------------------------------------------------------------------------
# 4. Cross-gel aggregation (run once, after every gel above is measured).
# -----------------------------------------------------------------------------
# Uncomment and point MANIFEST at a manifest.csv whose analysis_path rows are this
# session's gels. See HOWTO.md for authoring the manifest and the region caveat.
#
# MANIFEST="/path/to/manifest.csv"
# uv run validate_manifest.py "$MANIFEST"
# uv run aggregate_repeats.py "$MANIFEST" \
#     --selector "region_${REGION_START_MILLIMETRES}-${REGION_END_MILLIMETRES}mm_${REGION_NET_BASELINE}"

echo "[session] done."
