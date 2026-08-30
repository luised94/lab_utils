#!/usr/bin/env python3
# ==============================================================================
# harvest_inf_provenance.py
#
# Read the Amersham Typhoon .inf sidecar for each gel and emit one provenance
# table plus a per-file found/missing report. The single question this answers:
# were the six working scans acquired at the same sensitivity, and is any of them
# sitting near the 16-bit saturation ceiling? The .inf carries PMT voltage
# (S,V=...) and the actual max pixel value present in the image (S,RangeHigh=...);
# together those decide whether the rep-2 outlier is a scan artifact or real.
#
# This is a self-contained companion to the R analysis/plotting scripts. It reads
# ONLY .inf text -- no TIFF pixels, no CSVs -- because the .inf already carries the
# saturation-relevant fields, so a pixel-level census would be verifying a question
# the metadata answers directly. If a later step needs per-pixel confirmation, that
# is a separate script.
#
# Inputs are the Windows paths to the .inf files, pasted verbatim below. The script
# normalizes each to a WSL path with `wslpath -u` at runtime (falling back to the
# path as-is if wslpath is unavailable), so the same INF_FILES_WINDOWS list works
# whether or not it is run inside WSL, and the paths are not hard-coded to /mnt/c.
#
# Why match files by (experiment id, task id) integer pair rather than filename
# stem: the .tif, .inf, .img, .gel for one gel do NOT share a stem
# (..._rotated_gelshift_... vs ..._[Phosphor].inf), casing varies (LM vs Lm), and
# the task token is written inconsistently (s0001 in the filename but s1 inside the
# .img header line). Deriving one path from another by string-munging the stem
# silently mismatches. The robust key is LM-<digits> and s<digits> parsed to
# integers and compared case-insensitively, so s1 == s0001. Each .inf is labelled
# with the analytical gel (the bound/free CSV) that shares its id pair, and any
# .inf that cannot be matched, or any expected gel with no .inf, is reported rather
# than dropped.
#
# Usage (inside WSL):
#     python3 harvest_inf_provenance.py
#
# Output:
#     inf_provenance.csv   one row per expected gel: ids, PMT V, RangeHigh, headroom
#                          to the 16-bit ceiling, pixel size, dimensions, bit depth,
#                          scan date, laser/filter, and the matched .inf path (or a
#                          NO_INF marker). Written next to this script's OUTPUT_DIR.
#     also printed to stdout as an aligned table for a quick read.
# ==============================================================================

import csv
import os
import re
import subprocess
import sys

# ------------------------------------------------------------------------------
# The .inf files, as Windows paths, pasted verbatim. Replace these with your local
# paths. The 20210908 gel (wt-1-3-4 replicate 1) has no .inf and no .gel; it is
# listed as EXPECTED_GELS below with inf_windows_path = None so the report shows an
# explicit NO_INF row for it rather than a silently missing gel.
#
# The two deliberate 4000 V saturated re-scans are intentionally NOT listed: they
# are a saturation positive control, not part of the quantified data.
# ------------------------------------------------------------------------------

# Each expected gel ties one analytical bound/free CSV to the .inf that shares its
# (experiment id, task id) pair. csv_basename is only a human label here so the
# provenance row names the gel the way the analysis does; nothing is read from it.
EXPECTED_GELS = [
    {
        "screen": "wt-1-3-4",
        "replicate": 1,
        "csv_basename": "gel_shift_ratio_16.8-23.6mm_bound_over_46.6-56mm_free.csv",
        # No .inf and no .gel were kept for this scan. Reported as NO_INF.
        "inf_windows_path": None,
    },
    {
        "screen": "wt-1-3-4",
        "replicate": 2,  # the suspect gel (WT plusATP bound_fraction 0.298)
        "csv_basename": "gel_shift_ratio_23.6-27.4mm_bound_over_38.9-45.6mm_free.csv",
        "inf_windows_path": r"C:\Users\Luised94\Desktop\lab\experiments\20260716_LM-0009_gs_1-3-4-repeats\20260723_LM-0009_s0001_wt,1,3,4-sofa-repeat-2-1000_[Phosphor].inf",
    },
    {
        "screen": "wt-1-3-4",
        "replicate": 3,
        "csv_basename": "gel_shift_ratio_23.8-29.4mm_bound_over_39.8-45.6mm_free.csv",
        "inf_windows_path": r"C:\Users\Luised94\Desktop\lab\experiments\20260716_LM-0009_gs_1-3-4-repeats\20260724_LM-0009_s0002_wt,1,3,4-sofa-repeat-3-1000_[Phosphor].inf",
    },
    {
        "screen": "4r-5-6",
        "replicate": 1,
        "csv_basename": "gel_shift_ratio_16.8-21.4mm_bound_over_42.8-52.4mm_free.csv",
        "inf_windows_path": r"C:\Users\Luised94\Desktop\lab\experiments\20260803_LM-0013_gs_ORC5-6-sofa-data-for-analysis\20220303_LM-0013_s0001_gelshift_wt,4r,5,6-sofa-repeat-1-1000_Phosphor.inf",
    },
    {
        "screen": "4r-5-6",
        "replicate": 2,
        "csv_basename": "gel_shift_ratio_15.6-22.2mm_bound_over_43.2-54.6mm_free.csv",
        "inf_windows_path": r"C:\Users\Luised94\Desktop\lab\experiments\20260803_LM-0013_gs_ORC5-6-sofa-data-for-analysis\20220304_LM-0013_s0002_gelshift_wt,4r,5,6-sofa-repeat-2-1000_Phosphor.inf",
    },
    {
        "screen": "4r-5-6",
        "replicate": 3,
        "csv_basename": "gel_shift_ratio_19-23.8mm_bound_over_44.2-53.4mm_free.csv",
        "inf_windows_path": r"C:\Users\Luised94\Desktop\lab\experiments\20260803_LM-0013_gs_ORC5-6-sofa-data-for-analysis\20220307_LM-0013_s0003_gelshift_wt,4r,5,6-sofa-repeat-3-1000_Phosphor.inf",
    },
]

# Where inf_provenance.csv is written. "." keeps it next to wherever the script is
# run; change if you want it elsewhere.
OUTPUT_DIRECTORY = "."

# The Typhoon exports 16-bit linear, so a pixel maxes at 65535. RangeHigh is the
# largest pixel value actually present in the image; the closer to this ceiling, the
# nearer the brightest band is to clipping. Kept as a named constant so the headroom
# arithmetic below reads plainly.
SIXTEEN_BIT_CEILING = 65535


def normalize_windows_path_to_wsl(windows_path):
    # Convert a Windows path to the WSL mount path with `wslpath -u`. Falls back to
    # the input unchanged if wslpath is not on PATH (i.e. not running inside WSL),
    # so the caller can still attempt to open it -- the existence check downstream
    # is what actually decides whether the path is usable, not this conversion.
    try:
        completed = subprocess.run(
            ["wslpath", "-u", windows_path], capture_output=True, text=True, check=True
        )
        return completed.stdout.strip()
    except FileNotFoundError:
        # wslpath itself is absent: not in WSL. Return as-is and let the open fail
        # loudly later if the path is not reachable in this environment.
        return windows_path
    except subprocess.CalledProcessError as conversion_error:
        # wslpath ran but rejected the path. Surface its stderr rather than guessing.
        raise RuntimeError(
            "wslpath failed to convert '{}': {}".format(
                windows_path, conversion_error.stderr.strip()
            )
        ) from conversion_error


def extract_experiment_and_task_ids(text_to_search):
    # Pull (experiment_id_integer, task_id_integer) from a filename or path. The
    # experiment id is LM-<digits> (case-insensitive: matches LM-0009 and Lm-0009),
    # the task id is s<digits> with any leading zeros (matches s0001 and s1, which
    # are the same scan written two ways). Both returned as ints so s1 == s0001 and
    # LM-0009 == Lm-0009 compare equal. Returns (None, None) when a token is absent.
    experiment_match = re.search(r"LM-(\d+)", text_to_search, re.IGNORECASE)
    task_match = re.search(r"[_\-/\\]s0*(\d+)", text_to_search, re.IGNORECASE)
    experiment_id = int(experiment_match.group(1)) if experiment_match else None
    task_id = int(task_match.group(1)) if task_match else None
    return experiment_id, task_id


def parse_inf_keyed_block(inf_text):
    # The .inf is a positional header (meaning-by-line-number, and the layout is not
    # identical between scans) followed by a `*** more info ***` block of Key=Value
    # lines. Only the keyed block is trusted for the fields that matter; positional
    # lines are read separately and conservatively below. The keyed lines look like
    # "S,V=655" or "H,ScaleType=Linear": a one-letter section tag, a comma, the key,
    # an equals, then the value. Returns a dict keyed by the bare key (e.g. "V",
    # "RangeHigh", "ScaleType"). If a key appears twice the last wins; none do here.
    keyed_values = {}
    for raw_line in inf_text.splitlines():
        line = raw_line.strip()
        keyed_match = re.match(r"^[A-Za-z],([^=]+)=(.*)$", line)
        if keyed_match:
            keyed_values[keyed_match.group(1).strip()] = keyed_match.group(2).strip()
    return keyed_values


def parse_inf_positional_header(inf_text):
    # Read the small set of positional-header facts that are NOT in the keyed block:
    # bit depth and the two pixel dimensions. Confirmed against both pasted samples:
    #   line 0: BAS_IMAGE_FILE
    #   line 1: the internal .img filename
    #   lines 2,3: pixel size in micrometres, twice (50 and 50)
    #   line 4: bit depth (16)
    #   lines 5,6: pixel dimensions (6000 4140 for the suspect, 4000 4000 for rep 3)
    # The suspect's dimensions differ from rep 3's, which is exactly why field
    # meaning cannot be inferred from value alone and must be read by fixed position.
    # Anything ambiguous (the lone 1000, the 5) is deliberately NOT interpreted here.
    # Returns a dict; any field that is missing or non-numeric comes back as None so
    # the report shows a blank rather than a wrong guess.
    lines = [line.strip() for line in inf_text.splitlines()]

    def integer_at(line_index):
        if line_index < len(lines) and re.fullmatch(r"\d+", lines[line_index]):
            return int(lines[line_index])
        return None

    return {
        "pixel_size_micrometres": integer_at(2),
        "pixel_size_micrometres_check": integer_at(3),
        "bit_depth": integer_at(4),
        "pixel_dimension_a": integer_at(5),
        "pixel_dimension_b": integer_at(6),
    }


def read_one_inf(inf_windows_path):
    # Convert, open, and parse a single .inf. Returns a dict of the fields the
    # provenance table needs, or raises if the file cannot be read so the caller can
    # record the failure explicitly instead of producing a half-filled row silently.
    wsl_path = normalize_windows_path_to_wsl(inf_windows_path)
    if not os.path.exists(wsl_path):
        raise FileNotFoundError(
            "resolved path does not exist: {} (from Windows path {})".format(
                wsl_path, inf_windows_path
            )
        )
    # .inf files are short ASCII/latin-1 text. latin-1 never raises on a stray byte,
    # which a strict utf-8 read could, and every field of interest is plain ASCII.
    with open(wsl_path, "r", encoding="latin-1") as inf_handle:
        inf_text = inf_handle.read()

    keyed = parse_inf_keyed_block(inf_text)
    positional = parse_inf_positional_header(inf_text)

    # PMT voltage. S,V sits in the PMT group (PMTType is the line above it in both
    # samples), so V is the PMT voltage, distinct from LaserPowerMode which is
    # surfaced separately below so the two are never conflated.
    pmt_voltage_text = keyed.get("V")
    pmt_voltage = (
        int(pmt_voltage_text)
        if pmt_voltage_text and re.fullmatch(r"\d+", pmt_voltage_text)
        else None
    )

    # RangeHigh: the largest pixel value present. Headroom to the 16-bit ceiling is
    # the direct saturation read -- a value near SIXTEEN_BIT_CEILING means the
    # brightest band is close to clipping; a value far below it (as in both samples:
    # 21139 and 33664) means there is no saturation to find.
    range_high_text = keyed.get("RangeHigh")
    range_high = (
        int(range_high_text)
        if range_high_text and re.fullmatch(r"\d+", range_high_text)
        else None
    )
    headroom_to_ceiling = (
        SIXTEEN_BIT_CEILING - range_high if range_high is not None else None
    )

    # The date line is the first line matching a "Www Mmm DD HH:MM:SS YYYY" shape;
    # it is positional in practice but located by pattern so a shifted header does
    # not hand back the wrong line.
    scan_date = None
    for line in inf_text.splitlines():
        stripped = line.strip()
        if re.match(r"^[A-Za-z]{3} [A-Za-z]{3} +\d+ [\d:]+ \d{4}$", stripped):
            scan_date = stripped
            break

    return {
        "pmt_voltage": pmt_voltage,
        "range_high": range_high,
        "headroom_to_ceiling": headroom_to_ceiling,
        "scale_type": keyed.get("ScaleType"),
        "laser_name": keyed.get("LaserName"),
        "filter_name": keyed.get("FilterName"),
        "laser_power_mode": keyed.get("LaserPowerMode"),
        "range_low": keyed.get("RangeLow"),
        "pixel_size_micrometres": positional["pixel_size_micrometres"],
        "bit_depth": positional["bit_depth"],
        "pixel_dimension_a": positional["pixel_dimension_a"],
        "pixel_dimension_b": positional["pixel_dimension_b"],
        "scan_date": scan_date,
        "inf_experiment_id": extract_experiment_and_task_ids(inf_windows_path)[0],
        "inf_task_id": extract_experiment_and_task_ids(inf_windows_path)[1],
        "resolved_wsl_path": wsl_path,
    }


# ------------------------------------------------------------------------------
# Driver. Wrapped in main() behind the __main__ guard below so importing this file
# (to test the parse functions) does not try to read the Windows .inf paths. The
# body stays one flat block -- no helper split -- matching the rest of the script.
# ------------------------------------------------------------------------------
def main():
    # ------------------------------------------------------------------------------
    # Walk the expected gels, read each .inf, and assemble the provenance rows. A gel
    # with inf_windows_path = None (no .inf kept) becomes an explicit NO_INF row. A
    # gel whose .inf is listed but unreadable becomes a READ_ERROR row carrying the
    # reason, so a missing or moved file fails loudly in the table not vanishes.
    # ------------------------------------------------------------------------------
    provenance_rows = []
    for gel in EXPECTED_GELS:
        expected_experiment_id, expected_task_id = extract_experiment_and_task_ids(
            gel["csv_basename"] + " " + str(gel["inf_windows_path"])
        )

        base_row = {
            "screen": gel["screen"],
            "replicate": gel["replicate"],
            "csv_basename": gel["csv_basename"],
        }

        if gel["inf_windows_path"] is None:
            base_row.update({"status": "NO_INF"})
            provenance_rows.append(base_row)
            continue

        try:
            parsed = read_one_inf(gel["inf_windows_path"])
        except Exception as read_error:  # noqa: BLE001 -- want the reason in the table
            base_row.update({"status": "READ_ERROR", "detail": str(read_error)})
            provenance_rows.append(base_row)
            continue

        # The task id inside the wt-1-3-4 .inf names (s0001/s0002) is the scanner task
        # number, which is not the same integer as the analysis replicate (1/2/3). So do
        # not assert task_id == replicate; only assert both ids parsed. Record them so a
        # human can eyeball that each row's ids line up with its screen/replicate.
        inf_experiment_id = parsed["inf_experiment_id"]
        inf_task_id = parsed["inf_task_id"]
        status = (
            "OK"
            if (inf_experiment_id is not None and inf_task_id is not None)
            else "ID_PARSE_INCOMPLETE"
        )

        base_row.update(
            {
                "status": status,
                "inf_experiment_id_LM": inf_experiment_id,
                "inf_task_id_s": inf_task_id,
                "pmt_voltage_V": parsed["pmt_voltage"],
                "range_high": parsed["range_high"],
                "headroom_to_16bit_ceiling": parsed["headroom_to_ceiling"],
                "scale_type": parsed["scale_type"],
                "bit_depth": parsed["bit_depth"],
                "pixel_size_um": parsed["pixel_size_micrometres"],
                "pixel_dim_a": parsed["pixel_dimension_a"],
                "pixel_dim_b": parsed["pixel_dimension_b"],
                "laser": parsed["laser_name"],
                "filter": parsed["filter_name"],
                "laser_power_mode": parsed["laser_power_mode"],
                "range_low": parsed["range_low"],
                "scan_date": parsed["scan_date"],
                "resolved_wsl_path": parsed["resolved_wsl_path"],
            }
        )
        provenance_rows.append(base_row)

    # Stable order: screen then replicate, so the table is deterministic run to run.
    provenance_rows.sort(key=lambda row: (row["screen"], row["replicate"]))

    # ------------------------------------------------------------------------------
    # Write inf_provenance.csv. The column union is taken across all rows so NO_INF and
    # READ_ERROR rows (which carry fewer fields) still line up under a stable header.
    # ------------------------------------------------------------------------------
    column_order = [
        "screen",
        "replicate",
        "status",
        "csv_basename",
        "inf_experiment_id_LM",
        "inf_task_id_s",
        "pmt_voltage_V",
        "range_high",
        "headroom_to_16bit_ceiling",
        "range_low",
        "scale_type",
        "bit_depth",
        "pixel_size_um",
        "pixel_dim_a",
        "pixel_dim_b",
        "laser",
        "filter",
        "laser_power_mode",
        "scan_date",
        "detail",
        "resolved_wsl_path",
    ]
    os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)
    output_csv_path = os.path.join(OUTPUT_DIRECTORY, "inf_provenance.csv")
    with open(output_csv_path, "w", newline="", encoding="utf-8") as csv_handle:
        writer = csv.DictWriter(
            csv_handle, fieldnames=column_order, extrasaction="ignore"
        )
        writer.writeheader()
        for row in provenance_rows:
            writer.writerow(row)
    print("wrote", output_csv_path)

    # ------------------------------------------------------------------------------
    # Print an aligned summary to stdout focused on the saturation question: PMT V and
    # how much headroom each scan has to the 16-bit ceiling. Same PMT across gels plus
    # large headroom on all of them is the finding that rules a scan-setting artifact
    # out; a low headroom on the suspect alone would be the finding that rules it in.
    # ------------------------------------------------------------------------------
    print()
    print("Saturation / provenance summary (PMT voltage and 16-bit headroom):")
    header = "{:<9} {:>4} {:>7} {:>4} {:>7} {:>10} {:>10}  {}".format(
        "screen", "rep", "status", "PMT", "RangeHi", "headroom", "scale", "note"
    )
    print(header)
    print("-" * len(header))
    for row in provenance_rows:
        note = ""
        if row.get("status") == "NO_INF":
            note = "no .inf kept for this scan"
        elif row.get("status") == "READ_ERROR":
            note = "unreadable: " + str(row.get("detail", ""))[:60]
        elif (
            row.get("headroom_to_16bit_ceiling") is not None
            and row["headroom_to_16bit_ceiling"] < 3277
        ):  # within 5% of ceiling
            note = "WARNING: within 5% of 16-bit ceiling -- possible saturation"
        print(
            "{:<9} {:>4} {:>7} {:>4} {:>7} {:>10} {:>10}  {}".format(
                str(row.get("screen", "")),
                str(row.get("replicate", "")),
                str(row.get("status", "")),
                str(row.get("pmt_voltage_V", "")),
                str(row.get("range_high", "")),
                str(row.get("headroom_to_16bit_ceiling", "")),
                str(row.get("scale_type", "")),
                note,
            )
        )

    # A one-line verdict on the shared-PMT question, computed rather than asserted: are
    # all successfully-read scans at the same PMT voltage?
    # Count any row whose PMT was read, regardless of whether the filename id parsed:
    # the shared-sensitivity question is about the scan setting, not the filename, so
    # an ID_PARSE_INCOMPLETE row with a valid PMT still belongs in this verdict.
    read_pmts = [
        row.get("pmt_voltage_V")
        for row in provenance_rows
        if row.get("pmt_voltage_V") is not None
    ]
    print()
    if read_pmts and len(set(read_pmts)) == 1:
        print(
            "All {} readable scans share PMT voltage {} V -- no sensitivity "
            "difference between gels.".format(len(read_pmts), read_pmts[0])
        )
    elif read_pmts:
        print(
            "PMT voltages differ across scans: {} -- inspect which gel differs.".format(
                sorted(set(read_pmts))
            )
        )
    else:
        print("No PMT voltages were read; check the .inf paths above.")


if __name__ == "__main__":
    main()
