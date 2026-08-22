# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy>=2.5.2",
#     "tifffile>=2026.8.16",
# ]
# ///
r"""
Serve the gel band/region picker for one gel, locally, so the final friction step
is: point at the gel analysis directory, open the browser, choose a band or a
span, and the per-lane CSV is written into that directory by the same extractor
the command line uses.

    uv run serve_gel_picker.py '<gel_analysis_dir>'
    # then open http://localhost:8080

On startup it reads manual_lane_profiles.csv (the width-summed, plate-background-
subtracted profile plus each lane's ROI) and band_measurements.csv (lane identity,
the consensus band regions, per-lane apex peaks, and the fragility / occupancy /
smile metrics), and rasterises a display-scaled crop of the rotated TIFF for the
background. All of that is served as one JSON payload the page fetches once.

Finalize in the page POSTs the current selection; the server runs
extract_lane_values.py (a separate script, invoked as a subprocess, never
imported, so the number the browser writes is byte-for-byte the command-line
number) and returns the written CSV path and its rows.

Standard library HTTP server (no web framework). tifffile+numpy read the gel; the
PNG is encoded with stdlib zlib. Flat and procedural; emit_message / die are the
only helpers, duplicated verbatim per the house convention.
"""

import argparse
import base64
import csv
import datetime
import http.server
import json
import pathlib
import struct
import subprocess
import sys
import threading
import zlib

import numpy
import tifffile

# Display convention shared with measure_gel.py's overlays: dark bands on white,
# contrast set from robust percentiles so the MinIsWhite scan reads correctly.
DISPLAY_PERCENTILE_LOW = 1.0
DISPLAY_PERCENTILE_HIGH = 99.5
# The crop is downsampled to this width before PNG encoding; it is a background
# reference behind the profiles, not a measurement surface, so nearest-neighbour
# at a modest width keeps the payload small without misleading anyone.
GEL_CROP_TARGET_WIDTH_PIXELS = 480
EIGHT_BIT_MAXIMUM = 255
DEFAULT_PORT = 8080

PROFILE_CSV_FILENAME = "manual_lane_profiles.csv"
BAND_MEASUREMENTS_FILENAME = "band_measurements.csv"

ACCUMULATED_RUN_LOG_LINES = []


def emit_message(source_tag, message_text):
    timestamp_text = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    tagged_line = "[" + source_tag + "] " + message_text
    ACCUMULATED_RUN_LOG_LINES.append(timestamp_text + " " + tagged_line)
    print(tagged_line, file=sys.stderr)


def die(source_tag, message_text):
    emit_message(source_tag, "ERROR: " + message_text)
    sys.exit(1)


# =============================================================================
# Arguments and path resolution (directory-addressed, like the other scripts).
# =============================================================================

argument_parser = argparse.ArgumentParser(
    prog="serve_gel_picker.py",
    description="Serve the band/region picker for one gel and write the chosen CSV.",
    allow_abbrev=False,
)
argument_parser.add_argument(
    "gel_path", help="the gel analysis directory, or any file inside it"
)
argument_parser.add_argument(
    "--tiff",
    dest="tiff_override",
    default=None,
    help="path to the rotated TIFF; by default a single .tif in the dir or its parent",
)
argument_parser.add_argument("--port", type=int, default=DEFAULT_PORT)
argument_parser.add_argument("--profile-csv", dest="profile_csv_override", default=None)
argument_parser.add_argument(
    "--band-measurements", dest="band_measurements_override", default=None
)
parsed_arguments = argument_parser.parse_args()

given_path = pathlib.Path(parsed_arguments.gel_path)
if given_path.is_dir():
    gel_analysis_directory = given_path
elif given_path.is_file():
    gel_analysis_directory = given_path.parent
else:
    die("input", "path does not exist: " + str(given_path))

profile_csv_path = (
    pathlib.Path(parsed_arguments.profile_csv_override)
    if parsed_arguments.profile_csv_override
    else gel_analysis_directory / PROFILE_CSV_FILENAME
)
band_measurements_path = (
    pathlib.Path(parsed_arguments.band_measurements_override)
    if parsed_arguments.band_measurements_override
    else gel_analysis_directory / BAND_MEASUREMENTS_FILENAME
)
extractor_script_path = (
    pathlib.Path(__file__).resolve().parent / "extract_lane_values.py"
)
picker_html_path = pathlib.Path(__file__).resolve().parent / "gel_picker.html"
for required_path in (
    profile_csv_path,
    band_measurements_path,
    extractor_script_path,
    picker_html_path,
):
    if not required_path.is_file():
        die("input", "missing required file: " + str(required_path))

# Locate the rotated TIFF the ROIs were drawn on. Priority: an explicit override;
# then the export provenance the ImageJ macro wrote (authoritative); then the
# pipeline convention <stem>.tif beside the <stem>_gel_analysis directory. A blind
# *.tif glob is deliberately NOT used: an experiment folder holds several gels'
# TIFFs, and the alphabetically-first one is usually the wrong gel (a leading
# "2026.07.10..." sorts before "20260818..." because '.' precedes '0').
GEL_ANALYSIS_DIRECTORY_SUFFIX = "_gel_analysis"
expected_tiff_stem = gel_analysis_directory.name
if expected_tiff_stem.endswith(GEL_ANALYSIS_DIRECTORY_SUFFIX):
    expected_tiff_stem = expected_tiff_stem[: -len(GEL_ANALYSIS_DIRECTORY_SUFFIX)]
expected_tiff_filename = expected_tiff_stem + ".tif"

provenance_image_directory = None
provenance_image_title = None
provenance_path = gel_analysis_directory / "manual_export_provenance.txt"
if provenance_path.is_file():
    for provenance_line in provenance_path.read_text(encoding="utf-8-sig").splitlines():
        stripped_line = provenance_line.strip()
        lowered_line = stripped_line.lower()
        # Accept "key: value", "key = value", "key<tab>value", or "key value".
        if lowered_line.startswith("image_title"):
            provenance_image_title = stripped_line[len("image_title") :].lstrip(" \t:=")
        elif lowered_line.startswith("image_directory"):
            provenance_image_directory = stripped_line[len("image_directory") :].lstrip(
                " \t:="
            )

if parsed_arguments.tiff_override:
    tiff_path = pathlib.Path(parsed_arguments.tiff_override)
    if not tiff_path.is_file():
        die("input", "--tiff not found: " + str(tiff_path))
else:
    tiff_path = None
    tiff_search_candidates = []
    if provenance_image_directory and provenance_image_title:
        tiff_search_candidates.append(
            pathlib.Path(provenance_image_directory) / provenance_image_title
        )
    # The pipeline keeps <stem>.tif as a sibling of <stem>_gel_analysis, but also
    # check inside the analysis dir in case a copy lives there.
    tiff_search_candidates.append(
        gel_analysis_directory.parent / expected_tiff_filename
    )
    tiff_search_candidates.append(gel_analysis_directory / expected_tiff_filename)
    for candidate_path in tiff_search_candidates:
        if candidate_path.is_file():
            tiff_path = candidate_path
            break
    # Last resort: a TIFF whose name contains this gel's stem, never a blind first
    # match across unrelated gels in the same folder.
    if tiff_path is None:
        for search_directory in (gel_analysis_directory, gel_analysis_directory.parent):
            stem_matched_tiffs = sorted(
                found_tiff
                for found_tiff in search_directory.glob("*.tif")
                if expected_tiff_stem in found_tiff.name
            )
            if stem_matched_tiffs:
                tiff_path = stem_matched_tiffs[0]
                break
    if tiff_path is None:
        emit_message(
            "input",
            "no TIFF matching '" + expected_tiff_filename + "' found beside the gel "
            "directory; serving without the gel crop (pass --tiff to set it explicitly)",
        )

# =============================================================================
# Read the two CSVs (utf-8-sig strips the Excel BOM; every cell trimmed) so the
# lane/band joins cannot break invisibly on a hand-exported file.
# =============================================================================

with band_measurements_path.open(newline="", encoding="utf-8-sig") as band_csv_handle:
    band_measurement_rows = [
        {
            column_name: (
                cell_value.strip() if isinstance(cell_value, str) else cell_value
            )
            for column_name, cell_value in raw_row.items()
        }
        for raw_row in csv.DictReader(band_csv_handle)
    ]
with profile_csv_path.open(newline="", encoding="utf-8-sig") as profile_csv_handle:
    profile_rows = [
        {
            column_name: (
                cell_value.strip() if isinstance(cell_value, str) else cell_value
            )
            for column_name, cell_value in raw_row.items()
        }
        for raw_row in csv.DictReader(profile_csv_handle)
    ]

gel_id = band_measurement_rows[0]["gel_id"]

# Profiles and ROIs, grouped by lane, ordered by migration. raw_value is already
# the plate-background-subtracted width sum, so integrate-span sums it directly.
profile_values_by_lane_index = {}
roi_by_lane_index = {}
migration_millimetres_by_pixel = {}
plate_background_median = None
for profile_row in profile_rows:
    lane_index = int(profile_row["lane_index"])
    profile_values_by_lane_index.setdefault(lane_index, []).append(
        (
            int(profile_row["migration_position_pixels"]),
            int(round(float(profile_row["raw_value"]))),
        )
    )
    if lane_index not in roi_by_lane_index:
        roi_by_lane_index[lane_index] = {
            "roi_x": int(profile_row["roi_x"]),
            "roi_y": int(profile_row["roi_y"]),
            "roi_w": int(profile_row["roi_w"]),
            "roi_h": int(profile_row["roi_h"]),
        }
    migration_millimetres_by_pixel[int(profile_row["migration_position_pixels"])] = (
        float(profile_row["migration_position_millimetres"])
    )
    if plate_background_median is None:
        plate_background_median = float(profile_row["plate_background_median"])

sorted_lane_indices = sorted(profile_values_by_lane_index)
profiles_payload = {}
for lane_index in sorted_lane_indices:
    profile_values_by_lane_index[lane_index].sort(
        key=lambda migration_sample: migration_sample[0]
    )
    profiles_payload[str(lane_index)] = [
        summed_signal
        for (_migration_pixels, summed_signal) in profile_values_by_lane_index[
            lane_index
        ]
    ]

migration_row_count = len(profiles_payload[str(sorted_lane_indices[0])])
maximum_migration_pixel = max(migration_millimetres_by_pixel)
millimetres_per_pixel = (
    migration_millimetres_by_pixel[maximum_migration_pixel] / maximum_migration_pixel
)
roi_y_origin = roi_by_lane_index[sorted_lane_indices[0]]["roi_y"]

# Lane identity, from band_measurements (a ladder/empty lane absent there gets no
# identity; the page shows it as excluded rather than inventing one).
identity_by_lane_index = {}
for measurement_row in band_measurement_rows:
    lane_index = int(measurement_row["lane_index"])
    if lane_index not in identity_by_lane_index:
        identity_by_lane_index[lane_index] = {
            "well": int(measurement_row["well_number"]),
            "label": measurement_row["sample_label"],
            "role": measurement_row["analysis_role"],
        }

# Consensus band metadata (constant per band index) and the per-lane per-band
# apex + fragility used for the overlays and the consistency panel.
band_metadata_by_index = {}
for measurement_row in band_measurement_rows:
    band_index = int(measurement_row["canonical_band_index"])
    if band_index not in band_metadata_by_index:
        band_metadata_by_index[band_index] = {
            "index": band_index,
            "mm": float(measurement_row["canonical_position_millimetres"]),
            "window_start_pixels": int(measurement_row["window_start_pixels"]),
            "window_end_pixels": int(measurement_row["window_end_pixels"]),
            "occupancy": int(measurement_row["band_occupancy"]),
            "supported": measurement_row["band_is_well_supported"] == "yes",
            "spread_mm": float(measurement_row["cross_lane_spread_millimetres"]),
            "fragile_lanes": 0,
        }
    if measurement_row["baseline_agreement_status"] == "fragile":
        band_metadata_by_index[band_index]["fragile_lanes"] += 1
bands_payload = [
    band_metadata_by_index[band_index] for band_index in sorted(band_metadata_by_index)
]

per_lane_per_band_payload = {}
for measurement_row in band_measurement_rows:
    lane_index = int(measurement_row["lane_index"])
    band_index = int(measurement_row["canonical_band_index"])
    per_lane_per_band_payload.setdefault(str(lane_index), {})[str(band_index)] = {
        "area": float(measurement_row["reported_area"]),
        "apex": float(measurement_row["apex_height_above_baseline"]),
        "fragile": measurement_row["baseline_agreement_status"] == "fragile",
        "apex_pixels": int(measurement_row["apex_migration_pixels"]),
        "basis": measurement_row["reported_value_basis"],
    }

# =============================================================================
# Rasterise a display-scaled crop of the gel spanning all lane ROIs, and encode
# it as a PNG data URI with stdlib zlib (no Pillow/matplotlib dependency).
# =============================================================================

roi_left_edge_pixels = min(roi["roi_x"] for roi in roi_by_lane_index.values())
roi_right_edge_pixels = max(
    roi["roi_x"] + roi["roi_w"] for roi in roi_by_lane_index.values()
)
gel_crop_data_uri = None
if tiff_path is not None and tiff_path.is_file():
    gel_image = tifffile.imread(str(tiff_path))
    display_minimum = numpy.percentile(gel_image, DISPLAY_PERCENTILE_LOW)
    display_maximum = numpy.percentile(gel_image, DISPLAY_PERCENTILE_HIGH)
    gel_crop = gel_image[
        roi_y_origin : roi_y_origin + migration_row_count,
        roi_left_edge_pixels:roi_right_edge_pixels,
    ].astype(numpy.float32)
    display_scaled = numpy.clip(
        (gel_crop - display_minimum) / (display_maximum - display_minimum), 0.0, 1.0
    )
    # gray_r: invert so bands (high signal) read dark on a white gel.
    gray_reversed = (EIGHT_BIT_MAXIMUM * (1.0 - display_scaled)).astype(numpy.uint8)
    downsample_scale = gray_reversed.shape[1] / GEL_CROP_TARGET_WIDTH_PIXELS
    target_height_pixels = int(round(gray_reversed.shape[0] / downsample_scale))
    row_indices = numpy.linspace(
        0, gray_reversed.shape[0] - 1, target_height_pixels
    ).astype(int)
    column_indices = numpy.linspace(
        0, gray_reversed.shape[1] - 1, GEL_CROP_TARGET_WIDTH_PIXELS
    ).astype(int)
    downsampled = gray_reversed[numpy.ix_(row_indices, column_indices)]
    # Build an 8-bit grayscale PNG: each scanline prefixed with filter byte 0.
    raw_scanlines = bytearray()
    for scanline in downsampled:
        raw_scanlines.append(0)
        raw_scanlines.extend(scanline.tobytes())
    compressed_image_data = zlib.compress(bytes(raw_scanlines), 9)
    png_bytes = bytearray(b"\x89PNG\r\n\x1a\n")
    image_header = struct.pack(
        ">IIBBBBB", downsampled.shape[1], downsampled.shape[0], 8, 0, 0, 0, 0
    )
    for chunk_tag, chunk_data in (
        (b"IHDR", image_header),
        (b"IDAT", compressed_image_data),
        (b"IEND", b""),
    ):
        png_bytes += struct.pack(">I", len(chunk_data))
        png_bytes += chunk_tag + chunk_data
        png_bytes += struct.pack(">I", zlib.crc32(chunk_tag + chunk_data) & 0xFFFFFFFF)
    gel_crop_data_uri = "data:image/png;base64," + base64.b64encode(
        bytes(png_bytes)
    ).decode("ascii")
    emit_message(
        "payload",
        "gel crop %dx%d encoded from %s"
        % (downsampled.shape[1], downsampled.shape[0], tiff_path.name),
    )
else:
    emit_message("payload", "no TIFF resolved; serving without the gel crop")

payload = {
    "gel_id": gel_id,
    "gel_dir": str(gel_analysis_directory),
    "migration_row_count": migration_row_count,
    "millimetres_per_pixel": millimetres_per_pixel,
    "roi_y_origin": roi_y_origin,
    "roi_left_edge_pixels": roi_left_edge_pixels,
    "roi_right_edge_pixels": roi_right_edge_pixels,
    "plate_background_median": plate_background_median,
    "lanes": sorted_lane_indices,
    "profiles": profiles_payload,
    "roi": {
        str(lane_index): roi_by_lane_index[lane_index]
        for lane_index in sorted_lane_indices
    },
    "identity": identity_by_lane_index,
    "bands": bands_payload,
    "per_lane_per_band": per_lane_per_band_payload,
    "gel_crop": gel_crop_data_uri,
}
payload_json_bytes = json.dumps(payload).encode("utf-8")
picker_html_bytes = picker_html_path.read_bytes()

# =============================================================================
# HTTP server. GET / serves the page, GET /payload the data, POST /extract runs
# the extractor and returns the written CSV, POST /shutdown stops the server.
# =============================================================================


class GelPickerRequestHandler(http.server.BaseHTTPRequestHandler):
    def log_message(self, message_format, *message_arguments):
        # Route the access log through emit_message (stderr) rather than the
        # handler's own stdout writer, so stdout stays clean.
        emit_message("http", (message_format % message_arguments))

    def _respond(self, status_code, content_type, body_bytes):
        self.send_response(status_code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body_bytes)))
        self.end_headers()
        self.wfile.write(body_bytes)

    def do_GET(self):
        if self.path == "/" or self.path == "/index.html":
            self._respond(200, "text/html; charset=utf-8", picker_html_bytes)
        elif self.path == "/payload":
            self._respond(200, "application/json", payload_json_bytes)
        else:
            self._respond(404, "text/plain", b"not found")

    def do_POST(self):
        request_length = int(self.headers.get("Content-Length", 0))
        request_body = self.rfile.read(request_length) if request_length else b"{}"
        if self.path == "/shutdown":
            self._respond(200, "application/json", b'{"ok":true}')
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return
        if self.path != "/extract":
            self._respond(404, "text/plain", b"not found")
            return
        selection = json.loads(request_body)
        # Build the exact command line; the extractor writes into the gel dir.
        extractor_command = [
            sys.executable,
            str(extractor_script_path),
            str(gel_analysis_directory),
        ]
        if selection.get("mode") == "band":
            extractor_command += [
                "--band",
                str(int(selection["band_index"])),
                "--quantity",
                selection.get("quantity", "area"),
            ]
        elif selection.get("mode") == "region":
            extractor_command += [
                "--region",
                ("%g" % float(selection["start_millimetres"])),
                ("%g" % float(selection["end_millimetres"])),
            ]
            if selection.get("net_baseline") in ("none", "straight"):
                extractor_command += ["--net-baseline", selection["net_baseline"]]
        else:
            self._respond(
                400, "application/json", b'{"ok":false,"error":"bad selection"}'
            )
            return
        completed = subprocess.run(extractor_command, capture_output=True, text=True)
        if completed.returncode != 0:
            self._respond(
                200,
                "application/json",
                json.dumps(
                    {
                        "ok": False,
                        "error": completed.stderr.strip().splitlines()[-1]
                        if completed.stderr.strip()
                        else "extractor failed",
                        "command": extractor_command,
                    }
                ).encode("utf-8"),
            )
            return
        # Recover the written CSV: the extractor names it deterministically from
        # the selection, so re-derive the stem and read it back for the page.
        if selection["mode"] == "band":
            output_stem = "extract_band_%d_%s" % (
                int(selection["band_index"]),
                selection.get("quantity", "area"),
            )
        else:
            output_stem = "extract_region_%g-%gmm_%s" % (
                float(selection["start_millimetres"]),
                float(selection["end_millimetres"]),
                selection.get("net_baseline", "none"),
            )
        written_csv_path = gel_analysis_directory / (output_stem + ".csv")
        with written_csv_path.open(
            newline="", encoding="utf-8-sig"
        ) as written_csv_handle:
            written_rows = list(csv.DictReader(written_csv_handle))
        self._respond(
            200,
            "application/json",
            json.dumps(
                {
                    "ok": True,
                    "csv_path": str(written_csv_path),
                    "command": " ".join(extractor_command),
                    "rows": written_rows,
                }
            ).encode("utf-8"),
        )


http_server = http.server.ThreadingHTTPServer(
    ("127.0.0.1", parsed_arguments.port), GelPickerRequestHandler
)
emit_message("serve", "gel %s" % gel_id)
emit_message(
    "serve", "open http://localhost:%d  (Ctrl-C to stop)" % parsed_arguments.port
)
try:
    http_server.serve_forever()
except KeyboardInterrupt:
    emit_message("serve", "stopped")
