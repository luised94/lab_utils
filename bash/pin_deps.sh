#!/usr/bin/env bash
# Floor every live script's third-party deps at the versions in resolved_versions.txt
# (written by capture_versions.sh). Floor (>=), not pin (==): uv keeps resolving the
# latest compatible, so you are always testing the newest and the harness tells you
# the moment a numeric golden drifts; a >= floor also never becomes unbuildable the
# way an == pin eventually can. Run from images/, AFTER capture_versions.sh.
#
# After running, VERIFY:  uv run tests/run_pipeline_tests.py   (expect 13/13, NO re-bless)
set -euo pipefail

VERSIONS_FILE=resolved_versions.txt
[ -f "$VERSIONS_FILE" ] || { echo "run ./capture_versions.sh first"; exit 1; }

# Load name->version.
declare -A VERSION
while read -r name ver; do
  case "$name" in \#*|"") continue;; esac
  VERSION["$name"]="$ver"
done < "$VERSIONS_FILE"

# Each live script and the third-party deps it declares (stdlib-only scripts skipped).
python3 - "$VERSIONS_FILE" << 'PY'
import re, sys, pathlib

versions = {}
for line in pathlib.Path(sys.argv[1]).read_text().splitlines():
    line = line.strip()
    if not line or line.startswith("#"): continue
    name, ver = line.split()
    versions[name] = ver

scripts = [
    "analyze_gel.py", "plot_single_experiment.py", "aggregate_repeats.py",
    "serve_gel_picker.py",
    # stdlib-only (empty deps) -- nothing to floor, left untouched:
    # extract_lane_values.py, validate_sample_sheet.py, validate_manifest.py
]

for script in scripts:
    path = pathlib.Path(script)
    text = path.read_text()
    m = re.search(r"(# dependencies = \[)(.*?)(\])", text, re.S)
    if not m:
        print("SKIP (no dep block):", script); continue
    body = m.group(2)
    names = re.findall(r'"([A-Za-z0-9_.\-]+)', body)
    if not names:
        print("SKIP (empty deps):", script); continue
    new_lines = []
    for name in names:
        if name not in versions:
            print("WARN: no captured version for %s in %s; leaving unfloored" % (name, script))
            new_lines.append('#     "%s",' % name)
        else:
            new_lines.append('#     "%s>=%s",' % (name, versions[name]))
    new_block = m.group(1) + "\n" + "\n".join(new_lines) + "\n# " + m.group(3)
    text = text[:m.start()] + new_block + text[m.end():]
    path.write_text(text)
    print("floored %-26s -> %s" % (script, ", ".join("%s>=%s" % (n, versions.get(n,"?")) for n in names)))
PY

echo
echo "Done. VERIFY now:"
echo "    uv run tests/run_pipeline_tests.py     # expect 13/13, and NO re-bless lines"
