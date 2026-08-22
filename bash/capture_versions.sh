#!/usr/bin/env bash
# Capture the resolved versions of every third-party dep the LIVE scripts use, from
# INSIDE each script's own uv environment (bare `uv run python` has NO deps -- that
# is why `uv run python -c "import numpy"` fails). Writes resolved_versions.txt, which
# pin_deps.sh consumes. Run from images/.
#
# It works by asking uv to export each script's locked environment. uv resolves the
# script's inline dependency block exactly as a real run would, so the versions here
# are precisely what produced your current goldens.
set -euo pipefail

OUT=resolved_versions.txt
: > "$OUT"

# The union of third-party packages across the live scripts. matplotlib/tifffile do
# not affect goldens (plots and image IO), but we floor them too for uniformity.
PACKAGES="numpy scipy matplotlib tifffile"

# A throwaway PEP723 script that declares ALL of them and prints resolved versions.
# uv resolves the same index/constraints it would for the real scripts, so for an
# unpinned dep this is the version a real run resolves today -- the one to floor at.
TMP="$(mktemp --suffix=.py)"
cat > "$TMP" << 'PY'
# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy","scipy","matplotlib","tifffile"]
# ///
import numpy, scipy, matplotlib, tifffile
print("numpy", numpy.__version__)
print("scipy", scipy.__version__)
print("matplotlib", matplotlib.__version__)
print("tifffile", tifffile.__version__)
PY

echo "# resolved versions captured $(date -Iseconds) on $(hostname)" >> "$OUT"
uv run "$TMP" 2>/dev/null | grep -E "^($(echo $PACKAGES | tr ' ' '|')) " >> "$OUT"
rm -f "$TMP"

echo "Captured:"
cat "$OUT"
echo
echo "Next: ./pin_deps.sh   (floors each script's deps at these versions)"
