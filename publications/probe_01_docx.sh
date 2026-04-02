#!/usr/bin/env bash
# =============================================================================
# PROBE 01: .DOCX CITATION STORAGE FORMAT
# =============================================================================
# Destructive operations: NONE.
#   - Reads .docx (zip) and extracts to /tmp. Original file untouched.
#   - All output to stdout.
#
# Purpose: Determine how Zotero 8 stores citation data in this Word file.
#
# Assumptions tested:
#   A1. Citations live in field codes (<w:instrText> with "ZOTERO_ITEM").
#   A2. Citation JSON has citationItems[].uris[] containing item keys.
#   A3. URIs match pattern: http://zotero.org/users/.../items/ITEMKEY
#
# Alternative locations checked:
#   - word/document.xml    (field codes or content controls)
#   - customXml/           (custom XML parts - newer method)
#   - word/settings.xml    (document variables)
#   - word/footnotes.xml   (citations in footnotes)
#
# Expected outcome: Zotero data found in at least one location.
#   Output reveals storage method, URI format, and sample item keys.
# =============================================================================
set -euo pipefail

DOCX_PATH="/mnt/c/Users/Luised94/MIT Dropbox/Luis Martinez/Lab/publications-and-presentations/lemr_publication_bypass_orc4r/manuscript_orc4r_bypass_v16.docx"
PROBE_DIR=$(mktemp -d /tmp/zotero_probe_docx_XXXXXX)

echo "========================================"
echo "PROBE 01: .DOCX CITATION FORMAT"
echo "========================================"
echo "Temp dir: ${PROBE_DIR}"
echo ""

# --- Extract ---
echo "--- 1a. Archive file listing ---"
unzip -l "$DOCX_PATH" 2>/dev/null | tail -n +4 | head -30
echo "  ..."
echo ""

unzip -o -q "$DOCX_PATH" -d "$PROBE_DIR"

# --- Search for Zotero strings across all files ---
echo "--- 1b. Zotero string search ---"
for term in "ZOTERO" "zotero.org" "CSL_CITATION" "citationItems"; do
    HITS=$(grep -rl "$term" "$PROBE_DIR" 2>/dev/null || true)
    if [ -n "$HITS" ]; then
        while IFS= read -r f; do
            COUNT=$(grep -c "$term" "$f" 2>/dev/null || echo 0)
            echo "  '${term}' -> ${f#"$PROBE_DIR/"} (${COUNT}x)"
        done <<< "$HITS"
    else
        echo "  '${term}' -> not found"
    fi
done
echo ""

# --- Parse with Python ---
echo "--- 1c. Citation data extraction ---"
python3 - "$PROBE_DIR" << 'PYEOF'
import xml.etree.ElementTree as ET
import json, re, sys, os

extract_dir = sys.argv[1]
doc_path = os.path.join(extract_dir, "word", "document.xml")
ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}

if not os.path.isfile(doc_path):
    print("word/document.xml not found - unexpected .docx structure")
    sys.exit(0)

print(f"document.xml size: {os.path.getsize(doc_path)} bytes\n")

tree = ET.parse(doc_path)
root = tree.getroot()

# ---- Method 1: Reconstruct fragmented field codes ----
# Word splits long field codes across multiple w:instrText elements
# between fldChar begin/end pairs. We concatenate them.
elements = []
for elem in root.iter():
    tag = elem.tag.split('}')[-1] if '}' in elem.tag else elem.tag
    if tag == 'fldChar':
        ftype = elem.get(f'{{{ns["w"]}}}fldCharType', '')
        if not ftype:
            for attr in elem.attrib:
                if 'fldCharType' in attr:
                    ftype = elem.attrib[attr]
                    break
        elements.append(('fld', ftype))
    elif tag == 'instrText':
        elements.append(('instr', elem.text or ''))

fields = []
buf = None
for etype, val in elements:
    if etype == 'fld':
        if val == 'begin':
            buf = []
        elif val == 'end' and buf is not None:
            full = ''.join(buf)
            if 'ZOTERO' in full:
                fields.append(full)
            buf = None
    elif etype == 'instr' and buf is not None:
        buf.append(val)

print(f"Reconstructed ZOTERO field codes: {len(fields)}")

# ---- Extract keys from field codes ----
all_uris = []
sample_shown = False
for field in fields:
    m = re.search(r'\{.*\}', field, re.DOTALL)
    if not m:
        continue
    try:
        parsed = json.loads(m.group())
    except json.JSONDecodeError:
        continue

    if not sample_shown and 'citationItems' in parsed:
        ci = parsed['citationItems'][0]
        print(f"\nSample citationItem keys: {list(ci.keys())}")
        if 'uris' in ci:
            print(f"Sample URI: {ci['uris']}")
        if 'itemData' in ci:
            d = ci['itemData']
            print(f"itemData keys: {list(d.keys())}")
            print(f"  title: {str(d.get('title',''))[:80]}")
            print(f"  author: {d.get('author', [None])[:1]}")
            print(f"  issued: {d.get('issued','')}")
        sample_shown = True

    for ci in parsed.get('citationItems', []):
        all_uris.extend(ci.get('uris', []))

keys = set()
for uri in all_uris:
    m = re.search(r'/items/([A-Za-z0-9]+)$', uri)
    if m:
        keys.add(m.group(1))

print(f"\nTotal URIs: {len(all_uris)}")
print(f"Unique item keys: {len(keys)}")
if keys:
    sample = sorted(keys)[:5]
    print(f"Sample keys: {sample}")
    key_lengths = set(len(k) for k in keys)
    print(f"Key lengths: {key_lengths}")

# ---- Method 2: Check customXml ----
cxml_dir = os.path.join(extract_dir, "customXml")
if os.path.isdir(cxml_dir):
    print(f"\ncustomXml/ contents: {os.listdir(cxml_dir)}")
    for f in os.listdir(cxml_dir):
        fp = os.path.join(cxml_dir, f)
        if os.path.isfile(fp):
            with open(fp, 'r', errors='replace') as fh:
                content = fh.read()
            if 'zotero' in content.lower():
                print(f"  {f}: contains Zotero data ({len(content)} chars)")
                print(f"  preview: {content[:200]}")
else:
    print(f"\ncustomXml/ not found")

# ---- Method 3: Check settings.xml for doc variables ----
settings = os.path.join(extract_dir, "word", "settings.xml")
if os.path.isfile(settings):
    with open(settings, 'r', errors='replace') as fh:
        sc = fh.read()
    if 'zotero' in sc.lower():
        print(f"\nsettings.xml: contains Zotero data")
        matches = re.findall(r'<[^>]*[Zz]otero[^>]*>.*?</', sc, re.DOTALL)
        for m in matches[:2]:
            print(f"  {m[:200]}")
    else:
        print(f"\nsettings.xml: no Zotero data")

# ---- Write keys to file for probe_02 ----
keys_file = os.path.join(extract_dir, "zotero_item_keys.txt")
with open(keys_file, 'w') as f:
    for k in sorted(keys):
        f.write(k + '\n')
print(f"\nKeys written to: {keys_file}")

PYEOF

echo ""
echo "Temp dir preserved for probe_02: ${PROBE_DIR}"
echo "Cleanup when done: rm -rf ${PROBE_DIR}"
