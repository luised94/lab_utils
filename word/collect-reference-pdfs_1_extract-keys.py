#!/usr/bin/env python3
"""
Extract Zotero item keys from a .docx file.

Zotero stores citation data as JSON inside Word field codes (w:instrText elements).
Each citation contains URIs of the form:
    http://zotero.org/users/<userID>/items/<ITEMKEY>

This script unzips the .docx, parses word/document.xml, reconstructs fragmented
field codes, extracts item keys from URIs, and writes deduplicated keys to a file.

Usage:
    python3 collect-reference-pdfs_1_extract-keys.py <docx_path> [output_directory]

Output:
    <output_directory>/<docx_stem>_keys.txt   (one key per line, sorted)

Requires: python3 (standard library only - zipfile, xml.etree, json, re)
"""

import sys
import os
import zipfile
import tempfile
import xml.etree.ElementTree as ET
import json
import re

# =============================================================================
# CONFIGURATION
# =============================================================================

WORD_NAMESPACE = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
DOCUMENT_XML_PATH = "word/document.xml"
ZOTERO_FIELD_MARKER = "ZOTERO"
# As of Zotero 8, item keys are 8 alphanumeric characters.
# If this changes in a future version, the regex still captures them but
# the assertion below will warn about unexpected key lengths.
ITEM_KEY_PATTERN = re.compile(r"/items/([A-Za-z0-9]+)$")
CITATION_JSON_PATTERN = re.compile(r"\{.*\}", re.DOTALL)
EXPECTED_KEY_LENGTH = 8

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if len(sys.argv) < 2:
    print("Usage: python3 collect-reference-pdfs_1_extract-keys.py <docx_path> [output_directory]", file=sys.stderr)
    sys.exit(1)

docx_path = sys.argv[1]
output_directory = sys.argv[2] if len(sys.argv) >= 3 else os.getcwd()

if not os.path.isfile(docx_path):
    print(f"Error: File not found: {docx_path}", file=sys.stderr)
    sys.exit(1)

if not docx_path.endswith(".docx"):
    print(f"Error: Expected a .docx file, got: {docx_path}", file=sys.stderr)
    sys.exit(1)

if not os.path.isdir(output_directory):
    print(f"Error: Output directory not found: {output_directory}", file=sys.stderr)
    sys.exit(1)

docx_stem = os.path.splitext(os.path.basename(docx_path))[0]
output_keys_path = os.path.join(output_directory, docx_stem + "_keys.txt")

print(f"Input:  {docx_path}", file=sys.stderr)
print(f"Output: {output_keys_path}", file=sys.stderr)
print("", file=sys.stderr)

# =============================================================================
# EXTRACT document.xml FROM .DOCX
# =============================================================================

temporary_directory = tempfile.mkdtemp(prefix="zotero_extract_keys_")
document_xml_extracted_path = os.path.join(temporary_directory, DOCUMENT_XML_PATH)

with zipfile.ZipFile(docx_path, "r") as docx_zip:
    zip_contents = docx_zip.namelist()
    if DOCUMENT_XML_PATH not in zip_contents:
        print(f"Error: {DOCUMENT_XML_PATH} not found in .docx archive", file=sys.stderr)
        print(f"Archive contains: {zip_contents}", file=sys.stderr)
        sys.exit(1)
    docx_zip.extract(DOCUMENT_XML_PATH, temporary_directory)

document_xml_size = os.path.getsize(document_xml_extracted_path)
print(f"Extracted {DOCUMENT_XML_PATH} ({document_xml_size} bytes)", file=sys.stderr)

# =============================================================================
# PARSE XML - COLLECT fldChar AND instrText ELEMENTS IN DOCUMENT ORDER
# =============================================================================

tree = ET.parse(document_xml_extracted_path)
root = tree.getroot()

field_char_tag_suffix = "fldChar"
instr_text_tag_suffix = "instrText"
field_char_type_attribute = f"{{{WORD_NAMESPACE}}}fldCharType"

ordered_elements = []
for element in root.iter():
    tag_local = element.tag.split("}")[-1] if "}" in element.tag else element.tag

    if tag_local == field_char_tag_suffix:
        field_type = element.get(field_char_type_attribute, "")
        # Fallback: check all attributes for fldCharType (namespace variations)
        if not field_type:
            for attribute_name in element.attrib:
                if "fldCharType" in attribute_name:
                    field_type = element.attrib[attribute_name]
                    break
        ordered_elements.append(("field_char", field_type))

    elif tag_local == instr_text_tag_suffix:
        ordered_elements.append(("instr_text", element.text or ""))

print(f"Found {len(ordered_elements)} field_char + instr_text elements", file=sys.stderr)

# =============================================================================
# RECONSTRUCT FRAGMENTED FIELD CODES
# =============================================================================
# WHY: Word splits long strings across multiple <w:instrText> elements within
# a single field code run (delimited by fldChar begin/end). A Zotero citation
# JSON can be 500+ characters, so it's almost always fragmented. If we parsed
# individual instrText elements, we'd get broken JSON. Instead, we concatenate
# all instrText between each begin/end pair, then keep only those containing
# the ZOTERO marker.
#
# FRAGILITY: If Zotero moves away from field codes (e.g., to content controls
# or custom XML parts), this section will find 0 field codes. The assertion
# below catches this explicitly.

zotero_field_codes = []
current_buffer = None

for element_type, element_value in ordered_elements:
    if element_type == "field_char":
        if element_value == "begin":
            current_buffer = []
        elif element_value == "end" and current_buffer is not None:
            concatenated_text = "".join(current_buffer)
            if ZOTERO_FIELD_MARKER in concatenated_text:
                zotero_field_codes.append(concatenated_text)
            current_buffer = None
    elif element_type == "instr_text" and current_buffer is not None:
        current_buffer.append(element_value)

print(f"Reconstructed {len(zotero_field_codes)} Zotero field codes", file=sys.stderr)

# =============================================================================
# EXTRACT ITEM KEYS FROM CITATION JSON URIs
# =============================================================================

all_uris = []
parse_failures = 0

for field_code_text in zotero_field_codes:
    json_match = CITATION_JSON_PATTERN.search(field_code_text)
    if json_match is None:
        parse_failures += 1
        continue

    try:
        citation_data = json.loads(json_match.group())
    except json.JSONDecodeError:
        parse_failures += 1
        continue

    citation_items = citation_data.get("citationItems", [])
    for citation_item in citation_items:
        uris = citation_item.get("uris", [])
        all_uris.extend(uris)

unique_keys = set()
for uri in all_uris:
    key_match = ITEM_KEY_PATTERN.search(uri)
    if key_match:
        unique_keys.add(key_match.group(1))

sorted_keys = sorted(unique_keys)

# =============================================================================
# ASSERTIONS
# =============================================================================

if len(zotero_field_codes) == 0:
    print("ERROR: No Zotero field codes found in document.xml.", file=sys.stderr)
    print("  The Zotero Word plugin may have changed its storage format.", file=sys.stderr)
    print("  Run probe_01_docx.sh to inspect the .docx structure.", file=sys.stderr)
    sys.exit(1)

if len(sorted_keys) == 0:
    print("ERROR: Field codes found but no item keys extracted.", file=sys.stderr)
    print("  The URI format may have changed. Check citation JSON manually.", file=sys.stderr)
    sys.exit(1)

unexpected_key_lengths = {k: len(k) for k in sorted_keys if len(k) != EXPECTED_KEY_LENGTH}
if unexpected_key_lengths:
    print(f"WARNING: {len(unexpected_key_lengths)} keys have unexpected length (expected {EXPECTED_KEY_LENGTH}):", file=sys.stderr)
    for key, length in list(unexpected_key_lengths.items())[:5]:
        print(f"  {key} (length {length})", file=sys.stderr)
    print("  Proceeding anyway - keys may still be valid.", file=sys.stderr)

# =============================================================================
# WRITE OUTPUT
# =============================================================================

with open(output_keys_path, "w") as keys_file:
    for key in sorted_keys:
        keys_file.write(key + "\n")

# =============================================================================
# SUMMARY
# =============================================================================

print("", file=sys.stderr)
print("--- Summary ---", file=sys.stderr)
print(f"Zotero field codes:  {len(zotero_field_codes)}", file=sys.stderr)
print(f"Parse failures:      {parse_failures}", file=sys.stderr)
print(f"Total URIs:          {len(all_uris)}", file=sys.stderr)
print(f"Unique item keys:    {len(sorted_keys)}", file=sys.stderr)
print(f"Output written to:   {output_keys_path}", file=sys.stderr)

if parse_failures > 0:
    print(f"WARNING: {parse_failures} field codes could not be parsed", file=sys.stderr)

# =============================================================================
# CLEANUP
# =============================================================================

os.remove(document_xml_extracted_path)
# Remove the intermediate directories created by zipfile.extract
word_directory = os.path.join(temporary_directory, "word")
if os.path.isdir(word_directory):
    os.rmdir(word_directory)
os.rmdir(temporary_directory)
