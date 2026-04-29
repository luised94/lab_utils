# Collect Reference PDFs from a Word Document

Extract all Zotero-cited references from a `.docx` file and copy their PDFs into a single folder for sharing.

## The Problem

You have a Word document with Zotero citations and need to collect all the cited PDFs into one folder (e.g., to share with a collaborator or for AI-assisted review). Manually finding each PDF is tedious with 50-100 references.

## How It Works

The pipeline has three steps:

1. **Extract keys** - The Zotero Word plugin embeds citation metadata as JSON inside Word field codes (`<w:instrText>` elements in the `.docx` XML). Each citation contains a URI like `http://zotero.org/users/USERID/items/ITEMKEY`. The script unzips the `.docx`, parses the XML, reconstructs fragmented field codes, and extracts the 8-character item keys.

2. **Resolve paths** - Zotero's SQLite database (`zotero.sqlite`) maps item keys to attachment file paths. For linked files (managed by Attanger/ZotFile), paths are stored with an `attachments:` prefix relative to a base directory. The script queries the database, strips the prefix, and prepends the base attachment path to produce absolute file paths.

3. **Copy PDFs** - A PowerShell script (not `cp`) copies files to the target directory. PowerShell is required because Dropbox smart-sync stores files as cloud-only placeholders that must be "hydrated" (downloaded) before copying. PowerShell's file access triggers this hydration automatically. The script uses `-LiteralPath` throughout to handle Unicode characters (e.g., `ú`, `”`, `-`) and special characters (brackets, backticks) in filenames.

## Prerequisites

- **WSL** (Ubuntu) with `python3`, `sqlite3`, `unzip`
- **PowerShell** accessible via `powershell.exe` (WSL interop)
- **Zotero** with Better BibTeX and Attanger/ZotFile (linked file attachments)
- Zotero must be **closed** when running the pipeline (SQLite locking)

## Quick Start

```bash
# Full pipeline (dry run - no files copied)
bash collect-reference-pdfs_4_run-pipeline.sh \
    "/mnt/c/.../your_manuscript.docx" \
    "/mnt/c/.../output_pdfs" \
    --dry-run

# Full pipeline (actual copy)
bash collect-reference-pdfs_4_run-pipeline.sh \
    "/mnt/c/.../your_manuscript.docx" \
    "/mnt/c/.../output_pdfs"
```

Or run each step individually:

```bash
python3 collect-reference-pdfs_1_extract-keys.py  <docx_path> [output_dir]
bash    collect-reference-pdfs_2_resolve-paths.sh  <keys_file> [output_dir]
bash    collect-reference-pdfs_3_copy-pdfs.sh      <paths_tsv> <target_dir> [--dry-run]
```

## Machine-Specific Configuration

These paths are hardcoded and must be updated if your environment changes:

| Setting | Location | Current Value |
|---|---|---|
| Zotero SQLite | `_2_resolve-paths.sh` | `/mnt/c/Users/Luised94/Zotero/zotero.sqlite` |
| Zotero prefs.js | `_2_resolve-paths.sh` | `/mnt/c/Users/Luised94/AppData/Roaming/Zotero/Zotero/Profiles/*/prefs.js` |
| Fallback base path | `_2_resolve-paths.sh` | `/mnt/c/Users/Luised94/MIT Dropbox/Luis Martinez/zotero-storage` |

The base attachment path is read from `prefs.js` automatically when possible. The fallback is only used if `prefs.js` cannot be found or does not contain `extensions.zotero.baseAttachmentPath`.

## Output Files

| File | Description |
|---|---|
| `{stem}_keys.txt` | One Zotero item key per line (sorted, deduplicated) |
| `{stem}_keys_paths.tsv` | Columns: key, relative_path, absolute_path, exists_on_disk |
| `{stem}_keys_unresolved.tsv` | Columns: key, item_type, title - items cited but without a linked PDF |
| `_copy_results.tsv` | Per-file copy status: OK, SKIPPED, NOT_FOUND, HYDRATE_FAIL, COPY_FAIL |

## Known Limitations

- **Stale attachments**: Zotero may have multiple attachment entries for one item (e.g., from Attanger renames). The copy step deduplicates by file path, so this produces a NOT_FOUND for the stale entry but the real file still copies. These can be cleaned up in Zotero (right-click the broken attachment and delete).
- **No PDF attached**: Some cited items (books, websites, datasets) may not have a PDF. These appear in the unresolved file.
- **Zotero must be closed**: The SQLite database is locked while Zotero is running.
- **Field code format**: Tested with Zotero 8 on Word. If Zotero changes its Word integration method (e.g., switching from field codes to content controls), the extraction script will need updating - it will report 0 field codes found as a clear signal.

## Probe Scripts

The `probe_*.sh` scripts were used during development to verify assumptions about the `.docx` format and SQLite schema. They are read-only diagnostics and can be re-run to debug issues:

- `probe_00_preflight.sh` - Dependencies, file existence, Zotero process check
- `probe_01_docx.sh` - Inspects how Zotero stores citations in the `.docx` XML
- `probe_02_sqlite.sh` - Inspects SQLite schema and linked file path format

## Troubleshooting

**"0 Zotero field codes found"** - The `.docx` citation format may have changed. Run `probe_01_docx.sh` to see where Zotero is now storing citation data. Check `customXml/`, `settings.xml`, and content controls.

**"Keys resolved: 0"** - Zotero may be open (locking the database), or the item keys in the `.docx` don't match the local library (e.g., the document was created on a different machine). Run `probe_02_sqlite.sh` to verify.

**"HYDRATE_FAIL"** - Dropbox couldn't download the file. Check your internet connection and Dropbox sync status. The file may have been deleted from Dropbox.

**"NOT_FOUND"** - The attachment record exists in Zotero's database but the file is missing from disk. Usually a stale attachment from a rename. Check if a duplicate attachment was copied successfully.

**Wrong number of PDFs** - Some items have duplicate attachment entries. The copy step deduplicates by file path, but duplicate entries with different filenames (from renames) will both attempt to copy.
