# design_template.py -- Annotated template for factorial experiment manifests
#
# Copy this file into your experiment folder as design.py, then edit.
# Generate the manifest with: lw exp manifest <id>
#
# The file must define a single dict named DESIGN. It is loaded via exec(),
# so standard Python (imports, variables, helper functions) is allowed.
# This is trusted user code -- no sandboxing.

DESIGN = {
    # ----------------------------------------------------------------
    # experiment_id (required)
    # Must match the experiment you are generating for.
    # The manifest command checks this and aborts on mismatch.
    # ----------------------------------------------------------------
    "experiment_id": "LM-0050",
    # ----------------------------------------------------------------
    # expected_samples (required)
    # Total number of samples: len(extras) + factorial_after_excludes.
    # The command aborts if the generated count does not match.
    # Calculate this by hand before running -- it catches design errors.
    #
    # Example: 2 extras + (5 strains x 2 DNAs x 3 salt - 4 excluded) = 28
    # ----------------------------------------------------------------
    "expected_samples": 28,
    # ----------------------------------------------------------------
    # extras (optional)
    # Non-factorial samples positioned first (positions 1, 2, ...).
    # Common uses: input lanes, size markers, positive/negative controls.
    #
    # Each entry is a dict. Required key: "label".
    # Optional keys: "strain", "dna", "amount", or any custom key.
    # All keys are stored in the conditions JSON column.
    #
    # If "strain" is provided, it must be registered in the database.
    # If "label" is provided in an extra, it overrides the strain label
    # from the database. This lets you describe the sample's role
    # (e.g., "input") rather than the strain identity.
    # If no "label" is given but "strain" is, the strain's database
    # label is used as fallback.
    # ----------------------------------------------------------------
    "extras": [
        {"label": "input", "strain": "LY100", "amount": "5% total"},
        {"label": "input", "strain": "LY456", "amount": "5% total"},
    ],
    # ----------------------------------------------------------------
    # categories (required)
    # Dict of category_name -> list of levels.
    # The factorial cross product is generated from all categories.
    #
    # Special category names:
    #   "strain" -- values are strain IDs (or None). Stored in
    #               samples.strain_id as a foreign key. Labels are
    #               resolved from the strains table.
    #   "dna"    -- values are free text. Stored in samples.dna.
    #
    # All other category names are stored in the conditions JSON column.
    #
    # Use None for no-protein controls. The manifest renders these
    # with strain_id = NULL and label = "no protein".
    #
    # Order within each list determines display order (after sorting
    # by sort_order). Put the most common/control condition first.
    # ----------------------------------------------------------------
    "categories": {
        "strain": ["LY100", "LY456", "LY789", "LY012", None],
        "dna": ["ARS1", "A-B2-"],
        "kglut": ["250", "300", "350"],
    },
    # ----------------------------------------------------------------
    # sort_order (required)
    # List of category names defining sort priority for factorial rows.
    # Must contain exactly the same keys as categories -- no more, no less.
    # The command validates this and aborts on mismatch.
    #
    # Samples are sorted by each category in order, using the position
    # of each value in its category list (not alphabetically).
    # So if sort_order is ["strain", "dna", "kglut"], all samples for
    # the first strain come first, then within that by dna order, etc.
    # ----------------------------------------------------------------
    "sort_order": ["strain", "dna", "kglut"],
    # ----------------------------------------------------------------
    # exclude (optional)
    # List of lambda functions. Each receives a dict of
    # {category_name: value} for one factorial combination.
    # Return True to drop that combination from the manifest.
    #
    # Common patterns:
    #
    # Drop no-protein controls at non-baseline salt concentrations:
    #   lambda r: r["strain"] is None and r["kglut"] != "250"
    #
    # Drop a specific strain + DNA combination:
    #   lambda r: r["strain"] == "LY789" and r["dna"] == "A-B2-"
    #
    # Drop all combinations for a strain (simpler to remove from
    # categories, but useful when toggling):
    #   lambda r: r["strain"] == "LY012"
    #
    # Errors in exclude lambdas are caught and reported with the
    # offending row, so you will see which combination caused the
    # problem.
    # ----------------------------------------------------------------
    "exclude": [
        lambda r: r["strain"] is None and r["kglut"] != "250",
    ],
}


# ----------------------------------------------------------------
# CALCULATING expected_samples
#
# Full cross product:
#   len(strain) x len(dna) x len(kglut) = 5 x 2 x 3 = 30
#
# Excluded combinations:
#   strain=None with kglut=300: 1 x 2 x 1 = 2
#   strain=None with kglut=350: 1 x 2 x 1 = 2
#   Total excluded: 4
#   (Only kglut != "250" excluded, so None+ARS1+250 and None+A-B2-+250
#    survive. Each non-250 kglut with dna=ARS1 and dna=A-B2- is dropped.)
#
#   Wait -- the exclude drops rows where strain is None AND kglut != "250".
#   None appears once in strain. kglut != "250" means kglut in {300, 350}.
#   So excluded rows = 1 (None) x 2 (dna) x 2 (non-250 kglut) = 4.
#
# Factorial after excludes: 30 - 4 = 26
# Extras: 2
# Total: 2 + 26 = 28
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# MINIMAL EXAMPLE (no extras, no excludes)
#
# DESIGN = {
#     "experiment_id": "LM-0010",
#     "expected_samples": 6,
#     "categories": {
#         "strain": ["LY100", "LY456", None],
#         "dna": ["ARS1", "A-B2-"],
#     },
#     "sort_order": ["strain", "dna"],
# }
#
# Produces 3 x 2 = 6 samples. No extras, no excludes.
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# NOTES
#
# - All strains in categories and extras must be registered with
#   `lw strain add` before running `lw exp manifest`.
#
# - Strains should have labels set (`lw strain update LY100 label "WT"`)
#   for readable manifests. If a strain has no label, its ID is used
#   as the label instead (no error, but harder to read).
#
# - If samples already exist for the experiment, the command aborts.
#   Use --force to delete existing samples and regenerate:
#       lw exp manifest LM-0050 --force
#
# - The manifest command writes manifest.csv to the experiment folder.
#   This CSV is for reference (loading into a spreadsheet, printing
#   for the bench). The database is the source of truth.
#
# - Positions are 1-indexed. Extras come first (1..N), then factorial
#   samples (N+1..M).
#
# - The conditions JSON column stores all category values for each
#   sample, enabling cross-experiment SQL queries like:
#       SELECT * FROM samples
#       WHERE json_extract(conditions, '$.kglut') = '300';
# ----------------------------------------------------------------
