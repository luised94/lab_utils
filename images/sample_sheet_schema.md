# Sample sheet schema

One row per lane, for one gel. Authored in Excel, exported as sample_sheet.csv
into the gel's <stem>_gel_analysis directory. The repo copy of this file is the
source of truth; a copy pasted into an Excel tab is a convenience mirror only.

Validated by validate_sample_sheet.py, which hard-fails on structural or
self-contradictory sheets and warns on soft issues. Run it after the ImageJ export
macro has produced manual_lane_profiles.csv in the same directory:

    uv run validate_sample_sheet.py <gel_analysis_dir>

The directory (or any file inside it) is the only argument in the common case; the
sample sheet and profile CSV are found inside it by their standard names.

## The four column families

Columns are grouped by what varies them and when. Identity and disposition are
required on every sheet. Experimental class is required. Provenance is required.
Biological factor columns are permitted but omitted when they do not apply (this
is normal: factors are per-experiment).

### Identity (required)

    lane_index    Image position, left to right, 1..comb. Matches the profile CSV.
                  Assigned by the ImageJ macro from ROI position; means nothing
                  biological. You do not choose this; it is the image coordinate.
    well_number   Load order / sample-sheet order, a permutation of 1..comb, each
                  once. This is the stable identity of a sample. YOU author this to
                  match how you loaded the gel. If you loaded left-to-right it
                  equals lane_index; if you loaded in reverse it is comb+1-lane_index.

The validator reports the observed relationship (identity, full_inversion, or
permutation). It never applies a flip. Read that line every run: a gel you meant
to invert that reports "permutation" means a well_number was mistyped.

### Disposition (required)

    lane_content    What was physically loaded.
                    Vocabulary: sample | ladder | empty
    analysis_role   What the analysis does with the lane.
                    Vocabulary: reference | measured | excluded

Legal combinations:
    sample  -> analysis_role in {reference, measured}
    ladder  -> analysis_role = excluded
    empty   -> analysis_role = excluded
Exactly one lane in the sheet has analysis_role = reference, and it must be a
sample. The reference is the lane you normalize against; pick one input.

### Experimental class (required)

    condition_type  The design category of the lane.
                    Vocabulary: input | experiment | positive_control |
                                negative_control | not_applicable

Legal combinations (this is the impossible-combination guard):
    sample          -> any of input | experiment | positive_control | negative_control
    ladder or empty -> not_applicable   (a non-analyte lane has no design category)

Positive and negative controls are distinct values, neither is required, and more
than one of each is allowed. Their absence is a WARNING, not a failure, because a
clean control is not always achievable, but you should notice when it is missing.

### Provenance (required)

    sample_label    Human name of the sample (e.g. MCM_46_001). On sample rows it
                    must be a real label. On ladder/empty rows it must be the
                    reserved token not_applicable.
    prep_source     Where the material came from: a prep id (LM-0008) or a year for
                    pre-schema material (2015). Same rule: real on samples,
                    not_applicable on non-samples.
    notes           Free text. Optional; may be blank.

Why not_applicable instead of blank on non-sample rows: blank is ambiguous (did you
mean "does not apply" or "forgot to fill it?"). The reserved token not_applicable
means "deliberately nothing here", so the validator can treat a blank on those rows
as a mistake and fail it. Use the token, not an empty cell, and not a descriptive
word like "ladder".

### Biological factors (optional, add when relevant)

Independent variables you analyze or plot against. Add columns only when the
experiment varies them. Examples: temperature_c, ph, buffer, salt_mM,
mutation_present, suppressor_present, dna_sequence, protein_concentration_nM.
When present they should use controlled vocabularies or unit-separated values, so a
later group-by does not fragment on "300 mM" vs "300mM". Absent columns are fine;
the validator does not require them and does not invent defaults for them.

## Not in this sheet

Technical replicate is a per-gel fact, not a per-lane one, so it lives in the
manifest (the replicate column), not here. One sample sheet describes one gel.

## Authoring tips

- Put Excel data-validation dropdowns on lane_content, analysis_role, and
  condition_type so a typo cannot be entered; the validator is the backstop.
- Author well_number explicitly (type it, or fill-down a formula then eyeball it),
  and confirm the reported relationship line matches what you intended.
- Cross-check the load order against the physical image and your labelled tubes
  before exporting.
