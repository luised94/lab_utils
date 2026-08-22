# Load-bearing files, and what is archived

One screen, kept current by hand. The point of this note is to stop someone
moving, renaming, or deleting a load-bearing file six months from now without
realising the pipeline (or the regression harness) depends on it. If you move
anything in the LIVE list, run the harness immediately: red means you moved
something load-bearing.

    uv run tests/run_pipeline_tests.py     # from images/; expect all steps PASS

ASCII only, to match the repo.

===============================================================================
LIVE - do NOT move or rename without updating every referrer
===============================================================================

Data-side pipeline (CSV in, per-sample numbers out), in order:

    export_lane_profiles.ijm    ImageJ macro; the LIVE producer of the seed input
                                manual_lane_profiles.csv. Not Python, not exercised
                                by the harness, but the front of the whole pipeline.
                                Nothing downstream runs without its output.
    validate_sample_sheet.py    gate the sheet against the profiles
    analyze_gel.py              detrend, detect bands, measure -> band_measurements
                                (the manual measure stage; contains measure_gel.py's
                                math inlined verbatim, per the no-sharing rule)
    extract_lane_values.py      region/band mode -> per-lane values
    plot_single_experiment.py   blank-correct -> single_experiment values
    validate_manifest.py        gate the replicate manifest
    aggregate_repeats.py        stack replicates -> aggregate

Picker (live, server-based):

    serve_gel_picker.py         local picker server; POST /extract shells the
                                extractor. Computes the extractor path as a SIBLING
                                (pathlib __file__ .parent / extract_lane_values.py),
                                so serve_gel_picker.py and extract_lane_values.py
                                MUST stay co-located. Moving one without the other
                                breaks serve even if the harness step still launches.
    gel_picker.html             the picker page; fetch('/payload'), POST selection.
                                serve computes its path as a sibling too.

Regression harness:

    tests/run_pipeline_tests.py         the harness (drives the 7 scripts + serve)
    tests/make_replicate_fixture.py     one-shot authoring tool for the n>1 fixture
    tests/golden_values.json            captured goldens + seed-fixture hash
    tests/fixtures/                      seed inputs (two gels + manifest +
                                         failure_modes/malformed_sample_sheet.csv)

Anything the harness names is load-bearing BY TEST: aggregate_repeats,
analyze_gel, extract_lane_values, plot_single_experiment, validate_manifest,
validate_sample_sheet, serve_gel_picker, gel_picker.html.

===============================================================================
LIVE but NOT harness-tested - move with care, the harness will not catch it
===============================================================================

    make_payload.py             builds a static, server-free copy of the picker
                                (payload.json) for a frozen, shareable page. An
                                OFFLINE alternative to serve, not part of the live
                                serve+picker path (the live picker fetches /payload
                                from the server, not the file). Kept as a live
                                utility; the harness does not exercise it, so a
                                break here will not turn the harness red.

===============================================================================
ARCHIVED - superseded, kept for reference/reuse, NOT run (see archive/README.md)
===============================================================================

    archive/measure_gel.py              old AUTOMATED measure stage (Stage 2).
                                        Superseded by the ImageJ macro + analyze_gel
                                        (manual path). Its math lives on, inlined in
                                        analyze_gel.py; analyze_gel's comments cite
                                        it at archive/measure_gel.py.
    archive/validate_gel_image.py       old TIFF-quality validator (Stage 1).
    archive/prepare_gel_image.py        old image prep / path normalisation (Stage 0).
    archive/inventory_tiff_tags.py      TIFF tag dumper used by the old Stage 1.
    archive/sensitivity_harness.py      the OLD harness (drove Stage 0->1->2 by
                                        subprocess). Superseded by
                                        tests/run_pipeline_tests.py.
    archive/validate_manifest_sheet.py  old manifest validator. Strict subset of
                                        validate_manifest.py (older schema, no
                                        is_example support, no template writing).

Nothing live imports or execs any archived file. analyze_gel.py MENTIONS
archive/measure_gel.py in comments only (provenance of copied math), never at
runtime. archive/ was reached with git mv, so history is preserved; recovery is
`git log --follow archive/<file>`.
