# Deposit finalization - renv lockfile + bare-run prerequisites

## renv lockfile
- The deposit library must capture exactly the packages the scripts load:
  readxl, tidyverse (ggplot2/dplyr/tidyr/readr/...), and nlme. base `stats`
  needs no entry. nlme was restored/recorded in Thread 1 (C0d) so its consumers
  (C9 loading diagnostic, C13 kglut interaction, C17 ATPase models) had it
  before this final step; this commit only VERIFIES and freezes.
- Finalize:
    renv::snapshot()        # writes/refreshes renv.lock from the loaded library
    renv::status()          # MUST report the project in sync, nothing missing
- Confirm nlme is present and pinned:
    grep -A3 '"nlme"' renv.lock   # expect a version + repository record
- nlme and stats ship with R, so a clean install resolves them from the R
  distribution / CRAN; readxl removed the old Java/rJava (xlsx) dependency in
  Thread 1 (C3), so no system Java is required for the bare run.

## Clean-room run-through (do this in a fresh checkout)
1. renv::restore()
2. source the scripts in dependency order (see run commands below).
3. Confirm each prints "=== SCRIPT COMPLETE ===" and that the in-script
   assertions (the project's tests, per I1) did not stop().

## KNOWN BARE-RUN GOTCHA (flagged, NOT silently changed) - loading scripts only
- The two LOADING scripts place an Arial/Cairo FONT GATE ABOVE the data load.
  In a font-less clean room they stop BEFORE any compute, so loading stats never
  run. This is pre-existing Thread-1 infrastructure; the ATPase scripts in this
  thread have NO such gate and run font-free.
- Two options, both explicit (we do NOT change behavior silently):
  (a) PREREQUISITE (recommended for the deposit as-is): install a font providing
      Arial (e.g. ttf-mscorefonts-installer, or map to Liberation Sans) before
      running the loading scripts. Document this as a run prerequisite.
  (b) ONE-LINE RELOCATION (only if a maintainer chooses): move the font-gate
      block to immediately ABOVE the loading scripts' PLOT section (below all
      stats compute + stats-CSV writes). Then loading STATISTICS run font-free
      and only plotting requires the font. This is a real behavior change to the
      failure mode and must be made deliberately by a maintainer, with a
      DECISION_LOG entry - it is FLAGGED here, not performed.

## Files the deposit should contain after this thread
- 4 R scripts (2 loading, 2 ATPase) + renv.lock + .Rprofile/renv/ activate.
- STATISTICAL_METHODS.md, STATISTICAL_SELF_STUDY.md, DEPLOYMENT_NOTES.md,
  DECISION_LOG.md (consolidated).
- Data CSVs co-located for the script-relative resolution path.
