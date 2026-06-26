[Thread 0 / planning]
- Delta-free ATPase framing chosen over TOST: no biologically defensible
  equivalence margin; shared catalytic mutation makes non-restoration
  mechanistically expected. Cite Lakens (TOST) as considered-and-rejected.
- nlme (not lme4, not hand-rolled) for mixed models: ships with R, citable
  (Pinheiro & Bates), dependency-light, reviewer-friendly. lme4 adds
  transitive deps; hand-rolling REML on unbalanced ATPase is indefensible.
- Holm (not Bonferroni) over pre-specified families; family scope per I2.
- source()-only path support; Rscript/interactive deliberately unsupported.

[Thread 1 appends here ...]
- C1: script directory detected by scanning sys.frame(i)$ofile over
  seq_len(sys.nframe()) and taking the last non-NULL (robust to nested
  source()); normalizePath()+dirname() gives SCRIPT_DIRECTORY. Absence of any
  ofile is the Rscript/interactive case -> clear stop(). Same inline block
  duplicated in all four scripts (procedural; no helper, per style contract).
- C1: per-script resolution chain = script-relative (SCRIPT_DIRECTORY/<file>)
  -> MC_DROPBOX_PATH/<orig relative path> -> stop() naming both absolute
  paths. Early hard-stop on empty MC_DROPBOX_PATH removed; emptiness is now
  tolerated and only fatal if the script-relative branch also misses.
- C1 (atpase-analyze): base directory chosen by probing the FIRST registry
  file's relative path under each candidate base; OUTPUT_DIRECTORY follows the
  chosen base (SCRIPT_DIRECTORY for Zenodo, consolidated_analysis for Dropbox).
  dir.create(OUTPUT_DIRECTORY) moved to after resolution.
- C1 (atpase-plot): processed_data.csv resolved script-relative ->
  MC_DROPBOX_PATH consolidated_analysis -> stop(); OUTPUT_DIRECTORY set to the
  resolved input's dirname so plot writes beside the data it read. getwd()
  dependency removed.
- C3: as.integer() coercion in the sequential-index assert confirmed to keep
  firing despite readxl returning imagej_index as double. Dimension (nrow/ncol)
  checks MOVED above the colnames<- assignment so an unexpected column count
  yields the clear assertion instead of an opaque length-mismatch error.
  Sheet numbering unchanged (xlsx sheetIndex= and readxl sheet= both 1-based).
- C4 RATIONALE CORRECTION: the handoff calls this a "guard-ordering crash"; it
  is actually a vacuously-PASSING (ineffective) guard -- loading_data[[col]] on
  a missing column returns NULL, is.infinite(NULL)=logical(0), any()=FALSE, so
  stopifnot passed without checking anything. Fix is unchanged (move derived-
  calc mutate above the guard); it makes the guard effective. No crash existed.
- C4: relocated the premature "...computed." message to after the mutate;
  removed the now-orphaned comment-only preamble.
- C5 (350mm): removed duplicate trailing message "Labels mapped and factor
  order applied." (redundant with "Factor ordering applied.").
- C5 (kglut): reworded premature message("Plotting completed...") (emitted
  before the global-baseline plot was built) to "Faceted-by-kglut plot
  constructed."; also normalized a non-ASCII +/- glyph in a comment.
- C0d: nlme restored/recorded via renv (renv::install("nlme")+snapshot());
  no script edit. nlme is a recommended package shipping with R, so this is a
  record step for the deposit library, not a build.

[Thread 2 appends here ...]
[Thread 3 appends here ...]
