You are implementing Thread 3 -- the FINAL thread -- of a multi-thread plan to
add statistical analysis to R scripts for an ORC4R suppressor screen (PLOS ONE
publication, Zenodo deposit). Threads 1 (infrastructure), 2a (350mm loading),
and 2b (kglut loading) are DONE and verified against real data. You do the
ATPase statistics AND the two documentation deliverables. The full plan and
global invariants are in the attached IMPLEMENTATION_HANDOFF.md -- read it first,
in full. Then read the attached CONSOLIDATED DECISION_LOG (Threads 0+1+2a+2b) --
it is the source material for the docs and carries decisions you must honor and
consolidate.

Attached:
- IMPLEMENTATION_HANDOFF.md                       (plan -- authoritative)
- DECISION_LOG.md                                 (CONSOLIDATED 0+1+2a+2b --
  read ALL of it; the docs are built by consolidating it)
- 260411_orc4r-screen_analyze-sofa-atpase.R       (ATPase compute; updated by
  Thread 1's xlsx->readxl migration -- you edit this)
- 260411_orc4r-screen_plot-sofa-atpase.R          (ATPase plot -- you edit this)
- the ATPase output CSVs from the Thread-1 run    (combined_raw_data,
  processed_data, summary_data -- real data shape, esp. the UNBALANCED blocks)

YOUR SCOPE: commits C15, C16, C17, C18a, C18b (ATPase statistics), then
C19 (STATISTICAL_METHODS.md), C20 (self-study reference), C21 (renv lockfile
finalize + bare-run notes). This is the last thread; after it the deposit is
complete.

STYLE CONTRACT (non-negotiable): procedural, no helpers, no abstractions,
explicit descriptive names, "nothing too clever." Dependencies: base R stats
+ nlme only, nothing hand-rolled. Plots show EXACT p-values, never stars.

GLOBAL INVARIANTS (full text in handoff/log; the ones that bite here):
- I1 "green" = script source()s clean AND its assertions pass on real data.
- I2 Holm family = WT-contrasts within a SINGLE timepoint. NEVER pooled across
     timepoints. State it in each test.
- I4 n=3 HONESTY: reuse the consolidated-log caveat wording. Normality
     untestable at the per-cell replicate counts; Wilcoxon companions (where
     used) DECORATIVE; justification = blocked design + assay track record.
- I5 Mixed models that stand in for a reported paired test are diagnostics.
     NOTE: C17's pooled mixed model is NOT a stand-in diagnostic -- per the
     review delta D2 it is the PRIMARY home for the 4R-proximity claim. Do not
     mislabel it diagnostic-only.

CARRY-FORWARD (binding, from prior threads):
1. VACUOUS-GUARD LESSON (Thread 1): df[["missing_col"]] returns NULL and
   any()/is.* on NULL passes silently. EVERY guard naming a column MUST first
   stopifnot("col" %in% names(df)) before testing values.
2. DATA SHAPE you must verify (do not assume): ATPase blocks are
   UNBALANCED and INCOMPLETE -- WT and 4R ride along in nearly every experiment
   (high replicate counts); each mutant appears in only its ~3 experiments. The
   response is percent_adp_corrected (adp/(adp+atp), then minus No_ORC t=90
   background). Confirm from the attached processed_data CSV: the real per-
   sample replicate counts, the timepoints present (expect 0/15/45/90), the
   experiment/block column name, and which samples co-occur in which blocks.

*** ATPase NEGATIVE-VALUE HANDLING (user decision -- DIFFERS FROM LOADING ON
PURPOSE; encode exactly this) ***
Loading FLOORED negatives because they were background-subtraction artifacts
that ALSO broke fold/ratio metrics via near-zero denominators. ATPase does NOT
floor. Rationale (put in the log + methods doc): percent_adp_corrected is
analyzed by MIXED MODEL + per-timepoint DIFFERENCES -- both subtraction-based,
negative-safe, NO division by a near-zero quantity. Small negatives are EXPECTED
noise around a true zero in the near-zero 4R/suppressor group (the script
already documents at its background-subtraction step that suppressors can read
below background), and that noise is DATA the model must absorb. Flooring would
erase exactly the variance the model needs and would bias the near-zero group's
mean upward (which, unlike loading, is NOT cleanly one-directional safe here).
Therefore implement:
  (a) DETECT + REPORT: after background subtraction, scan percent_adp_corrected
      < 0; PRINT every offending row (sample, replicate, timepoint, raw value);
      write them to an output CSV. This is transparency, not handling.
  (b) KEEP RAW: negatives flow UNCHANGED into the mixed model and the
      per-timepoint differences. NO floor, NO new floored column, NO pmax.
  (c) SANITY TRIPWIRE (not a floor gate): stopifnot that no negative is
      implausibly large (choose + document a bound, e.g. percent_adp_corrected
      > -0.15, mirroring the existing WT-t0 background sanity threshold). A wild
      negative signals an upstream sign/load error and must STOP; small near-
      zero scatter passes silently. Membership-checked before value-checked.
  (d) DIVISION-COLUMN CHECK: confirm (the analyze script was scanned in
      planning and has NO fold-over-4R / percent-of-WT column -- the only
      division is adp/(adp+atp), bounded [0,1] BEFORE subtraction, so it is
      safe). VERIFY this still holds in the file you receive; if any division-
      based descriptive column exists, mark it descriptive-only and keep it out
      of every claim. Do NOT introduce one.

DELTA-FREE NULL FRAMING (handoff F5 / log Thread 0): the ATPase claim is the
conjunction of (a) "not restored to WT" = WT-vs-suppressor SIGNIFICANT
(rejectable; the load-bearing half) and (b) "near the 4R floor across all
timepoints" = (suppressor - 4R) difference + CI, POOLED across timepoints via
the mixed model (NOT a p>0.05 equivalence claim; justified by all suppressors
sharing 4R's catalytic mutation). No equivalence margin delta. TOST was
considered and rejected (cite Lakens in the docs).

COMMITS (one at a time; tight reporting; verify against real data):
- C15 Assumption diagnostics. Per-cell {n, mean, sd} + pooled/grouped residual
  Shapiro (DECORATIVE, I4). EXPLICITLY note the unequal-n -> variance
  assumptions and which tests are valid. Emit to CSV.
- C16 Per-timepoint WT-contrasts (PRIMARY/plotted; the "not restored to WT"
  half). Block = experiment. *** WT MUST BE SUBSET to each suppressor's SHARED
  blocks (review delta D1): a "WT vs Orc5" paired test uses only the
  experiments where BOTH appear; WT's full replicate count does NOT strengthen
  these. REPORT per-contrast n explicitly. *** Holm within-timepoint (I2).
  Because blocks are unbalanced/incomplete, confirm the pairing per contrast
  and state usable n; if a suppressor x timepoint cell is missing, say so.
- C17 Pooled mixed model (PRIMARY home for 4R-proximity, D2 -- NOT a
  diagnostic). nlme: percent_adp_corrected ~ sample * timepoint,
  random = ~1|experiment. Report: the suppressor-vs-4R contrast + CI POOLED
  across timepoints (the 4R-proximity result); the sample:timepoint interaction
  (do curves diverge -- WT pulls away over time); the sample main effect.
  timepoint as FACTOR (accumulation saturates, not linear) -- default; document
  the numeric/trend alternative, do not adopt without justification. Wrap
  lme + intervals() in tryCatch -> non-convergence is a reported finding, not a
  crash. EXPECT (carry-forward from loading's CI reality): suppressor-vs-4R CIs
  may be tight near zero (supporting "near floor") while WT-vs-suppressor is
  clearly significant; report honestly whatever the numbers are.
- C18a analyze script: write a results CSV (per-timepoint WT-contrasts +
  per-timepoint n; pooled mixed-model terms incl. suppressor-vs-4R contrast+CI).
- C18b plot script: READ the results CSV (do NOT recompute -- anti-drift) and
  annotate per-timepoint EXACT p above each cluster (WT contrasts). The plot
  script currently RECOMPUTES the summary independently; switch its stats source
  to the C18a CSV so the two scripts cannot drift. EXACT p, never stars.

DOCUMENTATION (consolidation, not reconstruction -- built from the log):
- C19 STATISTICAL_METHODS.md (Zenodo-facing, terse/defensive). One section per
  analysis (350mm, kglut, ATPase): design, block structure, test, assumption
  status at n=3, Holm family scope (I2), the delta-free justification (ATPase),
  the normalization firewall and WHICH column feeds WHICH claim (kglut), the
  negative-handling decisions AND why they DIFFER between loading (floor) and
  ATPase (no floor). Describe what the code ACTUALLY did, citing the real run
  results already in the log (e.g. kglut salt-sensitivity interaction n.s.,
  the firewall preventing a denominator-confounded claim). Cite Holm 1979;
  Pinheiro & Bates (nlme); Lakens (TOST, considered-rejected); a blocked-design
  reference.
- C20 Self-study reference (expansive, user-facing). Conceptual background; the
  affirming-the-null / absence-of-evidence discussion; equivalence/TOST as the
  considered-and-rejected alternative and WHY; the normalization-firewall
  cautionary tale (per-kglut WT-pinning -> degenerate interaction; how the
  global baseline fixes it; how the WT denominator confounded the naive "ORC4R
  collapses at 350" reading); a self-test question bank; the primary sources to
  read. This one may be richer/longer than the methods doc.
- C21 Finalize + verify the deposit renv lockfile (nlme captured). Note the
  bare-run gotcha already in the log: the loading scripts' Arial/Cairo font gate
  sits ABOVE data load, so a font-less clean-room run stops before any compute;
  document the prerequisite (or, if in scope and low-risk, note the one-line
  relocation that would let stats run font-free -- but do NOT silently change
  behavior; flag it).

NO DELIBERATE RED THIS THREAD. The only sanctioned red in the whole project was
C13 (Thread 2b, done). Everything here lands green.

DECISION_LOG: append under "[Thread 3 appends here]" every non-obvious choice
(per-timepoint test details, the WT shared-block subsetting and per-contrast n,
the negative-detection bound, the C17 modeling choices and the real
suppressor-vs-4R + interaction results, any convergence finding, doc structure
decisions). Output the FULL consolidated log at the end.

WHAT I CANNOT DO: I (the user) run the scripts, not you. You cannot execute to
confirm green; reason explicitly about why each sources clean, then give me at
the end: the two edited ATPase scripts, the two doc files, the lockfile notes,
the exact run commands (analyze BEFORE plot), and the exact outputs to paste
back (console trail incl. the negative-detection table and per-contrast n's; the
augmented + results CSVs; the C17 pooled suppressor-vs-4R contrast/CI and the
interaction term; confirmation the PDF rendered with exact p-values).

If anything in the handoff or log is wrong against the ACTUAL ATPase data or
code -- STOP and raise it before editing. Data and code are ground truth. In
particular: if the block structure is MORE incomplete than "each mutant in ~3
experiments" (e.g. a suppressor missing an entire timepoint), that changes which
per-timepoint contrasts are computable and must be surfaced before C16.

Start with your read-back: the actual ATPase data shape from the processed_data
CSV (per-sample replicate counts, timepoints, experiment/block column, which
samples co-occur in which blocks), the full list of negative percent_adp_corrected
rows with counts, confirmation there is no division-based descriptive column,
and your plan for the C16 WT shared-block subsetting. Then proceed commit by
commit.
