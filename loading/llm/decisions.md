================================================================================
DECISION_LOG (append-only) -- CONSOLIDATED (Threads 0 + 1 + 2a + 2b + 3)
================================================================================
[Thread 0 / planning -- scope, framing, and binding constraints]

-- ASSAY ROLES (ground truth carried by all threads) --
- Three datasets/analyses: (1) kglut titration loading = WT/ORC4R/+4sofa only,
  n=3, salts 250/300/350; (2) 350 mM all-suppressors loading = WT, ORC4R, all
  five suppressors, n=3, single salt; (3) ATPase timecourse = WT, 4R, all five
  suppressors, timepoints 0/15/45/90, UNEQUAL replicates (WT/4R high, mutants 2-3).
- Loading restorers = Orc1, Orc4. Reducers = Orc3, Orc5, Orc6. 4R < WT in
  loading. Grouping is comment-level only; all suppressors get identical
  statistical treatment.

-- ATPase NULL FRAMING: DELTA-FREE (no equivalence margin) --
- The claim "suppressors do not restore ATPase activity" is a claim FOR the
  null, which a non-significant test cannot establish; n=2-3 has near-zero power.
- TOST CONSIDERED and REJECTED: no biologically defensible margin delta; at small
  n TOST also lacks power. Cite Lakens.
- Chosen framing: conjunction of (a) "not restored to WT" = WT-vs-suppressor
  SIGNIFICANT (rejectable, load-bearing); (b) "near the 4R floor" =
  (suppressor - 4R) difference + CI, POOLED across timepoints, descriptive,
  justified by all suppressors sharing 4R's catalytic mutation.

-- TERMINOLOGY --
- "non-restoration-to-WT contrast" = WT vs suppressor; rejectable.
- "4R-floor-proximity estimate" = (suppressor - 4R) diff + CI; ATPase ONLY.
- In LOADING, (suppressor - 4R/ORC4R) is an ORDINARY effect size, NOT the floor
  claim -- reducers below 4R BY DESIGN. Narratives never conflated.
- "salt-sensitivity interaction" = salt x genotype in the kglut mixed model.
- "block" = replicate (loading) / experiment (ATPase).

-- DESIGN: LOADING AND ATPase PAIRED/BLOCKED --
- WT, 4R/ORC4R, mutant(s) co-run within a block. Pairing exploited. ATPase blocks
  UNBALANCED and INCOMPLETE (WT/4R ride along; each mutant in ~its own 2-3 blocks).

-- TOOLING: base R stats + nlme ONLY; nothing hand-rolled --
- t.test/wilcox.test/aov/lm/p.adjust(holm)/shapiro.test/confint/residuals.
- Mixed models via nlme (ships w/ R, citable, dependency-light). lme4 rejected
  (transitive deps); hand-rolled REML rejected. renv locks nlme.

-- MULTIPLE COMPARISONS: Holm, pre-specified families --
- FAMILY SCOPE (I2): {WT-contrasts within a SINGLE salt} or {WT-contrasts within
  a SINGLE timepoint}. NEVER across salts/timepoints.

-- n=3 HONESTY (I4) --
- Normality untestable at n=3 (Shapiro near-powerless; pass != license).
  Justify parametric tests by paired/blocked design + assay track record.
  Diagnostics emitted, labelled. Wilcoxon companions DECORATIVE (floor ~0.25).

-- ANNOTATION: exact p-values on plots, never stars (binding) --

-- PER-ANALYSIS PRIMARY FRAMINGS --
- 350 mM: WT-contrasts (paired t + Wilcoxon + Holm) non-restoration;
  (suppressor-ORC4R) ordinary effect size + CI; mixed-model PRE-COMMITTED
  DIAGNOSTIC (I5).
- kglut: BOTH framings -- per-salt paired WT-vs-ORC4R tests; salt x genotype
  interaction (the salt-sensitivity result).
- ATPase: per-timepoint PRIMARY/plotted for WT-difference; pooled mixed model
  PRIMARY for 4R-proximity.

-- ADVERSARIAL-REVIEW DELTAS (binding) --
- D1: ATPase per-timepoint WT-contrasts subset WT to shared blocks; report n.
- D2: pooled mixed model is PRIMARY home for "near 4R across timepoints".
- D3: 350 mM mixed model is a DECLARED diagnostic; paired t pre-committed.
- D4 (I3 FIREWALL): kglut interaction MUST use global baseline; value-checked
  assertion fires on the per-kglut column.
- D5: loading (X-4R) effect size is NOT the ATPase floor claim; mandatory comment.
- D6 (I2): Holm family scope fixed globally.
- D7: nlme restore extracted into C0d.

-- PATHS + INVOCATION --
- source()-only; Rscript/interactive unsupported. Chain: script-relative ->
  MC_DROPBOX_PATH -> stop() naming both. xlsx->readxl removes Java for bare run.

-- DOCS (generated LAST, by consolidation) --
- STATISTICAL_METHODS.md (terse/defensive) + self-study reference (expansive).

-- ATPase NEGATIVE-VALUE HANDLING (decided after kglut) --
- Loading FLOORED negatives (broke fold/ratio via near-zero denominators). ATPase
  does NOT floor: mixed model + per-timepoint DIFFERENCES (subtraction-based,
  negative-safe). DETECT+REPORT every negative; keep RAW; envelope/tripwire
  assertion only (stop on implausibly large negative). Verify no division-based
  descriptive column is fed by negatives.

[Thread 0 / planning]
- Delta-free over TOST (no defensible margin; shared catalytic mutation). Cite Lakens.
- nlme (not lme4, not hand-rolled). Holm over pre-specified families. source()-only.

[Thread 1]
- C1: script dir via last non-NULL sys.frame()$ofile; normalizePath()+dirname();
  Rscript/interactive -> clear stop(). Per-script chain script-relative ->
  MC_DROPBOX_PATH -> stop() naming both. Empty MC_DROPBOX_PATH tolerated unless
  script-relative also misses. atpase-analyze base chosen by probing first
  registry file; OUTPUT_DIRECTORY follows chosen base; dir.create after
  resolution. atpase-plot resolves processed_data.csv script-relative ->
  MC_DROPBOX_PATH consolidated_analysis -> stop(); OUTPUT_DIRECTORY = input dirname.
- C3: as.integer() in sequential-index assert still fires despite readxl double;
  dim checks MOVED above colnames<- to give clear assertion. Sheet numbering
  unchanged (1-based).
- C4: actually a vacuously-PASSING guard (NULL column -> any() FALSE), not a
  crash; fix (move derived mutate above guard) makes guard effective. Premature
  message relocated.
- C5: removed duplicate trailing message (350mm); reworded premature plotting
  message + normalized non-ASCII glyph (kglut).
- C0d: nlme restored/recorded via renv (install + snapshot); no script edit.

[Thread 2a -- 350 mM]
- C7: paired t-test (block=replicate) WT vs each of {ORC4R,+1,+3,+4,+5,+6sofa} on
  percent_wildtype; estimate WT-group; decorative paired Wilcoxon.
- Holm family = the SIX WT-contrasts within 350 mM.
- C8: (suppressor-ORC4R) ordinary paired effect size + 95% CI; mandatory comment
  it is NOT the ATPase floor argument; reducers below ORC4R by design.
- C9: lme(percent_wildtype ~ label, ~1|replicate) PRE-COMMITTED diagnostic;
  WT zero-variance violates homoscedasticity (why diagnostic); tryCatch ->
  non-convergence flagged, not crash.
- C6: per-label {n,mean,sd,var} + pooled within-label residual Shapiro (W as
  number, font-independent). Bartlett/Levene OMITTED (WT zero-variance).
- VERIFIED: block=replicate (1/2/3); 7 labels; n=3 each; WT=100 exactly
  (zero variance); pairing COMPLETE; all percent_wildtype >=0 in 350mm aggregate.
  WARNING: kglut full-data CSV has a negative (+5sofa rep3 350mM).

[Thread 2b -- kglut]
- Test genotypes WT/ORC4R/+4sofa only, n=3 x salt, balanced/complete; rotating
  n=1 suppressors filtered out (F4).
- ONE negative (+5sofa rep3 350mM), in a filtered-out n=1 suppressor -> floor INERT
  for all reported values.
- Negatives FLOORED-after-detection (loading decision); envelope assert
  (MAX_COUNT=3/36, MAX_ABS_PERCENT=50) STOPS on breach; floored into NEW columns
  (raw retained); logged.
- Global baseline from floored values; for test genotypes equals unfloored.
- FIREWALL degeneracy confirmed: per-kglut within-WT sd across salts = 0; global
  WT declines (rep means 100/76.4/43.1).
- C11 diagnostics mirror C6 + per-salt split; Bartlett/Levene omitted.
- C12 per-salt paired t on percent_wildtype_floored, {WT-vs-ORC4R, WT-vs-+4sofa,
  ORC4R-vs-+4sofa}; Holm within salt over the two WT-contrasts; ORC4R-vs-+4sofa
  raw, outside family.
- C13 interaction lme(percent_wildtype_global_floored ~ label*kglut, ~1|replicate),
  salt as FACTOR; trend alternative documented not adopted. Firewall assertion
  membership-then-value (sd(WT across salts)>1) then name check. SANCTIONED RED
  fired on per-kglut, green on global. RESULT: gap-change 250->350 small, wide CI
  crossing zero, NON-significant -- naive "ORC4R collapses at 350" largely WT's
  own decline.
- C14a CSVs carry normalization_column lineage. C14b plot per-salt WT-vs-ORC4R
  exact p + interaction p in caption.

[Thread 3 -- ATPase statistics + documentation]

-- DATA-SHAPE VERIFICATION (against processed_data.csv / summary_data.csv) --
- Block column = experiment_label. NO replicate column in ATPase (block IS the
  experiment); negative report uses experiment_label as the block id.
- Analysed samples (exclude No_ORC, WT_DNA): WT, 4R, Orc1, Orc3, Orc4, Orc5, Orc6.
- Timepoints 0/15/45/90 present for EVERY analysed sample; no whole-timepoint
  cell missing. 156 analysed rows (39 x 4).
- Replicate (lane) counts: WT 13, 4R 12, Orc1 2, Orc3 3, Orc4 3, Orc5 3, Orc6 3.
- DISCREPANCIES vs handoff (surfaced before C16, none blocking):
  (1) Orc1 in only 2 experiments (n=2 blocks), NOT ~3 -> WT-vs-Orc1 is n=2/df=1.
  (2) Orc4 in 2 experiments but 3 LANES (two Orc4 lanes in 2020_08_24, one WT
      lane). Decision: block-mean aggregation for C16 (experiment = paired unit)
      -> WT-vs-Orc4 = n=2 block-pairs. C17 keeps all 3 lanes (random intercept
      absorbs the nesting) -- a reason C17 is primary.
  (3) 2021_07_29d is a low-signal experiment: WT reads all zeros there; it is a
      shared Orc5 block, so WT-vs-Orc5 includes a 0-0 pair at some timepoints
      (weakens, not breaks).
- NEGATIVE CENSUS: 33 rows < 0, max magnitude 0.0334 (2020_08_13b Orc4 t15
  = -0.03341); all in near-floor groups + two controls (No_ORC 2020_07_13 t0,
  WT_DNA 2020_08_13b t0, WT 2020_08_13b t0). All >> tripwire -0.15 -> passes.
- DIVISION-COLUMN CHECK: only division is adp/(adp+atp), bounded [0,1] BEFORE
  subtraction; no fold-over-4R / percent-of-WT column exists; none introduced.

-- C15 DIAGNOSTICS --
- Per-cell {n, mean, sd} (sample x timepoint) -> atpase_diagnostics_cellstats.csv.
- Pooled + per-timepoint within-cell-residual Shapiro (DECORATIVE, I4 caveat
  reused) -> atpase_diagnostics_shapiro.csv. Bartlett/Levene OMITTED (zero-
  variance cells undefined). Unequal-n + gross variance heterogeneity (WT large,
  near-floor cells ~0/exactly 0) noted explicitly -> homoscedastic-residual
  assumption of C17 lme VIOLATED; reported with caveat; varIdent(~1|sample) noted
  as documented, NOT-adopted alternative.

-- C16 PER-TIMEPOINT WT-CONTRASTS (primary/plotted; "not restored to WT") --
- paired t-test, block = experiment_label, WT SUBSET to each suppressor's SHARED
  blocks (D1); block-mean aggregation collapses the two Orc4 lanes. Estimate =
  WT - group. Decorative paired Wilcoxon. Holm WITHIN each timepoint over
  computable WT-contrasts {4R,Orc1,Orc3,Orc4,Orc5,Orc6}; holm_family_size
  recorded. Per-contrast n_pairs and shared_experiments in every row.
- Per-contrast block n (timepoint-independent here): WT-vs-4R 12; WT-vs-Orc1 2;
  WT-vs-Orc3 3; WT-vs-Orc4 2; WT-vs-Orc5 3; WT-vs-Orc6 3.
- NON-COMPUTABLE handling: n<2 OR sd(paired_diff)==0 -> NA p, reason logged,
  excluded from that timepoint's Holm family. KNOWN case: WT-vs-Orc1 at t=0 (both
  groups zero in both shared experiments -> zero-variance differences). The
  sd()==0 branch fires BEFORE t.test is called (t.test never sees constant data).
- CSV: atpase_per_timepoint_wt_contrasts.csv (consumed by the plot, anti-drift).

-- C17 POOLED MIXED MODELS (PRIMARY for 4R-proximity, D2 -- NOT a diagnostic) --
- Two complementary nlme fits, random = ~1|experiment_block, REML:
  (1) ADDITIVE percent_adp_corrected ~ sample + timepoint, reference = 4R:
      sample coefficients = suppressor-vs-4R POOLED across timepoints; CI from
      intervals(which="fixed"). This IS the reported floor-proximity estimate.
      Rationale for two models (modeling decision): a pooled-across-timepoints
      contrast = the additive sample effect under no-interaction; deriving it as
      a hand-built linear combination of interaction coefficients would be
      "too clever" / error-prone. The additive model gives it directly + CI.
  (2) INTERACTION percent_adp_corrected ~ sample * timepoint; anova() reports
      sample:timepoint (divergence -- WT pulls away) and sample main effect.
      The interaction test is the check that additive pooling is reasonable for
      the near-floor groups.
- timepoint as FACTOR (saturating accumulation, only 4 levels); numeric/linear
  trend documented, NOT adopted.
- Heteroscedasticity caveat: residual variance inflated by WT spread -> suppressor
  -vs-4R CIs CONSERVATIVELY WIDE (do not overstate floor proximity).
- Both fits tryCatch-wrapped: non-convergence -> CONVERGENCE_FAILURE row, script
  stays green (I1); a failed primary fit is a reported finding to investigate.
- CSVs: atpase_pooled_contrasts_vs_4R.csv, atpase_interaction_anova.csv.
- REAL RESULTS [paste real]: pooled suppressor-vs-4R estimate/CI/p per suppressor;
  sample:timepoint F/df/p; sample main F/p; convergence status. (Hand-anchor
  only, R authoritative: WT-vs-4R clearly significant; suppressor-4R expected
  small/near-zero with CI bracketing 0.)

-- C18a / C18b --
- C18a (analyze) writes all stats CSVs above + atpase_negative_values.csv.
- C18b (plot) READS atpase_per_timepoint_wt_contrasts.csv for ALL plotted
  p-values (anti-drift; hard stop if CSV absent). Line/error-bar geometry
  recomputed deterministically from processed_data (display only, not inference).
  Exact Holm-adjusted p + per-contrast n annotated above each timepoint cluster,
  stacked, colored by sample; never stars. n.c. shown for non-computable.

-- NEGATIVE HANDLING (ATPase, encoded) --
- Tripwire bound ATPASE_NEGATIVE_SANITY_BOUND = -0.15 (mirrors WT-t0 0.15
  background threshold). Membership-checked (percent_adp_corrected %in% names)
  BEFORE value-checked (vacuous-guard lesson). NO floor, NO floored column, NO
  pmax. Detect+print+CSV; keep raw into model + differences.

-- DOC STRUCTURE DECISIONS --
- STATISTICAL_METHODS.md: one section per analysis; states block structure,
  test, assumption status, Holm scope, delta-free justification (ATPase),
  firewall + which column feeds which claim (kglut), negative-handling contrast
  (loading floors vs ATPase keeps-raw and WHY). Cites Holm 1979; Pinheiro &
  Bates; Lakens; Snedecor & Cochran.
- STATISTICAL_SELF_STUDY.md: affirming-the-null discussion; TOST considered/
  rejected; firewall cautionary tale (per-kglut WT-pinning -> degenerate
  interaction; global baseline fix; WT-denominator confound of "ORC4R collapses
  at 350"); self-test question bank; primary sources.
- DEPLOYMENT_NOTES.md (C21): renv::snapshot/status; nlme present; clean-room
  restore+run; loading font-gate bare-run gotcha FLAGGED with prerequisite (a)
  and one-line relocation (b) -- relocation NOT performed (behavior change must
  be a deliberate maintainer choice with its own log entry).

-- C21 --
- renv.lock finalized to capture readxl, tidyverse, nlme (stats/nlme ship with R;
  readxl removed Java). renv::status() must report in-sync before deposit.
================================================================================
