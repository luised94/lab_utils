================================================================================
DECISION_LOG (append-only) -- THREAD 0 / PLANNING BLOCK
Place this block at the TOP of the consolidated log, above Threads 1/2a/2b.
================================================================================
[Thread 0 / planning -- scope, framing, and binding constraints]

-- ASSAY ROLES (ground truth carried by all threads) --
- Three datasets/analyses: (1) kglut titration loading = WT/ORC4R/+4sofa only,
  n=3, salts 250/300/350; (2) 350 mM all-suppressors loading = WT, ORC4R, and
  all five suppressors (+1/+3/+4/+5/+6sofa), n=3, single salt; (3) ATPase
  timecourse = WT, 4R, all five suppressors, timepoints 0/15/45/90, UNEQUAL
  replicates (WT/4R high, mutants n=3).
- Loading restorers = Orc1, Orc4. Reducers = Orc3, Orc5, Orc6. 4R < WT in
  loading. Grouping is comment-level only; all suppressors get identical
  statistical treatment.

-- ATPase NULL FRAMING: DELTA-FREE (no equivalence margin) --
- The claim "suppressors do not restore ATPase activity" is a claim FOR the
  null, which a non-significant test cannot establish (absence of evidence is
  not evidence of absence; n=3 has near-zero power so non-significance is
  near-guaranteed regardless of truth).
- Equivalence testing (TOST) was CONSIDERED and REJECTED: it requires a
  pre-specified margin delta, and no biologically defensible delta exists for
  this assay; at n=3 TOST also lacks power to reject the equivalence null.
- Chosen framing instead: the claim is the conjunction of two ordinary,
  rejectable/descriptive parts -- (a) "not restored to WT" = WT-vs-suppressor
  SIGNIFICANT (rejectable test, load-bearing); (b) "near the 4R floor" =
  (suppressor - 4R) difference + CI (descriptive, NOT a p>0.05 equivalence
  claim). Justification for (b): all suppressor complexes carry the SAME
  catalytic mutation as 4R, so no restoration of hydrolysis is mechanistically
  possible. This shared-mutation argument is what makes the delta-free framing
  rigorous without an arbitrary margin. Cite Lakens (TOST) as the
  considered-and-rejected alternative in the methods doc.

-- TERMINOLOGY (used in code comments + both docs) --
- "non-restoration-to-WT contrast" = WT vs suppressor; a rejectable test;
  load-bearing half of both the ATPase and loading claims.
- "4R-floor-proximity estimate" = (suppressor - 4R) difference + CI; ATPase
  ONLY; rests on the shared catalytic mutation.
- In LOADING, the (suppressor - 4R/ORC4R) comparison is an ORDINARY effect
  size, NOT the floor-proximity claim -- reducers fall below 4R BY DESIGN. The
  two narratives must never be conflated (binding firewall, see deltas below).
- "salt-sensitivity interaction" = salt x genotype term in the kglut mixed
  model; quantifies whether the WT-ORC4R gap changes with salt.
- "block" = replicate (loading) or experiment (ATPase).

-- DESIGN: LOADING AND ATPase ARE PAIRED/BLOCKED --
- WT, 4R/ORC4R, and the mutant(s) are co-run within one replicate/experiment;
  block = replicate (loading) / experiment (ATPase). Pairing is exploited
  (more powerful than unpaired). ATPase blocks are UNBALANCED and INCOMPLETE
  (WT/4R ride along in nearly every experiment -> higher replicate counts;
  each mutant appears in only its ~3 blocks).

-- STATISTICAL TOOLING: base R stats + nlme ONLY; nothing hand-rolled --
- Everything possible in base `stats`: t.test (paired/Welch), wilcox.test,
  aov/lm, p.adjust(method="holm"), shapiro.test, confint, residuals/qqnorm.
- Mixed models via `nlme` (ships with R, citable Pinheiro & Bates,
  dependency-light, reviewer-friendly). lme4 rejected (transitive deps);
  hand-rolling REML rejected (indefensible on unbalanced ATPase). renv locks
  nlme. Both base stats and nlme ship with R -> bare-install runnable.

-- MULTIPLE COMPARISONS: Holm, pre-specified families --
- Holm (uniformly more powerful than Bonferroni, same assumptions) over
  PRE-SPECIFIED contrast families, NOT all pairwise.
- FAMILY SCOPE (binding, I2): family = {WT-contrasts within a SINGLE salt} or
  {WT-contrasts within a SINGLE timepoint}. NEVER pooled across salts or across
  timepoints (those are separate questions, not one family).

-- n=3 HONESTY (binding, I4) --
- At n=3, normality cannot be meaningfully tested (Shapiro at n=3 has almost no
  power; a passing test is NOT a license to assume normality). Justification
  for parametric tests rests on the paired/blocked design + the assay's
  established track record, NOT on a normality test. Diagnostics are still
  emitted for transparency, explicitly labeled. Wilcoxon companions are
  DECORATIVE (n=3 paired floors two-sided p ~ 0.25), included for transparency,
  never leaned on.

-- ANNOTATION: exact p-values on plots, never stars (binding) --

-- PER-ANALYSIS PRIMARY FRAMINGS (set in planning, refined by review) --
- 350 mM: WT-contrasts (paired t + Wilcoxon companion + Holm within family) as
  the non-restoration claim; (suppressor-ORC4R) ordinary effect size + CI;
  mixed-model cross-check as a PRE-COMMITTED DIAGNOSTIC (I5).
- kglut: BOTH framings, comprehensive by design -- (a) per-salt independent
  paired WT-vs-ORC4R tests (simple, direct reviewer answer per salt), and
  (b) salt x genotype mixed-model interaction (the salt-sensitivity result).
  Reviewer point being answered: "Salt sensitivity: provide quantitative WT vs
  ORC4R comparison."
- ATPase: per-timepoint analysis PRIMARY/plotted for the WT-difference
  (restoration) half; pooled mixed model PRIMARY for the 4R-proximity half
  ("near 4R across all timepoints" is a pooled claim, so it leans on the pooled
  contrast, not cherry-picked per-timepoint CIs).

-- ADVERSARIAL-REVIEW DELTAS (binding constraints folded in pre-implementation) --
- D1 (Lens 1): ATPase per-timepoint WT-contrasts must subset WT to each
  suppressor's SHARED blocks (incomplete design); report per-contrast n. WT's
  full replicate count does NOT strengthen these contrasts.
- D2 (Lens 4): the pooled mixed model is the PRIMARY home for the "near 4R
  across all timepoints" claim; per-timepoint 4R-proximity is secondary, to
  avoid cherry-pick exposure.
- D3 (Lens 3 / I5): the 350 mM mixed model is a DECLARED consistency
  diagnostic; the paired t-test is the pre-committed reported result;
  disagreement triggers investigation, never selection of the friendlier p.
- D4 (Lens 2 / I3 NORMALIZATION FIREWALL): the kglut salt-sensitivity
  interaction MUST use the GLOBAL-baseline normalization (WT@250mM), NOT the
  per-kglut column where WT is pinned to 100 at every salt (which makes the
  interaction mechanically degenerate). C13 carries a runtime assertion that
  FIRES if the per-kglut column is fed to the model -- and (post Thread 1) is
  value-checked, not merely membership-checked, so it cannot pass vacuously.
- D5 (Lens 5 / narrative firewall): the loading (X-4R) effect size is NOT the
  ATPase floor-proximity claim; reducers below 4R are expected. A mandatory
  in-code comment enforces the distinction.
- D6 (Lens 6 / I2): Holm family scope fixed globally (see above).
- D7 (STEP-3 blast radius): nlme restore extracted into its own early commit
  (C0d) because C9/C13/C17 consume nlme before the final lockfile commit; the
  lockfile commit was serving two meanings (make-available vs finalize-deposit)
  and was split.

-- PATHS + INVOCATION --
- source()-only path support; Rscript/interactive deliberately unsupported
  (they do not set 'ofile', so the script directory cannot be resolved).
  Resolution chain: script-relative (Zenodo co-located) -> MC_DROPBOX_PATH
  (working location) -> stop() naming BOTH attempted absolute paths. Both the
  original data location and the Zenodo copy remain live; xlsx->readxl migration
  removes the Java/rJava dependency for bare-install runnability.

-- DOCS (generated LAST, by consolidation not reconstruction) --
- STATISTICAL_METHODS.md (Zenodo-facing, terse/defensive) + a separate
  self-study reference (expansive: concepts, the affirming-the-null discussion,
  TOST as considered-and-rejected, primary sources [Holm 1979; Pinheiro & Bates;
  blocked-design refs; Lakens on equivalence], a self-test question bank).
  Both built in the FINAL commit from the accumulated DECISION_LOG so they
  describe what the code ACTUALLY did. Methods doc must LABEL which normalization
  column feeds which claim.

-- ATPase NEGATIVE-VALUE HANDLING (decided after kglut, recorded here for the
   consolidated log; differs from loading ON PURPOSE) --
- Loading floored negatives (background-subtraction artifacts that also broke
  fold/ratio metrics via near-zero denominators). ATPase does NOT floor:
  it is analyzed by mixed model + per-timepoint DIFFERENCES (subtraction-based,
  negative-safe; no division), small negatives are EXPECTED noise around a true
  zero in the near-zero 4R/suppressor group, and that noise is data the model
  must absorb. Rule: DETECT + REPORT every negative (sample, replicate,
  timepoint, raw value), keep RAW (no floor, no new column), and use the
  envelope assertion as a SANITY TRIPWIRE only (stop on an implausibly large
  negative that signals an upstream sign/load error, not on small near-zero
  scatter). Thread 3 must also verify no division-based descriptive column
  (e.g. fold-over-4R, percent-of-WT) is fed by negatives; any such column is
  descriptive-only and kept out of all claims (same near-zero-denominator
  pathology seen in the 350 mM summary fold columns).
================================================================================
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

[Thread 2a]
-- TEST PER CONTRAST (C7) --
- WT-contrasts (non-restoration-to-WT claim): paired t-test, block = replicate,
  WT vs each of {ORC4R,+1,+3,+4,+5,+6sofa}. Implemented via
  t.test(WT_vec, group_vec, paired=TRUE) on percent_wildtype, both arrange()d
  by replicate so pairs align by block. estimate reported as WT - group.
- Decorative paired Wilcoxon companion per contrast (I4): two-sided signed-rank
  floors at p=0.25 at n=3; suppressWarnings around it; never leaned on.

-- HOLM FAMILY SCOPE (I2) --
- Family = the SIX WT-contrasts within the SINGLE 350 mM condition. p.adjust(
  method="holm") applied across exactly these six. No salts/timepoints exist in
  this script, so no cross-family pooling is possible or permitted.

-- C8 EFFECT-SIZE FRAMING (narrative firewall / F6) --
- (suppressor - ORC4R) reported as ORDINARY paired effect size + 95% CI via
  t.test(group, ORC4R, paired=TRUE)$estimate/$conf.int. Mandatory in-code
  comment states this is NOT the ATPase floor-proximity argument and that
  reducers (+3/+5/+6sofa) fall below ORC4R BY DESIGN. CI column lineage:
  percent_wildtype (per-replicate WT=100 normalization).

-- C9 MIXED-MODEL DIAGNOSTIC (I5) --
- nlme::lme(percent_wildtype ~ label, random = ~1|replicate). PRE-COMMITTED
  consistency check; reported result remains the C7 paired t-test. WT is factor
  reference, so label coefficients = (group - WT). Fit wrapped in tryCatch:
  non-convergence -> flagged diagnostic row (converged=FALSE), NOT a crash
  (a failed diagnostic is itself an investigation flag; keeps I1 green).
  Documented that WT zero-variance violates lme's homoscedastic-residual
  assumption -> exactly why this is a diagnostic, not the reported test.

-- C6 DIAGNOSTIC STRUCTURE (C11/Thread 2b mirrors this) --
- Two artifacts: (1) per-label table {n, mean, sd, variance};
  (2) one-row pooled within-label-residual Shapiro {W, p, n, caveat}.
- QQ emitted as the numeric Shapiro statistic, NOT a font-gated plot, so
  diagnostics survive a font-less bare run. Bartlett/Levene deliberately
  OMITTED: WT zero-variance makes them undefined.
- Thread 2b (C11) must MIRROR this but ADD a per-salt split (per-salt residuals
  + pooled residuals), per handoff.

-- n=3 CAVEAT WORDING (I4, reused verbatim in 2b) --
  "with 3 biological replicates per condition, normality CANNOT be meaningfully
   tested (Shapiro-Wilk at n=3 has almost no power and is NOT a license to
   assume normality). Justification rests on the PAIRED/BLOCKED design plus the
   assay's established track record -- NOT on a normality test passing."

-- CSV COLUMN SCHEME (C10a) --
- wt-contrasts CSV: comparison, group_label, claim, test, holm_family, n_pairs,
  estimate_wt_minus_group, ci_lower/upper_wt_minus_group, t_statistic, df,
  raw_p_value, holm_adjusted_p_value, wilcoxon_p_value_decorative.
- effect-size CSV: suppressor, quantity(=firewall note), n_pairs,
  estimate_group_minus_orc4r, ci_lower, ci_upper.
- diagnostics-per-label CSV + diagnostics-normality CSV + mixed-model CSV.
- summary CSV augmented: per-label wt_contrast_* and effect_* columns
  (WT row + ORC4R row carry NA where a quantity doesn't apply).
- full-data CSV augmented: per-row group-level wt_contrast raw/holm p.
- Joins on explicit as.character keys to avoid factor/character join errors.

-- PLOT (C10b) --
- Exact Holm-adjusted p above each non-WT bar (geom_text), never stars.
  Caption names the Holm family + n. Annotation df from C7 results (upstream),
  display-only, font-gate-independent for compute.

-- FLAGGED, OUT OF SCOPE (for C21/bare-run owner) --
- The Arial/Cairo gate sits at the TOP of the loading script, above data load.
  So a font-less environment stops before ANY compute, including stats. I kept
  all stats compute + stats-CSV writes upstream of the Plot block (survives a
  ggsave failure) but could NOT put compute upstream of the font gate without
  relocating it (Thread-1 infrastructure, not in my commit list). If a future
  bare run must produce stats without fonts, relocate the font gate to just
  above the Plot section.

================================================================================
VERIFIED FACTS -- 350 mM all-suppressors (Thread 2a, confirmed against real CSVs)
================================================================================
- Blocking variable name in code/data: `replicate` (NOT `gel`); values 1,2,3.
- 7 labels present: WT, ORC4R, +1sofa, +3sofa, +4sofa, +5sofa, +6sofa.
- Replicate count per condition: 3 each. Total 21 rows.
- Pre-normalized values live in `percent_wildtype`; WT = 100 exactly in every
  replicate -> WT within-group variance = 0 (matters for C6 var-tests + C9 lme).
- PAIRING COMPLETE: WT, ORC4R and every suppressor co-occur in EVERY replicate;
  all 6 WT-contrasts and all 5 effect-sizes use full n=3. No incomplete block.
- All percent_wildtype values >= 0 in the 350 mM aggregate file (the existing
  >=0 guard passes here).
  *** WARNING TO THREAD 2B: the kglut TITRATION context CSV
  (loading_processed_full_data.csv) contains a NEGATIVE percent_wildtype
  (+5sofa, replicate 3, 350 mM kglut: net_intensity -2456.8 -> -17.58 %WT).
  That sub-zero pathology is NOT in this 350 mM aggregate file but MAY appear
  in the kglut script -- the >=0 guard and any fold/ratio Inf/NaN guard matter
  there, and a negative value makes percent_difference/fold metrics nonsensical
  for that row. Verify before running paired tests on kglut. ***
- ASSUMPTIONS THAT DID NOT HOLD CLEANLY:
  (1) Variance homogeneity is NOT assessable across groups (WT zero-variance).
      Standard Bartlett/Levene undefined -> omitted by design.
  (2) Normality untestable at n=3 (I4) -- not an assumption met, an assumption
      we explicitly decline to test; justified by paired design + track record.
  Thread 2b inherits BOTH caveats for kglut (also n=3, also WT-pinned).
- EXPECTED p-values (hand-calc; R authoritative): all 6 WT-contrasts raw p in
  [~0.0002, ~0.008], all significant after Holm (Holm p <= ~0.013). Wilcoxon
  companions all = 0.25 (floor). Effect sizes (supp - ORC4R): +4 ~ +17, +1 ~ +15
  (restorers above 4R); +3 ~ +2.7; +5 ~ -2.0, +6 ~ -2.9 (reducers below 4R).
  REPLACE WITH REAL VALUES once the user pastes the run output.
================================================================================
[Thread 2b / kglut titration loading -- C11..C14b]
- DATA SHAPE VERIFIED (kglut full-data CSV): blocking var = `replicate`
  (1/2/3). Salts 250/300/350. Test genotypes = WT/ORC4R/+4sofa ONLY, n=3 per
  genotype x salt, BALANCED and COMPLETE (every replicate has all three at
  every salt). Rotating suppressors +1sofa(rep1)/+3sofa(rep2)/+5sofa(rep3) are
  n=1, filtered out, never enter any test (F4 confirmed).
- NEGATIVE-VALUE CENSUS: exactly ONE negative in the whole dataset --
  +5sofa, replicate 3, 350 mM: net_intensity -2456.769, percent_wildtype
  -17.579 (and percent_wildtype_global -7.561). It is a low-signal background-
  subtraction artifact. CRITICALLY it sits in +5sofa, a rotating n=1 suppressor
  that is filtered out of all tests AND the interaction -> the floor is INERT
  for every reported value; it is applied for transparency/robustness only.
- NEGATIVE HANDLING = floor-after-detection (user decision). Implemented EARLY
  (pre-C12): detect+PRINT every offending row; envelope assert with documented
  thresholds NEGATIVE_VALUE_MAX_COUNT = 3 (of 36 rows) and
  NEGATIVE_VALUE_MAX_ABS_PERCENT = 50 (vs WT = 100) -- a breach STOPS the
  script (systematic failure must not be floored away); floor into NEW columns
  percent_wildtype_floored / percent_wildtype_global_floored (raw retained).
  Safety argument: floor can only move a value UP toward 0, so it cannot
  manufacture a restorer effect or widen a WT-mutant gap. Floored rows logged
  to loading_kglut-titration_floored-rows.csv.
- GLOBAL-BASELINE FROM FLOORED VALUES (carry-forward 1d decision): YES, C13 uses
  percent_wildtype_global_floored so a single artifact cannot distort the salt-
  sensitivity baseline. For the test genotypes this equals the unfloored values
  (none negative), so the choice is consequential only in principle here.
- FIREWALL DEGENERACY CONFIRMED (I3): per-kglut column pins WT = 100 at every
  salt (within-WT sd across salts = 0). Global column: WT declines with salt
  (sd of 9 WT values ~ 26: rep means 100/76.4/43.1). So genotype:salt
  interaction is mechanically degenerate on per-kglut, well-posed on global.
- C11 DIAGNOSTICS: mirror 2a C6 (per-cell {n,mean,sd,variance} + pooled-residual
  Shapiro) PLUS a per-salt residual split. Bartlett/Levene OMITTED (WT zero-
  variance undefined). Residual = value - (label x salt) cell mean. I4 caveat:
  normality untestable at n=3; Shapiro DECORATIVE; justification = paired/blocked
  design + assay track record. CSVs: diagnostics-cellstats, diagnostics-shapiro.
- C12 PER-SALT TESTS (SIMPLE framing, direct reviewer WT-vs-ORC4R answer):
  paired t-test (block = replicate), on percent_wildtype_floored, three
  contrasts per salt {WT-vs-ORC4R, WT-vs-+4sofa, ORC4R-vs-+4sofa}. Wilcoxon
  signed-rank companions DECORATIVE (n=3 paired floor p ~ 0.25). Pairing guard
  asserts equal length, n=3, AND replicate-order identity before each test.
- HOLM FAMILY (I2): WT-contrasts {WT-vs-ORC4R, WT-vs-+4sofa} WITHIN each salt
  (m=2), NEVER pooled across salts. ORC4R-vs-+4sofa is NOT a WT contrast ->
  reported raw, OUTSIDE the family (it is the narrative-firewall ordinary
  effect size, F6: a loading deviation by design, NOT the ATPase floor claim).
- EXPECTED C12 (to confirm on run): WT-vs-ORC4R raw/Holm p ~ 250: 0.027/0.053,
  300: 0.059/0.119, 350: 0.0016/0.0032; WT-vs-+4sofa raw ~ 250: 0.21, 300:
  0.20, 350: 0.0024; ORC4R-vs-+4sofa raw ~ 250: 0.003, 300: 0.010, 350: 0.036.
  [REPLACE with real run output.]
- C13 SALT-SENSITIVITY INTERACTION (STRUCTURED framing, REPORTED result, NOT a
  diagnostic per I5 note): nlme::lme(percent_wildtype_global_floored ~
  label * kglut, random = ~1|replicate, REML). salt as FACTOR (default; no
  trend assumption). Trend/numeric (linear over 250/300/350) alternative
  documented and NOT adopted (only 3 levels; decline need not be linear;
  factor is the conservative choice). WT-vs-ORC4R gap-change terms
  labelORC4R:kglut300 and labelORC4R:kglut350 ARE the salt-sensitivity result
  (estimate/CI/p). lme + intervals() wrapped in tryCatch -> non-convergence is
  reported as a finding, not a code bug.
- C13 FIREWALL ASSERTION: membership-checked first (NULL-column lesson), then
  load-bearing VALUE check sd(WT values across salts) > 1 (per-kglut -> 0 ->
  fires), then NAME check (must be a global column). SANCTIONED DELIBERATE RED:
  pointing interaction_response_column at "percent_wildtype_floored" fires the
  VALUE assertion (proving the guard checks values, not just a name); reverting
  to "percent_wildtype_global_floored" passes. [PASTE both states from run.]
- EXPECTED C13 (to confirm on run, carry-forward #3 reality): on the GLOBAL
  baseline WT and ORC4R decline together, so the WT-ORC4R gap change 250->350
  is SMALL (labelORC4R:kglut350 ~ +3 %WT-units, gap slightly NARROWS) with a
  WIDE CI crossing zero and a NON-significant p. This is a valid, publishable
  answer: once WT's own salt decline is accounted for, ORC4R's loading deficit
  does not meaningfully worsen with salt; the dramatic per-kglut "ORC4R
  collapses at 350" is largely WT itself declining. Do NOT chase significance.
  [REPLACE estimate/CI/p with real run output.]
- C14a CSVs: loading_kglut-titration_per-salt-contrasts.csv (C12) and
  loading_kglut-titration_salt-interaction.csv (C13), each carrying a
  normalization_column lineage field (per-kglut feeds C12; global feeds C13).
- C14b PLOT: faceted_by_kglut_plot annotated with per-salt WT-vs-ORC4R EXACT p
  (raw + Holm-within-salt) above each panel; salt-sensitivity interaction
  estimate/CI/p bound into the plot caption. EXACT p-values, never stars.
- DISCREPANCY FLAGGED (non-blocking): DECISION_LOG.md and the 2a 350mm stat
  CSVs were not in the Thread 2b attachment set; CSV column schemes here MIRROR
  the handoff-described 2a scheme rather than being copied from a seen file.
===================================================================
