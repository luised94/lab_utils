You are implementing Thread 2b of a multi-thread plan to add statistical
analysis to R scripts for an ORC4R suppressor screen (PLOS ONE publication,
Zenodo deposit). Threads 1 (infrastructure) and 2a (350mm loading statistics)
are DONE. You do the kglut TITRATION loading assay. The full plan and global
invariants are in the attached IMPLEMENTATION_HANDOFF.md -- read it first, in
full. Then read the attached DECISION_LOG (Threads 0+1+2a): it carries 2a's
decisions AND a VERIFIED FACTS block whose warnings are load-bearing for you.

Attached:
- IMPLEMENTATION_HANDOFF.md                      (plan -- authoritative)
- DECISION_LOG.md                                (Threads 0+1+2a -- read all of it)
- orc4r-screen_loading-kglut-titration.R         (the ONLY script you edit)
- the kglut data/output CSVs                     (full data -- you need these
  to count negatives and confirm shape; see "VERIFY" below)
- the 2a 350mm output CSVs                        (context: shows the pattern
  your C11/C12 mirror, and the CI-width reality you should expect)

YOUR SCOPE: commits C11, C12, C13, C14a, C14b -- the kglut titration, end to
end. Do NOT touch the 350mm script (Thread 2a, done) or ATPase (Thread 3).
Edit only the kglut script.

STYLE CONTRACT (non-negotiable): procedural, no helpers, no abstractions,
explicit descriptive names, "nothing too clever." Dependencies: base R stats
+ nlme only, nothing hand-rolled. Plots show EXACT p-values, never stars.

GLOBAL INVARIANTS (full text in handoff; the ones that bite here):
- I1 "green" = script source()s clean AND its assertions pass on real data.
- I2 Holm family = WT-contrasts within a SINGLE salt. NEVER pooled across the
     three salts (250/300/350) -- each salt is its own family. State it.
- I3 NORMALIZATION FIREWALL: the salt-sensitivity interaction (C13) MUST use
     percent_wildtype_global (WT-at-250mM baseline), NOT the per-kglut column
     where WT is pinned to 100 at every salt (which makes the interaction
     degenerate). C13 carries a runtime assertion that FIRES if the per-kglut
     column is fed to the model. This assertion must be membership-checked AND
     value-checked (see carry-forward #1) -- it must confirm it is reading the
     GLOBAL column, not merely that some column exists.
- I4 n=3 HONESTY: reuse 2a's exact caveat wording (in the DECISION_LOG) --
     normality untestable at n=3; justification = paired/blocked design + assay
     track record, not a passing normality test. Wilcoxon companions DECORATIVE
     (floor p=0.25), never leaned on.
- I5 Mixed models that stand in for a reported paired test are PRE-COMMITTED
     DIAGNOSTICS. NOTE: C13 is DIFFERENT -- the salt-sensitivity interaction is
     itself a REPORTED result (it answers the reviewer's salt-sensitivity
     question), not a diagnostic. It is the per-salt paired tests (C12) that are
     the simple framing and the interaction (C13) that is the structured one;
     BOTH are reported (handoff: comprehensive by design). Do not mislabel C13
     as diagnostic-only.

CARRY-FORWARD (from Threads 1 + 2a -- these change how you write code):
1. THE VACUOUS-GUARD LESSON (Thread 1). df[["missing_col"]] returns NULL and
   any(is.infinite(NULL)) is FALSE, so a stopifnot can pass without checking
   anything. EVERY guard that names a column MUST first assert existence:
       stopifnot("the_column" %in% names(the_data))
   before testing values. This applies with FULL FORCE to the I3 firewall: it
   must confirm it is bound to percent_wildtype_global specifically.
2. 2a's diagnostic structure (C6) that your C11 MIRRORS: per-label table
   {n, mean, sd, variance} + a one-row pooled-residual Shapiro {W, p, n,
   caveat}. Bartlett/Levene OMITTED (WT zero-variance makes them undefined).
   Your C11 mirrors this BUT ADDS a per-salt split: per-salt residuals AND
   pooled residuals. Reuse 2a's CSV column scheme and caveat wording verbatim
   where they apply, so the two loading analyses read as one coherent method.
3. THE CI-WIDTH REALITY (2a's real output). At 350mm, every (suppressor-ORC4R)
   effect-size CI crossed zero (n=3, high spread). EXPECT the same here: the
   per-salt WT-vs-ORC4R differences may be significant (WT is far from ORC4R)
   but the salt x genotype INTERACTION (C13) may well be NON-SIGNIFICANT with a
   wide CI. That is a VALID, publishable answer to the reviewer ("we quantified
   the salt sensitivity; here is the estimate and its uncertainty"). Do NOT try
   to force significance, do NOT switch tests to chase a smaller p. Report the
   estimate + CI + p honestly whatever it is.

*** NEGATIVE-VALUE HANDLING (user decision -- encode exactly this) ***
The kglut data contains at least one negative percent_wildtype (2a flagged
+5sofa, rep 3, 350mM = -17.58%, from a negative net_intensity). The user has
determined negatives are QUANTIFICATION ARTIFACTS (background-subtraction in
low-signal lanes), and the chosen handling is FLOOR-AFTER-DETECTION, NOT silent
floor and NOT exclusion (exclusion would break the balanced n=3 pairing
asymmetrically; flooring preserves the design and can only move a value UP
toward zero, so it cannot manufacture a restorer effect or exaggerate a
WT-mutant gap -- that is the argument that makes it safe; put it in the log).
Implement in this order, as its own commit step EARLY (before C12 tests run):
  (a) DETECT: before normalization/tests, scan for net_intensity < 0 and the
      resulting percent_wildtype < 0. PRINT a table of EVERY such row
      (sample, replicate, salt, raw net_intensity, raw percent value). No
      silent handling -- the user must see the full list to confirm they are
      all small-magnitude artifacts and not a systematic failure.
  (b) ENVELOPE ASSERT: stopifnot that negatives are FEW and SMALL relative to
      WT signal (choose an explicit, documented threshold, e.g. count below a
      stated bound AND |value| below a stated bound). A systematic problem must
      STOP the script loudly, not get floored away. Membership-checked.
  (c) FLOOR into a NEW explicitly-named column:
      percent_wildtype_floored <- pmax(percent_wildtype, 0). Do NOT overwrite
      the raw column -- it stays for transparency and the methods doc.
  (d) Run the per-salt tests (C12) and -- decide and document -- whether the
      global-baseline column for C13 is also computed from floored values.
      (Recommend: yes, derive percent_wildtype_global from the floored values
      so a single artifact does not distort the salt-sensitivity baseline; but
      surface this choice explicitly, it is a real decision.)
  (e) LOG every floored row in the DECISION_LOG and in an output CSV; the
      methods-doc note is drafted in the log for Thread 3 to consolidate.
This detection step ALSO answers the open question "is it only one negative?"
at runtime -- report the count.

NARRATIVE FIREWALL (handoff F6): kglut is WT/ORC4R/+4sofa only. The
(+4sofa - ORC4R) comparison, if computed, is an ORDINARY effect size -- NOT the
ATPase floor-proximity claim. Same firewall comment 2a used; reuse it.

VERIFY BEFORE YOU EDIT (do not assume):
- Open the kglut data CSV and confirm ACTUAL shape: genotypes present
  (expect WT/ORC4R/+4sofa ONLY -- confirm the n=1 rotating-suppressor rows are
  filtered out and never enter a test, per handoff F4); the three salts
  (250/300/350); replicate counts per genotype x salt cell (expect n=3);
  the blocking variable name (2a found `replicate`, values 1/2/3 -- confirm
  here).
- COUNT the negatives (detection step above) and report the full list before
  proceeding.
- Confirm BOTH normalization columns exist and their real names: the per-kglut
  (WT=100 within each salt) and the global-baseline (WT-at-250mM). Before
  fitting C13, demonstrate the degeneracy: show WT is FLAT across salts in the
  per-kglut column and NOT flat in the global column. That observation is the
  justification for the firewall -- capture it for the log.
- Confirm pairing completeness within each salt's blocks (WT/ORC4R/+4sofa
  co-occur in every replicate at every salt). Report any incomplete cell;
  it changes which pairs are usable and whether the balanced interaction holds.

WORK ONE COMMIT AT A TIME. For each:
- state what you compute and which CLAIM/framing it supports (C12 = per-salt
  simple framing, direct reviewer WT-vs-ORC4R answer; C13 = structured
  salt-sensitivity interaction, the reviewer's salt-sensitivity question
  quantified),
- state the Holm family (I2: within-salt) and n=3 caveat (I4) where they apply,
- make the edit, show the changed block,
- state how you confirmed it (what you read in the data; the membership +
  value-checked guard).
Keep it tight: a few sentences per commit.

SANCTIONED DELIBERATE RED -- C13 ONLY (the only one in the whole project).
When you implement the I3 firewall: FIRST point the assertion at the per-kglut
column, confirm it FIRES (red = the guard works), then repoint to
percent_wildtype_global and confirm green. Show BOTH states explicitly. The red
is the evidence the guard is real. Everywhere else, every commit lands green.

C13 MODELING CHOICES (document each):
- nlme: percent_wildtype_global ~ genotype * salt, random = ~1|replicate.
- salt as FACTOR (3 levels, no monotonic/trend assumption) is the default;
  NOTE the ordered/numeric (trend) alternative in the log but do not adopt it
  unless justified. The genotype:salt interaction term IS the salt-sensitivity
  result; report estimate, CI, p for the WT-vs-ORC4R gap change with salt.
- If the model has convergence trouble (2a's lme converged despite WT-pinning,
  but the global column un-pins WT so residual structure differs), wrap in
  tryCatch and report -- a failure is a real finding, not a code bug.

DECISION_LOG: append under "[Thread 2b appends here]" every non-obvious choice
(per-salt test details, within-salt Holm, the negative-detection threshold and
the floored-row list, whether global-baseline is derived from floored values,
salt-as-factor vs numeric, the firewall observation, C13 estimate/CI/p and
whether the interaction was significant). Output the FULL updated log at the
end for Thread 3.

WHAT I CANNOT DO: I (the user) run the script, not you. You cannot execute to
confirm green; reason explicitly about why it sources clean, then give me at
the end: the edited kglut script, the exact command to run it, and the exact
outputs to paste back (console trail including the negative-detection table and
BOTH C13 firewall states; the augmented CSVs; the C13 interaction estimate/CI/p;
confirmation the PDF rendered with exact p-values). Those outputs are an input
to Thread 3.

If anything in the handoff or carry-forward is wrong against the ACTUAL data --
STOP and raise it before editing. Data and code are ground truth. In
particular: if the negative count or magnitude exceeds the envelope (systematic,
not a stray artifact), STOP and report -- do not floor your way past it. And if
the per-kglut vs global-baseline degeneracy is not what we assumed, that
overturns C13 and must be surfaced.

Start with your read-back: the actual kglut data shape (genotypes, salts,
per-cell replicate counts, blocking-variable name), the FULL list of negative
rows with counts, confirmation of both normalization columns and a demonstration
of the per-kglut WT-flatness degeneracy, and your plan for the C13 firewall
assertion. Then proceed commit by commit.
