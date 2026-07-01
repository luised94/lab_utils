# Self-study reference - the statistics behind the ORC4R suppressor screen

A teaching companion to STATISTICAL_METHODS.md. The methods doc says what the
code did; this one explains why, including the traps we deliberately avoided.

## 1. The designs are blocked - exploit it
WT, the mutant baseline (ORC4R / 4R), and the test genotypes are run TOGETHER in
one replicate (loading) or one experiment (ATPase). The block soaks up gel-to-
gel / day-to-day shifts in absolute signal. A paired/blocked analysis compares
genotypes WITHIN a block, which is far more powerful than an unpaired comparison
that would drown the genotype effect in between-block scatter. In base R this is
`t.test(x, y, paired = TRUE)`; in a model it is a random intercept for the block
(`random = ~1|replicate` or `~1|experiment_label`).

The ATPase blocks are additionally UNBALANCED and INCOMPLETE: the controls ride
along everywhere (so WT/4R accumulate many replicates) but each mutant sits in
only its own ~2-3 experiments. The crucial consequence: a "WT vs Orc5" contrast
can only use the experiments where BOTH appear. WT's 13 replicates do NOT make
the Orc5 contrast stronger - pairing requires co-membership. We subset WT to the
shared blocks and report the per-contrast n honestly (it is 2 or 3, not 13).

## 2. Affirming the null / absence of evidence
The headline biological message for ATPase is partly a claim FOR a null:
"suppressors don't restore hydrolysis." A non-significant test does NOT establish
this - "absence of evidence is not evidence of absence." Worse, at n=2-3 the
power is near zero, so a non-significant result is almost guaranteed regardless
of truth. You cannot turn p>0.05 into "they're the same."

### 2a. The equivalence/TOST alternative (considered, rejected)
The principled way to argue "practically equal" is equivalence testing - TOST
(two one-sided tests; Lakens 2017): pre-specify a margin ë and reject "the
difference is at least ë" from both sides. We rejected it because (i) there is no
biologically defensible ë for this assay, and (ii) at n=2-3 TOST has no power to
reject the equivalence null either, so it would be theatre.

### 2b. What we did instead - the delta-free conjunction
We split the claim into two ORDINARY pieces:
  (a) NOT RESTORED TO WT - a rejectable WT-vs-suppressor test (the load-bearing
      half; if WT is clearly above the suppressor, "not restored" is supported
      by a real rejection, not by a failure to reject);
  (b) NEAR THE 4R FLOOR - the (suppressor - 4R) difference + CI, pooled across
      timepoints. This is descriptive, NOT a p>0.05 equivalence claim. Its
      rigour comes from biology, not from a margin: every suppressor complex
      carries 4R's catalytic mutation, so restoration of hydrolysis is
      mechanistically impossible. The estimate+CI simply quantifies how close to
      4R they sit; a tight CI near zero is consistent with "at the floor."
This is honest: the rejectable half does the lifting; the descriptive half is
labelled descriptive and leans on the shared-mutation argument.

## 3. Pooled vs per-timepoint, and why the pooled model owns the floor claim
"Near 4R ACROSS ALL TIMEPOINTS" is a pooled statement. If we read it off four
separate per-timepoint CIs we invite cherry-picking ("look, the t=45 CI is
tightest"). So the floor-proximity estimate comes from ONE pooled mixed-model
contrast. We use two fits:
  - an ADDITIVE model (sample + timepoint, reference = 4R) whose sample
    coefficients ARE the pooled suppressor-vs-4R differences (+ CI from
    intervals());
  - an INTERACTION model (sample * timepoint) whose anova answers a different
    question: do the curves diverge? WT should pull away over time (significant
    sample:timepoint driven by WT); the suppressor curves should stay flat near
    4R, which is exactly what justifies pooling them additively.
Per-timepoint WT-contrasts remain the plotted, load-bearing "not restored"
result; per-timepoint floor-proximity is secondary by design.

### 3a. timepoint as factor, not a slope
Accumulation of product saturates; it is not linear in time, and with only four
timepoints a linear term would be an unjustified shape assumption. Factor coding
makes no shape assumption and is the conservative default. A numeric/linear (or
polynomial) trend is the documented alternative we did not adopt.

### 3b. The heteroscedasticity you must disclose
WT cells have large variance; near-floor cells have ~0 (some exactly 0). A
homoscedastic lme assumes one residual variance for all - clearly false here.
We report the homoscedastic model with the caveat and note varIdent() as the
heteroscedastic alternative. Direction of the bias matters: pooling WT's big
variance into the residual makes the suppressor-vs-4R SEs (and CIs) WIDER, i.e.
conservative for a "near floor" claim. We never claim tighter than the data
support.

## 4. The normalization firewall - a cautionary tale (kglut)
Normalization can silently destroy the very effect you want to test.

The per-kglut loading column pins WT = 100 at EVERY salt. So WT has zero variance
across salts by construction. Ask "does the WT-ORC4R gap change with salt?" on
that column and the WT side literally cannot move - the genotype x salt
interaction is MECHANICALLY DEGENERATE. You would be testing an artifact of the
normalization, not biology.

The fix is a GLOBAL baseline: normalize everything to WT at 250 mM, so WT is free
to decline as salt rises (and it does: WT replicate means ~100 / 76 / 43 across
250/300/350). Now the interaction is well-posed. The result: once WT's own salt
decline is accounted for, ORC4R's loading deficit does NOT meaningfully worsen
with salt - the gap-change term is small with a wide CI crossing zero,
non-significant. The dramatic per-kglut reading "ORC4R collapses at 350 mM" was
largely the WT DENOMINATOR collapsing. A denominator approaching zero inflates
every ratio above it; that is the same near-zero-denominator pathology that also
broke the loading fold columns.

We hard-wired a guard that is VALUE-checked, not just name-checked: it requires
the WT-across-salt standard deviation to exceed 1, so it fires if the degenerate
column is ever wired into the interaction. (Pointing it at the per-kglut column
during development made it fire on purpose - proof the guard checks values, not a
column name.) Lesson for re-use: a guard that only checks `"col" %in% names(df)`
passes vacuously on a missing or wrong column; assert on the VALUES.

## 5. Negatives: floor (loading) vs keep-raw (ATPase) - same lab, opposite calls
Loading FLOORS small negatives: they are background-subtraction artifacts AND
they poison ratio/fold metrics through near-zero denominators; flooring up toward
zero is safe there and the one offending row is a filtered-out n=1 suppressor.

ATPase does NOT floor. It is analysed by subtraction only (mixed model +
per-timepoint differences); there is no division by a near-zero quantity. Small
negatives are expected noise around a true zero in the near-floor group, and that
noise is the variance the model needs. Flooring would erase it AND bias the
near-zero group's mean UPWARD - which is not one-directional safe here (it would
nudge suppressors away from the floor, weakening exactly the claim under test).
So ATPase DETECTS + REPORTS every negative, keeps them raw, and uses a sanity
TRIPWIRE (stop only on an implausibly large negative, < -0.15) to catch a sign/
load error without touching honest scatter. Same data philosophy (transparency),
opposite handling, because the downstream math differs.

## 6. n=3 honesty
Shapiro-Wilk at n=3 (or n=2) has almost no power: a "pass" is not evidence of
normality. We decline to lean on it and justify parametric tests by the
paired/blocked design plus the assay's track record. Diagnostics are emitted for
transparency and labelled DECORATIVE. The paired Wilcoxon signed-rank companion
floors at p÷0.25 two-sided at n=3 - it literally cannot reach significance, so it
is shown for completeness and never used to decide anything.

## 7. Multiple comparisons - Holm and family scope
Holm (1979) is a step-down procedure: uniformly more powerful than Bonferroni
under the same (no) assumptions, controlling the family-wise error rate. The hard
part is defining the FAMILY. We fixed it narrowly: WT-contrasts within ONE salt,
or WT-contrasts within ONE timepoint. We never pool across salts or timepoints,
because "is the gap different at 250?" and ".at 350?" are separate questions, not
one family - pooling them would over-correct and is conceptually wrong.

## 8. Why nlme and not lme4 or hand-rolled
nlme ships with R (bare-install runnable), is citable (Pinheiro & Bates 2000),
and is dependency-light. lme4 adds transitive dependencies; hand-rolling REML on
an unbalanced/incomplete design is indefensible. For the reported paired results
we used base `stats::t.test`; nlme is only for the mixed models.

## 9. Self-test question bank
1. Why does subsetting WT to shared blocks LOWER n for a mutant contrast, and why
   is that correct rather than wasteful?
2. WT-vs-Orc1 at t=0 is "not computable." What exactly fails, and why is dropping
   it from the Holm family (rather than counting it) the right call?
3. You're handed a non-significant WT-vs-suppressor ATPase result at n=2. Does it
   support "not restored"? Does it support "restored"? Explain both.
4. Write the TOST you would run if a defensible ë existed, and state why n=2-3
   still defeats it.
5. The per-kglut interaction returns a clean-looking p. Why is it meaningless?
   What single quantity, computed on WT only, proves the degeneracy?
6. The suppressor-vs-4R CI is wide. Given the heteroscedasticity, is that
   evidence against "near floor," or an expected conservative artifact? Which?
7. Why does ATPase keep raw negatives while loading floors them? Name the one
   downstream operation that flips the safety argument.
8. In the interaction model, which term encodes "WT pulls away over time," and
   why would you EXPECT it to be significant while the suppressor curves stay
   flat?
9. Why timepoint-as-factor rather than a linear slope? When would numeric be
   defensible?
10. A reviewer asks for "the WT vs ORC4R salt-sensitivity comparison." Which two
    framings answer it, which normalization feeds which, and what's the punchline
    once WT's own decline is removed?

## 10. Primary sources to read
- Holm (1979) - the correction and the family concept.
- Pinheiro & Bates (2000), *Mixed-Effects Models in S and S-PLUS* - random
  effects, REML, intervals(), unbalanced data.
- Lakens (2017) - equivalence testing / TOST and choosing a margin (the
  alternative we rejected).
- Snedecor & Cochran (1989), *Statistical Methods* - paired comparisons and
  (in)complete block designs.
- Wasserstein & Lazar (2016), "The ASA statement on p-values" - context for the
  affirming-the-null discussion.
