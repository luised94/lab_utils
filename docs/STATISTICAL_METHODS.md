# Statistical methods - ORC4R suppressor screen

This document describes the statistics actually performed by the deposited R
scripts. Tools: base R `stats` for everything except mixed models; `nlme`
(Pinheiro & Bates 2000; ships with R) for mixed models. Nothing is hand-rolled.
All plots annotate EXACT p-values, never significance stars.

Three analyses: (1) kglut titration loading; (2) 350 mM all-suppressors
loading; (3) ATPase timecourse. All three are PAIRED/BLOCKED designs.

Multiple comparisons throughout use Holm's step-down correction (Holm 1979) over
PRE-SPECIFIED families. Family scope is fixed and narrow: WT-contrasts within a
SINGLE salt, or WT-contrasts within a SINGLE timepoint. Families are NEVER
pooled across salts or across timepoints - those are distinct questions, not one
family.

n=3 honesty (loading) / n=2-3 honesty (ATPase mutants): at these per-cell
replicate counts normality cannot be meaningfully tested (Shapiro-Wilk has
almost no power and a pass is not a license). Parametric tests are justified by
the paired/blocked design and the assays' established track record, not by a
normality test. Shapiro-Wilk statistics and Wilcoxon signed-rank companions are
emitted for transparency and are explicitly DECORATIVE (an n=3 two-sided paired
signed-rank test floors at p÷0.25). They are never load-bearing.

--------------------------------------------------------------------------------
## 1. kglut titration loading

Design: genotypes WT / ORC4R / +4sofa only, n=3 biological replicates, salts
250 / 300 / 350 mM, fully crossed and balanced. Block = replicate. n=1
rotating-suppressor rows are filtered out and never enter any test.

Negative handling: one low-signal background-subtraction artifact (+5sofa,
replicate 3, 350 mM) was floored to 0 after detection and logging; flooring can
only move a value up toward 0 and the offending row is a filtered-out n=1
suppressor, so the floor is inert for every reported value. Loading floors
because its ratio/fold metrics are broken by near-zero denominators (see 4).

Per-salt contrasts (direct reviewer answer): paired t-tests (block = replicate),
WT-vs-ORC4R, WT-vs-+4sofa, and ORC4R-vs-+4sofa, on the per-kglut normalized
column (WT pinned to 100 within each salt). Holm WITHIN each salt over the two
WT-contrasts; ORC4R-vs-+4sofa is reported raw, outside the family (it is an
ordinary effect size, not a WT-restoration contrast).

Salt-sensitivity interaction (the salt x genotype question): `nlme::lme(
percent_wildtype_global ~ genotype * salt, random = ~1|replicate)`, salt as a
FACTOR (no linear-trend assumption; only three levels; decline need not be
linear - the numeric/trend alternative was considered and not adopted).

NORMALIZATION FIREWALL (which column feeds which claim): the per-kglut column
pins WT=100 at every salt, so within-WT variance across salts is zero and a
genotype x salt interaction computed on it is MECHANICALLY DEGENERATE. The
interaction therefore uses percent_wildtype_global (WT-at-250mM baseline) only.
A runtime assertion in the code is value-checked (it requires the WT-across-salt
standard deviation to exceed 1), so it fires if the per-kglut column is ever
fed to the interaction - it cannot pass vacuously. RESULT: on the global
baseline WT and ORC4R decline together with salt; the WT-ORC4R gap change from
250350 mM is small with a wide CI crossing zero and is NON-SIGNIFICANT. The
naive per-kglut reading ("ORC4R collapses at 350 mM") is largely WT itself
declining; the firewall prevents a denominator-confounded claim.

--------------------------------------------------------------------------------
## 2. 350 mM all-suppressors loading

Design: WT, ORC4R, and five suppressors (+1/+3/+4/+5/+6sofa), n=3, single salt,
pairing COMPLETE (all genotypes co-occur in every replicate). Block = replicate.

WT-contrasts (the non-restoration-to-WT claim): paired t-tests (block =
replicate) of WT vs each of the six non-WT genotypes, on percent_wildtype
(per-replicate WT=100). Holm over the six WT-contrasts within the single 350 mM
condition. Decorative paired Wilcoxon companions reported alongside.

Effect sizes: (suppressor - ORC4R) paired difference + 95% CI. This is an
ORDINARY effect size. By design the reducer suppressors (+3/+5/+6sofa) fall
BELOW ORC4R in loading; this is NOT the ATPase 4R-floor-proximity argument and
the two narratives are kept strictly separate.

Mixed-model cross-check: `nlme::lme(percent_wildtype ~ genotype,
random = ~1|replicate)` is a PRE-COMMITTED consistency diagnostic; the reported
result is the paired t-test. Disagreement triggers investigation, never
selection of the friendlier p. (WT's zero within-group variance violates the
lme homoscedastic-residual assumption - exactly why the lme is a diagnostic, not
the reported test here.)

Assumption status: variance homogeneity is not assessable (WT zero-variance 
Bartlett/Levene undefined, omitted); normality untestable at n=3 (declined, see
header).

--------------------------------------------------------------------------------
## 3. ATPase timecourse

Design: WT, 4R, and five suppressors, timepoints 0/15/45/90 min. Response =
percent_adp_corrected = ADP/(ADP+ATP), then minus the No_ORC t=90 background
within each experiment. Block = experiment_label. The block structure is
UNBALANCED and INCOMPLETE: WT (n=13 experiments) and 4R (n=12) ride along in
nearly every experiment, while each mutant appears in only its own experiments
(Orc3/Orc5/Orc6 = 3; Orc1 = 2; Orc4 = 2 experiments / 3 lanes, two lanes in one
experiment). Every sample has all four timepoints; no timepoint is missing.

Delta-free null framing (no equivalence margin): the claim "suppressors do not
restore ATPase activity" is a conjunction of two ordinary parts -
  (a) NOT RESTORED TO WT = WT-vs-suppressor SIGNIFICANT (rejectable; the
      load-bearing half), addressed per timepoint;
  (b) NEAR THE 4R FLOOR ACROSS ALL TIMEPOINTS = (suppressor - 4R) difference +
      CI, POOLED across timepoints by the mixed model (NOT a p>0.05 equivalence
      claim). Justification: all suppressor complexes carry the SAME catalytic
      mutation as 4R, so restoration of hydrolysis is mechanistically
      impossible; this shared-mutation argument is what makes the delta-free
      framing rigorous without an arbitrary margin.
Equivalence testing (TOST; Lakens 2017) was CONSIDERED and REJECTED: no
biologically defensible equivalence margin exists, and at these replicate counts
TOST would lack power to reject the equivalence null anyway.

(a) Per-timepoint WT-contrasts (primary/plotted): paired t-tests (block =
experiment_label). WT is SUBSET to each suppressor's SHARED experiments - a
"WT vs Orc5" test uses only the experiments where both appear; WT's full
replicate count does NOT strengthen these contrasts. The two Orc4 lanes in the
shared experiment are aggregated to one block mean so the experiment is the
paired unit. Per-contrast n is reported in every row (WT-vs-4R n=12;
WT-vs-Orc1 n=2; WT-vs-Orc4 n=2 blocks; WT-vs-Orc3/Orc5/Orc6 n=3). Holm WITHIN
each timepoint over the computable WT-contrasts {4R, Orc1, Orc3, Orc4, Orc5,
Orc6}; the family size used is recorded. WT-vs-Orc1 at t=0 is NOT computable
(both groups read zero in both shared experiments  zero-variance paired
differences) and is reported as such and excluded from that timepoint's family.
Decorative Wilcoxon companions reported.

(b) Pooled mixed models (primary home for the 4R-floor-proximity claim - NOT a
diagnostic). Two complementary `nlme::lme(..., random = ~1|experiment_label)`
fits:
  - ADDITIVE, reference = 4R: percent_adp_corrected ~ sample + timepoint. The
    sample coefficients are the suppressor-vs-4R differences POOLED across
    timepoints, with 95% CIs from `intervals()`. These are the reported
    floor-proximity estimates; a small estimate with a CI bracketing zero
    supports "near the 4R floor."
  - INTERACTION: percent_adp_corrected ~ sample * timepoint; `anova()` reports
    the sample:timepoint term (do the curves diverge - does WT pull away over
    time) and the sample main effect.
timepoint is a FACTOR (accumulation saturates  not linear; the numeric/linear
trend alternative was documented and not adopted). Both fits are wrapped so that
non-convergence is reported as a finding, not a crash.

Assumption status / caveat: per-cell replicate counts are unequal and WT cells
carry large variance while near-floor mutant cells carry ~0 (several exactly 0),
grossly violating the lme homoscedastic-residual assumption. Bartlett/Levene are
undefined on zero-variance cells and omitted. The homoscedastic lme is reported
WITH this caveat; `varIdent(form = ~1|sample)` is a documented, not-adopted,
heteroscedastic alternative. Because the pooled residual variance is inflated by
WT's spread, the suppressor-vs-4R CIs are, if anything, CONSERVATIVELY WIDE -
they do not overstate floor proximity.

Negative-value handling (DIFFERS FROM LOADING ON PURPOSE - no flooring): ATPase
is analysed by mixed model + per-timepoint DIFFERENCES, both subtraction-based
and negative-safe, with NO division by a near-zero quantity (the only division,
ADP/(ADP+ATP), is bounded [0,1] BEFORE subtraction; there is no fold-over-4R or
percent-of-WT column). Small negatives are EXPECTED noise around a true zero in
the near-zero 4R/suppressor group, and that noise is data the model must absorb.
Flooring would erase exactly that variance AND bias the near-zero group mean
upward (not cleanly one-directional safe here, unlike loading). Therefore
negatives are DETECTED + REPORTED (every offending row written to
atpase_negative_values.csv) and kept RAW. A sanity tripwire stops the script
only on an implausibly large negative (< -0.15, mirroring the WT-t=0 background
threshold), which would signal an upstream sign/load error; the observed
negatives (max magnitude 0.033) pass silently.

--------------------------------------------------------------------------------
## References
- Holm S (1979). A simple sequentially rejective multiple test procedure.
  Scand J Stat 6:65-70.
- Pinheiro JC, Bates DM (2000). Mixed-Effects Models in S and S-PLUS. Springer.
  (the `nlme` package.)
- Lakens D (2017). Equivalence tests (TOST). Soc Psychol Personal Sci 8:355-362.
  (the considered-and-rejected equivalence alternative for the ATPase claim.)
- Blocked/paired design: Snedecor GW & Cochran WG (1989), Statistical Methods,
  8th ed. (randomized complete / incomplete block designs; paired comparisons).
