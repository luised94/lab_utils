# Implementation handoff: statistical analysis for ORC4R suppressor screen

This document hands off a fully-planned, adversarially-reviewed,
topologically-sorted implementation to be carried out across THREE sequential
threads. It is self-contained: the implementing thread has NO repo access
(files are attached, not cloned), so all cross-file blast-radius is
pre-computed here. Do not assume the implementer can grep siblings.

Style contract (non-negotiable): scripts stay PROCEDURAL. Configs/flags at top,
preprocessing, derivation, main logic, summary, output. NO helper functions, NO
abstractions, NO premature interfaces. Explicit descriptive variable names with
domain information. "Nothing too clever." Data in, data out.

Dependency contract: base R `stats` for everything except mixed models; `nlme`
(ships with R) for mixed models. NOTHING hand-rolled. renv manages packages.

Annotation contract: plots show EXACT p-values, never stars.

================================================================================
PROJECT FACTS (corrected and confirmed -- treat as ground truth)
================================================================================

F1. Four files, three logical analyses:
    - orc4r-screen_loading-kglut-titration.R     (loading, WT/ORC4R/+4sofa)
    - orc4r-screen_loading-all-suppressors-350mm.R (loading, all 5 supps + WT + 4R)
    - 260411_orc4r-screen_analyze-sofa-atpase.R   (ATPase compute)
    - 260411_orc4r-screen_plot-sofa-atpase.R      (ATPase plot, consumes above)

F2. Paths: scripts point at the data's CURRENT real location via
    MC_DROPBOX_PATH. For Zenodo, user copies data + script into the deposit
    repo where they are co-located, so the script must ALSO resolve
    script-relative. Resolution chain: script-relative -> MC_DROPBOX_PATH ->
    stop() naming BOTH attempted paths. Support ONLY source() invocation
    (not Rscript, not interactive paste); error clearly otherwise. Both the
    original location and the Zenodo copy remain live.

F3. ALL loading AND ATPase assays are paired/blocked: WT, 4R/ORC4R, and the
    mutant(s) are co-run within one replicate/experiment. Block = replicate
    (loading) or experiment_label (ATPase). WT and 4R appear in nearly every
    ATPase experiment (why they accumulate more replicates) -> the block
    structure is UNBALANCED and INCOMPLETE for ATPase (each mutant in only its
    own ~3 blocks).

F4. kglut titration = ONLY WT/ORC4R/+4sofa, n=3, salts 250/300/350. n=1
    rotating-suppressor rows are filtered out, never enter any test.

F5. ATPase null framing is DELTA-FREE. All suppressor complexes carry the same
    catalytic mutation as 4R, so no restoration of hydrolysis is
    mechanistically possible. Claim = conjunction of two ordinary claims:
    (a) "not restored to WT" = WT-vs-suppressor SIGNIFICANT (rejectable);
    (b) "near 4R across all timepoints" = (suppressor-4R) difference+CI,
    POOLED across timepoints via the mixed model (NOT per-timepoint, to avoid
    cherry-pick exposure). No equivalence margin delta is used.

F6. Loading restorers = Orc1, Orc4. Reducers = Orc3, Orc5, Orc6. 4R < WT.
    Comment-level only; all five suppressors get identical statistical
    treatment. CRITICAL: the loading (suppressor-4R) comparison is an ORDINARY
    effect size -- reducers fall BELOW 4R by design. This is NOT the ATPase
    floor-proximity argument. Do not conflate the two narratives.

================================================================================
GLOBAL INVARIANTS (every commit inherits these -- restate in each thread)
================================================================================

I1. "Green" for this project = the script source()s end-to-end without error
    AND its in-script stopifnot/assertion guards pass against real data. There
    is no separate test runner; the assertions ARE the tests.

I2. Holm correction family = {WT-contrasts within a SINGLE salt} or
    {WT-contrasts within a SINGLE timepoint}. NEVER across salts or across
    timepoints -- those are separate questions, not one family. Every test
    commit cites this.

I3. NORMALIZATION FIREWALL (kglut): per-kglut normalized values pin WT=100 at
    every salt, so the salt-sensitivity interaction computed on them is
    mechanically degenerate. The salt x genotype interaction MUST use
    percent_wildtype_global (WT-at-250mM baseline). C13 carries a runtime
    assertion that fires if the per-kglut column is fed to the interaction.

I4. n=3 HONESTY: every parametric test carries an explicit in-comment caveat
    that normality cannot be meaningfully tested at n=3; justification rests
    on assay track record + paired design, not on a passing normality test.
    Wilcoxon companions are labeled DECORATIVE (n=3 paired floors two-sided
    p ~ 0.25), included for transparency, never leaned on.

I5. MIXED MODELS ARE PRE-COMMITTED DIAGNOSTICS where a paired test is the
    reported result (C9). The reported result is fixed in advance; the model
    is a consistency check; disagreement triggers investigation, NEVER
    selection of the friendlier p.

I6. DECISION LOG: each thread appends its decisions (test chosen, assumption
    status, correction scope, modeling choices like factor-vs-numeric) to a
    running DECISION_LOG section and pastes it forward. The final docs
    (C19/C20) are CONSOLIDATION of this log, not reconstruction.

================================================================================
THREAD SPLIT (3 sequential threads, ~23 commits)
================================================================================

Thread boundaries chosen so each holds its invariants accurately (~6-10
commits each) and break at natural checkpoints.

--------------------------------------------------------------------------------
THREAD 1 -- Infrastructure + bug fixes + nlme  (C1-C5, C0d)  [6 commits]
--------------------------------------------------------------------------------
ATTACH: all four .R files (C1, C2 touch all; C3 edits atpase-analyze; C4,C5
        edit loading scripts). Attach all four so cross-file config (PLOT_CONFIG,
        factor_order, SAMPLE_COLORS) is visible even where not edited.

C1  SONNET  all 4   Dual-location path resolution, source()-only.
    Edge: none upstream. Blast radius: every data-load in every script.
    Green: each script source()s and finds data in both homes; stop() lists
    both attempted absolute paths when neither resolves.

C2  HAIKU   all 4   Embedded git-state reference comment block (manual-fill).
    Green: scripts source() unchanged; comment present. No runtime git calls.

C3  SONNET  atpase-analyze  library(xlsx)->readxl; read.xlsx(header=FALSE,
    sheetIndex=) -> readxl::read_excel(col_names=FALSE, sheet=).
    BLAST RADIUS (pre-mapped): the imagej_index integer check (~line 357-364)
    and even-row pairing (~line 366-372) and sequential-index assert. readxl
    returns a tibble (index may parse as double) -- the integer check MUST
    still fire. Re-confirm sheet indices match old sheetIndex numbering.
    Green: all sheets load; integer/pairing/sequential asserts pass.

C4  HAIKU   350mm   Fix guard-ordering crash. Move derived-calc mutate block
    (current lines 257-270) ABOVE the Inf/NaN guard (current 246-255).
    Relocate premature "computed" message (current 244). Welded-atomic: script
    does not run until both halves land.
    Green: script source()s end-to-end; Inf/NaN guard passes on now-existing
    columns.

C5  HAIKU   kglut,350mm  Cosmetic message-ordering cleanup (duplicated/premature
    message() calls). Green: scripts source(); messages logically ordered.

C0d HAIKU   (renv)  Restore + record nlme in the project library (EXTRACTED
    from the old C21 because C9/C13/C17 consume nlme before final lockfile).
    Green: library(nlme) loads.

>>> THREAD 1 ENDS AT THE RUN CHECKPOINT <<<
USER ACTION after Thread 1: run all four scripts. PASTE BACK:
  - full console output of each script (the message() trail + any error)
  - the existing output CSVs (summary + full-data) so Thread 2 sees real
    data shape and real values
  - confirmation that nlme loads
This output is REQUIRED input to Thread 2 (stats depend on real data shape).

--------------------------------------------------------------------------------
THREAD 2 -- Loading statistics: 350mm + kglut  (C6-C14b)  [~10 commits]
--------------------------------------------------------------------------------
ATTACH: both loading .R files (as updated by Thread 1), the DECISION_LOG from
        Thread 1, and the CSV outputs pasted at the Thread-1 checkpoint.
        Heaviest statistical thread; carries the normalization firewall (I3).

Phase 1 (350mm), depends on C4:
C6   SONNET  Assumption diagnostics + emission (residuals, QQ, Shapiro,
     variance check). I4 caveat in comments. Green: diagnostics CSV written.
C7   SONNET  WT-contrasts: paired t-test (block=replicate) WT vs each of 5
     suppressors AND WT vs ORC4R; Wilcoxon companions (I4 decorative); Holm
     within WT-family (I2). Exact p retained. Green: contrasts computed,
     Holm applied within family.
C8   SONNET  (suppressor-4R) ordinary effect size + CI. MANDATORY COMMENT
     (F6): NOT the ATPase floor argument; reducers below 4R by design.
     Green: CIs computed; comment present.
C9   HAIKU   Mixed-model consistency diagnostic (nlme, ~1|replicate).
     PRE-COMMITTED diagnostic (I5). Green: fit converges; flagged diagnostic.
C10a SONNET  Augment data+summary+results CSVs (estimate/CI/raw-p/Holm-p/
     test/n columns). Green: columns present in CSVs.
C10b SONNET  Plot: exact-p annotations, WT-contrast family. Independent of
     C10a (different file). Green: ggsave succeeds; p-values render.

Phase 2 (kglut), depends on C5 + global-baseline column already in code:
C11  SONNET  Assumption diagnostics (per-salt + pooled residuals).
C12  SONNET  Per-salt paired tests: WT-vs-ORC4R (direct reviewer answer),
     WT-vs-+4sofa, ORC4R-vs-+4sofa; Wilcoxon companions; Holm WITHIN each
     salt (I2). Green: per-salt contrasts; within-salt Holm.
C13  SONNET  Salt-sensitivity interaction. nlme: genotype*salt, ~1|replicate,
     ON percent_wildtype_global (I3 FIREWALL). Salt as FACTOR (no trend
     assumption); document the numeric/trend alternative.
     *** SANCTIONED DELIBERATE RED ***: first point the firewall assertion at
     the per-kglut column, observe it FIRE (red = guard works), then repoint
     to global-baseline, green. The red is evidence the guard works. This is
     the ONLY deliberate red in the whole implementation.
     Green: interaction fits on global-baseline col; firewall assertion passes.
C14a SONNET  CSVs: BOTH framings (per-salt + interaction). Label column
     lineage (which normalization feeds which claim). Green: both in CSVs.
C14b SONNET  Plot: per-salt WT-vs-ORC4R exact p primary; interaction p in
     bound caption text. Green: ggsave succeeds.

>>> THREAD 2 appends to DECISION_LOG and pastes it forward. <<<
Optional checkpoint: user may re-run loading scripts and paste outputs to
verify before Thread 3, but Thread 3 (ATPase) is independent of loading
results, so this is not strictly required.

--------------------------------------------------------------------------------
THREAD 3 -- ATPase statistics + documentation  (C15-C21)  [~9 commits]
--------------------------------------------------------------------------------
ATTACH: both ATPase .R files (analyze as updated by Thread 1's C3), the
        accumulated DECISION_LOG (Threads 1+2), and the ATPase output CSVs from
        the Thread-1 checkpoint. Carries the unbalanced-block subtlety (F3) and
        the pooled-vs-per-timepoint split (F5).

C15  SONNET  Assumption diagnostics (ATPase). Note unequal-n explicitly --
     affects variance assumptions. Green: diagnostics emitted.
C16  SONNET  Per-timepoint WT-contrasts. Block=experiment_label. WT SUBSET to
     each suppressor's SHARED blocks (F3: incomplete design -> a "WT vs Orc5"
     paired test uses only the 3 Orc5 experiments; WT's full replicate count
     does NOT strengthen these). REPORT per-contrast n explicitly. Holm
     within-timepoint (I2). Serves ONLY the "not restored to WT" half (F5a).
     Green: subsetting correct; per-contrast n reported.
C17  SONNET  Pooled mixed model -- PRIMARY home for the 4R-proximity claim
     (F5b). nlme: percent_adp_corrected ~ sample*timepoint, ~1|experiment_label.
     Report suppressor-vs-4R contrast + CI POOLED across timepoints. Also
     report sample:timepoint interaction (do curves diverge) + sample main
     effect. timepoint as FACTOR (accumulation saturates, not linear);
     document numeric alternative. Green: model fits; pooled contrast+CI
     computed.
C18a SONNET  atpase-analyze: write results CSV (per-timepoint tests + pooled
     mixed-model terms). Green: results CSV written.
C18b SONNET  atpase-plot: READ results CSV (NOT recompute -- anti-drift) +
     annotate per-timepoint exact p above each cluster (WT contrasts).
     Edge: depends on C18a. Green: plot reads CSV; per-timepoint p renders.

Phase 4 (docs), depend on the full accumulated DECISION_LOG:
C19  SONNET  STATISTICAL_METHODS.md (Zenodo, terse/defensive). One section per
     analysis: design, test, assumption status at n=3, correction scope (I2),
     delta-free justification (F5), exact tools (base stats + nlme). LABEL
     which normalization column feeds which claim. Describes what the code
     ACTUALLY did. Green: file written; covers all three analyses.
C20  SONNET  Self-study reference (expansive). Concepts, the affirming-the-null
     discussion, equivalence/TOST as considered-and-rejected alternative,
     primary sources (Holm 1979; Pinheiro & Bates nlme; blocked-design refs;
     Lakens on equivalence), self-test question bank. Green: file written.
C21  HAIKU   Finalize + verify deposit renv lockfile; bare-install run-through.
     Green: lockfile complete; clean-environment run succeeds.

================================================================================
ADVERSARIAL-REVIEW DELTAS ALREADY FOLDED IN (do not re-litigate)
================================================================================
- C16 uses shared-block subset; per-contrast n reported (Lens 1 consensus).
- C17 elevated to PRIMARY for pooled 4R-proximity (Lens 4 consensus); per-
  timepoint 4R-proximity is secondary, not headline.
- C9 constrained to declared diagnostic; paired test pre-committed (Lens 3).
- C13 gated on percent_wildtype_global with firewall assertion (Lens 2 / I3).
- C8 carries mandatory deviates-vs-floor comment (Lens 5 / F6).
- Holm family scope fixed globally (Lens 6 / I2).
- C0d extracted from C21 so nlme is available to its consumers (STEP-3
  blast-radius finding: C21 was serving two meanings).

================================================================================
DECISION_LOG (append-only; seed entries -- threads add to this)
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
[Thread 2 appends here ...]
[Thread 3 appends here ...]
