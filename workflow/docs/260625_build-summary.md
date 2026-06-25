# lw -- build summary

Execution record for the `lw` rebuild, following the commit-by-commit plan
(`lw_build_plan`) against the v2 handoff packet (`lw_handoff_v2`). Section
references (SS) point into the v2 packet.

**Final state:** `lw.py` (the tool), `test_lw.py` (29 tests, all passing),
`README.md`. Two tables, five commands plus `init`. Pure-transform / IO-edge
boundary intact throughout.

---

## Commit log (in order)

| # | Commit id | Tier | Scope | What landed | Refs |
|---|---|---|---|---|---|
| 1 | C1-C3 | haiku | core | Data definitions, config, paths; effects-as-data primitives; id normalization. The inert dependency roots. | SS6.1, SS6.2, SS6.3, SS6.6 |
| 2 | C4-C5 | sonnet | store | Two-table schema + `schema_notes`; `load_store_from_database`; whole-store `write_store_to_database` with the union-of-keys / `row.get` fix. | SS1, SS6.4 |
| 3 | C6-C7 | sonnet | store | Missing-DB precondition (returns failure effects); idempotent `init` that never resets the counter. | SS2, SS6.6, SS6.7 |
| 4 | C8-C9 | sonnet | new | Counter edge (read + bump-last); pure `new` transform returning the next counter value as an extra effects key. | SS6.2, SS6.6 |
| 5 | C10 | haiku | new | `new` edge wiring: mkdir + commit first, counter bump last; completes the decomposed opus unit. | SS6.6, SS6.7 |
| 6 | C11-C12 | sonnet/haiku | link | Bare-edge `link` transform with four no-effect rejection paths; wiring. | SS4, SS6.5, SS8 |
| 7 | C13-C14 | sonnet/haiku | query | `show` (unconditional listdir at edge, missing-folder rendered not raised) + `list` (type/status filters). | SS2, SS3, SS5 |
| 8 | C15 | sonnet | cli | `if/elif` dispatch, per-command parse, missing-DB precondition, `command_log` audit effect on success. | SS5, SS6.6, SS8 |
| 9 | C15a | -- | all | Naming + control-flow refactor: full descriptive names, no abbreviations, explicit control flow. Behavior byte-for-byte identical. | SS6.1 |
| 10 | C16 | sonnet | test | 17 tests: per-command happy / error / no-effect, `link`'s four rejections mandatory. Mutation-checked. | SS9 |
| 11 | C17 | -- | harden | Counter-file guard (`CounterFileError`); counter/DB collision precondition; magic values surfaced to constants; broken-pipe guard; actionable stderr. +5 tests. | SS3, SS6.4, SS6.6 |
| 12 | C18 | -- | cli | Descriptive usage / dispatch errors (condition -> usage -> next step) from one declarative table; `list` stray-positional caught. +7 tests. | SS5, SS8 |

---

## Tier outcome (as planned)

The plan predicted exactly one opus-tier unit (`new`), decomposed before
writing into two sonnet children (C8 counter edge, C9 pure transform) plus
one haiku child (C10 wiring). That is what happened: no irreducible opus
commit survived, which is the healthy outcome for a two-table / five-command
tool. Everything else was a thin trunk of forced order with four
parallelizable branches (the C8/C9 pair, C11, C13, C14) reconverging only at
dispatch (C15).

---

## What changed from the plan during execution

Three additions beyond the original C1-C16 plan, each at your direction:

1. **C15a (refactor).** Applied the descriptive-naming / explicit-control
   -flow / no-abbreviations conventions across the whole file as a pure
   refactor, verified identical via the full CLI smoke test.

2. **C17 (hardening).** Five agreed items (counter-file corruption guard,
   counter/DB collision precondition, magic-value extraction, flow-stage
   stderr hints, critical invariant comments) plus one real defect the
   testing surfaced (a `BrokenPipeError` traceback when piping output into
   `head`). Three further hardening ideas were deliberately **skipped** as
   out of scope and recorded in the README's Known sharp edges: input
   sanitizing, transactional folder/row/counter coupling, and
   concurrent-writer locking.

3. **C18 (CLI polish).** Replaced the stub usage strings with descriptive
   condition -> usage -> next-step messages, and caught `list`'s previously
   silent stray-positional case.

---

## Carried verbatim (not rewritten)

Per SS6, these were carried as proven and correct rather than redesigned:
the effects-as-data primitives and their contract (the dict keys stay
abbreviated as the carried contract even after the naming pass), the
folder-naming convention, id normalization, and the whole-store load/commit
with its union-of-keys fix.

---

## Test coverage (29 tests)

- **Per command:** happy path, error path, no-effect assertion.
- **`link`:** all four rejection paths (bad id, missing experiment, unknown
  relationship, duplicate edge) as distinct no-effect cases, plus
  reverse-direction-is-not-a-duplicate.
- **`show`:** summary + full listing, missing-folder rendered-not-raised,
  missing-experiment no-effect.
- **Hardening:** corrupt / empty / missing counter file; counter-DB
  disagreement protecting the older row; id-constant coupling.
- **Dispatch:** every usage / error message asserted for condition and
  next-step content.

The no-effect assertions were mutation-checked: injecting a counter bump on
the `new` failure path made exactly the right test fail, confirming the
effects-boundary assertions are not vacuous.
