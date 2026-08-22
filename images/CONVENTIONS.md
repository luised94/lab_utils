Coding conventions - lab_utils
Self-contained. Applies to everything in this repository. Written to be pasted
into a fresh context with no other background.
---
1. Paradigm
Data oriented. Data in, data out.
Flat procedural. No classes. No OOP. No clean-code patterns, no design patterns,
no agile ceremony. Plain data types and plain data structures.
The unit of organization is a readable block of inline logic, not a function.
No abstraction built for a single call site
No helper functions, wrappers, adapters, or indirection that exists to serve one
caller. Inline the logic where it is used.
Add abstraction when the second real call site exists. Not in anticipation of
one. Anticipated need is not need.
No indirection or nesting that does not pay for itself
One longer readable block beats three short ones that force jumping around to
follow the control flow. Scrolling is cheaper than chasing.
Permitted exceptions
Only these:
A tagged message emitter (writes to stderr, accumulates into the run log)
A fail-fast error exit
Both have many call sites, so both pay for themselves. Nothing else gets a
function without an argument for it.
---
2. Naming
Full descriptive complete names carrying domain information, and unit
information where applicable.
No abbreviations
The exception is tokens that behave as nouns in ordinary use - `usb`, `tiff`,
`json`, `csv`, `pmt`. These are words, not abbreviations. Everything else is
spelled out.
No single-letter names
Loop variables included. There is no exemption for `i`, `j`, `x`, `y`, `n`.
Names carry domain and unit
```
WORKTREE_ROOT                    not   ROOT
HEAD_OBJECT_NAME                 not   OBJ
RESOLVED_ANCESTOR_DIRECTORY      not   PHYS
lane_pitch_millimetres           not   pitch
integrated_intensity_counts      not   intensity
rolling_ball_radius_millimetres  not   radius
```
If a quantity has a unit, the unit is in the name. This is not decoration - it
is the mechanism that prevents a millimetre value being passed where a pixel
value is expected.
Different things get different names, across scopes
If two variables hold different things, their names must say which - even when
they live in different scopes and could technically share a name without
colliding. The reader should never have to determine which `count` is meant by
tracing scope boundaries.
Renaming
Never trust blind pattern substitution. Before renaming, check for the token
in:
prose and comments
user-facing messages
ticket and issue references
output format strings
filenames and directory names
documentation
Re-run the tests after.
---
3. Comments
Comments state why, not what.
A comment earns its place by recording one of:
the constraint that forced this
the failure it prevents
the alternative that was rejected, and why
a non-obvious platform or format detail
```
# Typhoon writes .img with origin at bottom-left, so row 0 in the file is the
# bottom of the gel. Without this flip, band row indices invert while lane
# ordering stays correct - producing entirely plausible numbers with every
# band label wrong.
```
Not:
```
# flip the array vertically
```
If the line's behaviour is unclear without a comment describing what it does,
the fix is usually a better name, not a comment.
---
4. Character set
ASCII only. In code, comments, output messages, and filenames.
No smart quotes, no em dashes, no Greek letters, no arrows, no box-drawing
characters. Write `micrometres`, not the symbol. Write `->` if an arrow is
genuinely needed in a message.
---
5. Scope discipline
State what is out of scope and stay inside it.
If something worth doing falls outside the agreed scope, name it and ask. Do
not fold it in because it was convenient to do while nearby.
Flag anything sitting on the boundary and get a decision before proceeding.
Boundary cases are where scope creep enters, and they are recognizable in
advance.
---
6. Verification and honesty
Verify by executing, not by asserting. Show the output.
A claim that the code works is not evidence that the code works. Run it. Paste
what it printed.
If a planned step turns out to be unnecessary or already satisfied, say so.
Do not manufacture a diff to fill the slot in the plan.
Report failures and mistakes found in your own work, including in work
already delivered. Own them plainly, state what was wrong and what it affected,
and fix it. No apology spiral, no self-flagellation - the correction is the
point.
Prior comments and documentation are evidence, not ground truth. Where code
does not do what it claims to do, say so and say where.
---
7. Working method
Before writing any code
On a new task, respond first with:
what you understand the task to be
the proposed strategy
open questions
Do not write code in that turn. Wait.
Read the code you were given before proposing anything.
After a plan is agreed
Run an adversarial pass over the plan before implementing:
what has not been considered
what would make this go wrong
what is already broken that this change will expose
Then implement.
---
8. Options and recommendations
Never present an unranked list.
When presenting alternatives:
rank them
score them
give a one-line rationale for each
end with one recommendation
Explain the reasoning behind design choices, including the ones rejected and
why. Name the tradeoff, not only the conclusion.
A recommendation without the cost of taking it is incomplete.
---
9. Delivery
Commits
Commit by commit where history matters. One concern per commit.
Each commit independently valid: it parses, it runs, the series bisects.
Commit messages state why, not what. The diff already says what.
Note explicitly when a commit is comment-only or a no-op.
Format
Prefer a patch series that can be applied over pasted code.
Include:
a checksum of the baseline patched against
the final file as a fallback
---
10. Python specifics for this repo
Managed with `uv`
Scripts self-contained via PEP 723 inline metadata (`# /// script`)
Lock or requirements file written once dependencies are settled
Messages to stderr, tagged with their source, so stdout stays clean for
piping
All emitted messages also accumulate in memory and are serialized into the
provenance JSON - the run log is part of the record, not lost scrollback
Config as a flat block of plain data at the top of the file
Physical lengths specified in millimetres in config, converted to pixels using
the pixel size read from the input file's own header
---
11. ImageJ macro specifics (.ijm)
The macro language does not fit the Python rules literally; apply the intent, not
the letter.
Descriptive names still apply: `source_width` not `sw`, `copy_coordinates` not
`coords`. The no-abbreviations rule is language-independent.
Comments still state why, not what: why batch mode, why copy-once-before-paste,
why a mode was removed.
The language has no local scoping and no typed data structures; top-level globals
and the `newArray` flat-array idiom are unavoidable and are not a style violation.
No abstraction for a single call site holds here too. `distribute_selection.ijm`
dropped its Grid layout because gel lanes tile in one Row (or Column); Grid had no
real call site and its perfect-square path placed a copy on another. Row + Column
are the two real uses.
ASCII only, like everything else.
Macros are run interactively in ImageJ/Fiji on the Windows side; they are not
exercised by the Python regression harness. A change to a macro is verified by
running it on a real image, and its placement math can be checked by porting the
coordinate function to a scratch script (that is how the Grid bug was found).
