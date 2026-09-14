# Session audit trail: gluing experiments, convention correction

**Status:** 2026-09-13. Audit record for this session's gluing work. The
mathematical content is in
`NOTE_2026-09-13_MATCHED_BLOCK_GLUING_BOUNDS.md`; this file records the
process corrections.

## 1. Duplication check (user-requested, second pass)

The gluing ("min over cross completions") mechanism duplicates the recorded
completion route: `NOTE_2026-09-12_COMPLETION_DISCREPANCY_FRAMEWORK.md`
(prior-art notice), `NOTE_2026-09-11_ORDER16_HUNT.md`, and the recorded
identity in `NOTE_2026-09-01_ORIGINAL_LIMIT_TWO_RAY.md` eq. (4). The audit
was repeated across all branches and worktrees in this session; no other
copy of the matched-block computation was found. The sub-family differs
from the recorded fixed-`A8` completion family: both diagonal blocks are
chosen (not fixed to one source), which is what allows the family to attain
the merged order's optimum at N=8,9,10.

## 2. Convention error found and corrected in-session

An intermediate kernel (`glue_literal.py` and derivatives) evaluated split
sums over a reduced `2^(N-1)` state set. The reduction is valid for a single
quadratic form (complement symmetry) but not for a sum whose two terms are
not separately symmetric; and the merged cross block is directed. Result:
one verification fragment reported a value (19) inconsistent with the
literal loop and the full-cube kernel (both 13) for the same 5+5 witness.
All quoted numbers were then recomputed with `glue_exact_raw.py`, which
evaluates the full cube and was validated against `m_2..m_6 = 1,3,4,4,5`
and against a literal double loop. The earlier JSONs and helper scripts
were removed from the tree; `raw_glue_*.json` supersede them.

## 3. What stands

- Exact matched-block gluing minima: 10 (4+4), 12 (4+5), 13 (5+5), each
  equal to the recorded `m_8, m_9, m_10`; full enumerations 2^16, 2^20,
  2^25 completions.
- Counts of optimal completions: 184, 10800, 1840.
- Proposition 1 and Proposition 2 of the note (both-sided bound; the
  recorded cross floor does not lower-bound the gluing minimax).
- No convergence claim; the doubling estimate remains open.
