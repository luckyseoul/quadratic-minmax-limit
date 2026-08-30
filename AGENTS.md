# Research continuity rules

This repository is a long-running proof project.  Avoid rediscovering old
branches or turning completed computations into new work.

## Before starting an attack

1. Read the current gate in `STATUS.md` and `HANDOFF.md`.
2. Search `evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md`, `src/`, `tests/`,
   and `evidence/` for the proposed object, parameter range, invariant,
   script family, and expected output.
3. State the one unresolved implication the proposed work would close and
   the existing result it strictly advances.
4. Do not launch the work if its distinguishing output is already recorded.

## Duplication gate

- A rerun is allowed only after identifying a concrete changed premise: a
  code defect, corrected theorem hypothesis, new parameter range, or new
  invariant.  Record that delta before running it.
- Do not rerun finite-prime, CP-SAT, MILP, SAT, orbit, OEIS, or literature
  searches merely to reconfirm a result.  Reuse their stored artifacts.
- If the canonical documents disagree, repair the contradiction before
  spending compute.
- Keep one active mathematical gate.  Closing a subcase does not authorize
  appending another search list; update the gate and reassess the proof cold.
- Prefer an argument that closes an infinite family over a wider finite
  census unless the finite census is the explicitly recorded gate.

## Result discipline

Label every result as exactly one of: proved theorem, exhaustive finite
certificate, open reduction, counterexample, or retracted claim.  Never
promote computation or a heuristic pattern into a theorem.

After a genuine advance, update the proposition-dedup audit and the canonical
status/handoff documents in the same commit.  Preserve failed routes and
counterexamples: they are part of the project memory, not clutter to delete.
