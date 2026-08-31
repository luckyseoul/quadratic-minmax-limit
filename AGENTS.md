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

## Current residual-(ii) gate after Proposition 15.742

Do not reopen `p=13,k=58`.  Proposition 15.742 combines the `M_2=0`
congruence with six multiplicative interval cuts and proves sharp elevated
and opposite row energies 31 and 82.  Their total `667` contradicts the
common-graph lower bound `707+26C`, so neither the 84-class quartic solve nor
the binary midpoint lift remains live at that row.

Residual (ii) is still open at critical `p=5,7`, at `p=11,k>=50`, at
`p=13,k>=60`, and in the generic `p>=17` later layers.  The preferred next
serious front is a structural cross-direction theorem for generic
`p>=17,t=3`, where 15.739 already forces higher even moments, the opposite
entry alphabet `{-1,0,1,2,3}`, and at p17 every cut into `[-26,-12]`.
Another independent cell catalog, one-direction floor, halving heuristic,
or longer complete-domain timeout does not advance that gate.  Before
launching a finite p11/p13-later computation, identify the invariant that
could extend beyond that one row or explain why the finite row is a genuine
base obstruction.

## Result discipline

Label every result as exactly one of: proved theorem, exhaustive finite
certificate, open reduction, counterexample, or retracted claim.  Never
promote computation or a heuristic pattern into a theorem.

After a genuine advance, update the proposition-dedup audit and the canonical
status/handoff documents in the same commit.  Preserve failed routes and
counterexamples: they are part of the project memory, not clutter to delete.
