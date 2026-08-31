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

## Current residual-(ii) gate after Proposition 15.745

Do not reopen `p=13,k=58` or `p=17,k=74`.  Proposition 15.742 combines the
`M_2=0` congruence with six multiplicative interval cuts and closes the
former row by the exact energy contradiction `667<707+26C`.  Proposition
15.743 extends the same common difference-Radon mechanism to the latter.
At `p=17`, the local mean formula and the global signed total are both
needed: together they force `P_L=4+k_L` in every hard direction.  Do not
silently assume that identity in a local cell model, and do not use the
normalized exact-star row `q=(2)^8` to obtain its own normalization.  First
glue the unspecialized exact-row sums to get `hT=18P_L-69`; common `hT`
makes the exact-star `P_L` common, and then `6P_L<=75` together with the
isolated-chart congruence `P_L=5 (mod 8)` forces `P_L=5`.  Only afterward
may one set `hT=21` and `q=(2)^8`.  The complete catalog of 698
translation-averaged nine-set cuts then makes hard excess one
infeasible, gives sharp excess-two and excess-three energies 70 and 119,
and forces every opposite row to `(-3)^8`, of energy 72.  The only partition
not already killed by the excess-one row satisfies
`119+9*72=767<1211+34C`, so it is impossible as well.

At `p=13,t=4,k=60`, Proposition 15.744 replays all residues and closes the
exceptional `u=3` profile by a six-root quartic contradiction after rebuilding
the edge-count-sensitive mass-14 models at `|H|=61`.  Its `b=10` premise is
also exact: a rank-78 restriction promotes contact-layer equality to the
pointwise complement triple, and a separate 1,716-variable punctured-lift
model excludes the two-unit `b=10` cell.  Do not replace that model with
Proposition 15.688, because the difference can be negative on the omitted
intersection layer.  Proposition 15.745
closes `u=0`: the 74-cut row bounds and common Radon energy give `C<=1` in
the sole difficult partition, while its seven parallel edges in six classes
give `C>=1`; its unique doubled parallel displacement then makes
the transverse multiplicities Boolean and bounds the elevated row in
`[-7,6]`, giving `695<719`.  Do not import the old `|H|=59` height-four
infeasibility, omit the collision-one sign audit, or call the full
`p=13,k=60` row closed.  Its exact remaining residues are `u=4,6`.

Residual (ii) is still open at critical `p=5,7`, at `p=11,k>=50`, in the
two `p=13,k=60` residues `u=4,6` and all later p13 layers, at `p=17,k>=76`,
and in later layers for primes `p>=19`; within the branch-B fourth shell, the
next prime is `p=29`.  Propositions 15.743--15.745 are finite certificates,
not an all-prime row theorem.  The preferred generic front is therefore a
structural version of its common-energy/cut mechanism that survives when
the number of distance bins grows.  Another independent coefficient-cell
catalog, one-direction floor, halving heuristic, or longer complete-domain
timeout does not advance that gate.  The next finite p13 implication is the
support-330 sharp mass-ten equality classification on `J(13,7)`, reusing
15.738's third-difference annihilator; it feeds `u=4` and partially prunes
`u=6`.  Its expected candidates are the 78 omitted-pair and 286 all-equal
triple supports, with only 70 anchored no-goods.  This is not a census of the
already-closed `u=0,3` rows.  Before launching a finite p11/p13-later
computation, identify the invariant that could extend beyond that one row or
explain why the finite row is a genuine base obstruction.

## Result discipline

Label every result as exactly one of: proved theorem, exhaustive finite
certificate, open reduction, counterexample, or retracted claim.  Never
promote computation or a heuristic pattern into a theorem.

After a genuine advance, update the proposition-dedup audit and the canonical
status/handoff documents in the same commit.  Preserve failed routes and
counterexamples: they are part of the project memory, not clutter to delete.
