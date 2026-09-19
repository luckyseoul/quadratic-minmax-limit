# Residue (6.20) finite census

**Status:** finite observation. Not a proof of multiplier two, not emptiness
of (6.20) for large optimizers, and not a change to the OPEN original limit.

**Implication checked.** For each tested signing, is the unshielded set
(6.20) empty, and if not, does the explicit Paley-skew `R` of Proposition
6.6 still satisfy the mixed-state diamond
`|Q_A(x)+Q_A(y)| + |x^T R y| <= 2√2 Φ(A)` on that set?

**Backend.** One NumPy process per signing, `ProcessPool` over 10 signings,
`OMP_NUM_THREADS=1`. GPU unused (driver down). Inner work is an
`8192×8192` GEMM at order 14. Wall 5.6 s. Script:
`scripts/residue_620_census.py`. Receipt:
`evidence/residue_620_census_20260919/census.json`.

Constants as in (6.17): `a=(√2-1)/π`, `ρ=a²≈0.01738`. Cube folded at
`x_0=+1`. Paley-skew `R` is the compressed tournament of the next prime
`q≡3 (mod 4)` with the degree-balancing transfers of Prop 6.6, conjugated
by a maximizer.

## Results

| signing | n | Φ | α | residue fraction | max slack / n^{3/2} on residue | Paley-R holds on residue? |
|---------|---|---|---|------------------|--------------------------------|------------------------------|
| C5-cycle (exact `m_5=4`) | 5 | 4 | 0.358 | 0.547 | +0.419 | no |
| C6 conference (exact `m_6=5`) | 6 | 5 | 0.340 | 0.908 | +0.535 | no (energy threshold vacuous) |
| exact `m_8=10` | 8 | 10 | 0.442 | 0.508 | +0.341 | no |
| C10 conference | 10 | 15 | 0.474 | 0.223 | +0.492 | no |
| C14 conference | 14 | 21 | 0.401 | 0.671 | +0.661 | no |
| random seed 1 | 8 | 14 | 0.619 | 0.015 | −0.159 | yes |
| random seed 2 | 8 | 16 | 0.707 | 0.0006 | −0.409 | yes |
| random | 11 | 25 | 0.685 | 0.0005 | −0.129 | yes |
| random | 12 | 28 | 0.674 | 0.0014 | −0.173 | yes |
| random | 13 | 28 | 0.597 | 0.0046 | +0.102 | almost (0.03% of residue) |
| random | 14 | 41 | 0.783 | 0.00002 | −0.572 | yes |

At these orders `ρ n < 1`, so the Hamming-from-maximizer cut is essentially
“not the maximizer”. The live filter is joint energy versus
`2√2 Φ − n^{3/2}`. That threshold is vacuous whenever
`α ≤ 1/(2√2)≈0.353` (the low-α regime of Prop 6.6); C6 is in that regime.

## What this does and does not show

1. On every **exact minimizer** tested (`n=5,6,8`), residue (6.20) is a
   majority of pairs, and Paley-`R` misses the diamond by a
   `Θ(n^{3/2})` amount. So the leftover set is not an empty finite
   artifact of the geometric cuts.
2. On **high-α random** signings the residue is sparse and Paley-`R`
   typically holds. That is the wrong ensemble for the original minimum.
3. Conference constructions used as *upper-bound* matrices behave like
   the exact small minimizers: large residue, Paley-`R` fails. They are
   not claimed to be asymptotic minimizers.
4. This does **not** close residue (6.20), does not produce a doubling
   estimate, and does not change `liminf`/`limsup`. It shows that any
   closure for actual minimizers must use `A`-dependent anticorrelation
   of `Q_A` with `x^T R y` on a large Hamming-central set; the global
   estimate (6.24) is the gap the theorem already named.

The 2026-09-18 cross-transfer sweep (rescued into
`evidence/cross_transfer_20260918/`) measured a different 2n object at
n=13 and is not a substitute for this census.

## Next unresolved implication

A structural bound
`|x^T R y| ≤ 2√2 M − |Q_A(x)+Q_A(y)| + n^{3/2} Ω(n)`
on (6.20) for an optimal `A`, or a proof that (6.20) is empty for all
sufficiently large optimizers. The finite exact points above make
emptiness at small n false. Multiplier three / 1:2 split remains a
separate open ray.

OpenAI referee `suggest_direction` (gpt-5.6-sol, 2026-09-19), after the
census: (1) optimizer switching inequalities on a central pair, isolating
`|x^T R y|` against joint energy; (2) dualize
`|x^T R y| − c E_A(x,y)` over Hamming-central pairs; (3) dyadic amortized
potential so doubling residues telescope even if (6.20) stays nonempty
pointwise. Avoid A-independent norms, more Paley/Hadamard/Gaussian
sweeps, and the disk `I²+C²`. These are directions, not theorems. Claude
was not consulted.
