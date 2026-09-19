# Two-half orientation at n=9: zero-error diamond fails

**Status:** finite CP-SAT certificate for one exact minimizer of order 9.
Not a doubling proof, not nonconvergence, original limit OPEN.

**Object.** For an optimal signing `A` of order `n`,
`B(A,R)=max_{x,y}(|Q_A(x)+Q_A(y)|+|x^T R y|)`. The zero-error form of
(6.13) is `min_R B(A,R) ≤ 2√2 Φ(A)`. Equivalently (Prop 6.5a)
`min_S max_U Φ(A^{F_S(U)}) ≤ √2 Φ(A)`. This is the directed-half-cut
reformulation of multiplier two, not a new target.

**Reused, not rerun.** Orders 5–8 were already solved in
`scripts/original_mo_two_half_geometry.py` /
`evidence/original_mo_two_half_geometry.json`:

| n | m_n | min B | 2√2 m_n | (B−2√2 m)/n^{3/2} | zero-error |
|---|-----|-------|---------|-------------------|------------|
| 5 | 4 | 16 | 11.31 | +0.419 | fail |
| 6 | 5 | 18 | 14.14 | +0.262 | fail |
| 7 | 9 | 22 | 25.46 | −0.187 | pass |
| 8 | 10 | 28 | 28.28 | −0.013 | pass |

**New.** An ILS exact witness of `m_9=12` (Phi replayed 12) and a complete
CP-SAT model on its 36 orientation bits. Minimize with 86 workers for 600 s
gave FEASIBLE `B=36`, bound 32. Feasibility then:

- `B≤32` **INFEASIBLE** (183.6 s, 2375 conflicts)
- `B≤34` **INFEASIBLE** (179.6 s, 2160 conflicts)
- `B=36` attained by ILS, independently replayed by `analyze_orientation`

Hence `min_R B(A,R)=36 > 2√2·12 ≈ 33.941`. Zero-error doubling **fails**
at this exact order-9 minimizer. Normalized excess `+0.076`.

Receipts: `evidence/two_half_n9_n10_20260919/{A9_exact.npy,n9_cpsat.json,n9_feas32.json,n9_feas34.json}`.
Backend: CP-SAT `num_search_workers=86`. GPU unused.

**What this does not show.** Paley-R is still the wrong construction
(yesterday's census). An A-dependent orientation exists and is optimal at
n=7,8 for the zero-error cut, then fails again at n=9. The signed excess
does not appear to be eventually negative. Dini error `n^{3/2}` still
covers the n=9 gap of 2, so this is compatible with (6.13) as stated.
It blocks any hope that “exact small minimizers already satisfy the
zero-error diamond for all n≥7.”

n=10 was not completed (ILS inner loop too slow; not a CP-SAT run).

**Next implication.** Either a Dini-only (not zero-error) bound on the
excess for large optimizers, or an optimizer switching certificate that
controls `2 min(|F|,|G|)` on opposite-sign halves (6.14a3) uniformly.
Do not rerun n=5–8.
