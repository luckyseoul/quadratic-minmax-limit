# Two-outlier near-pencil / three-spike phase gate

**Date:** 2026-09-02
**Status:** proved necessary parity gate and exact relative-alignment barrier;
residual (ii) remains open

This note joins, for the first time, the lower-endpoint two-outlier
near-pencil reduction with the explicit all-prime Kiss--Somlai all-bad
three-spike completion.  It is not a small-prime census and it does not claim
that a phase-compatible graph satisfies the full Max inequalities.

Put \(p=4r+1\).  At the lower branch-B endpoint,

\[
 |H|=h=4r^2+6r+5,
 \qquad Q=|H_{\rm opposite}|=2r^2+r,
 \qquad S_H(x)=8-3p.
\]

The target score has

\[
 N_-={h-(8-3p)\over2}=2r^2+9r
\]

negative edge features.  On the other hand, for every simple graph with
\(\partial H=D\),

\[
 (-1)^{N_-}
 =\prod_{uv\in H}C_{uv}x_ux_v
 =(-1)^Q\prod_{v\in D}x_v.
\]

Both \(N_-\) and \(Q\) have parity \(r\).  Therefore every common graph at
this endpoint must satisfy the exact phase condition

\[
 \boxed{\prod_{v\in D}x_v=+1.}
\tag{1}
\]

This already explains why the first aligned deterministic \(p=53\)
matching-plus-\(C_7\) construction could not be repaired by changing its
matching or cycle: its boundary phase was \(-1\), forcing every admissible
score to be \(3\pmod4\), whereas \(8-3p=-151=1\pmod4\).

## Line-product form

Let \(a\) be the pencil centre and \(b=\ell_1\cap\ell_2\).  For both the
ordinary and triple two-outlier types,

\[
 D=D_0\mathbin\triangle\ell_1(a)\mathbin\triangle\ell_1(b)
     \mathbin\triangle\ell_2(a)\mathbin\triangle\ell_2(b).
\]

Writing \(P(\ell)=\prod_{v\in\ell}x_v\), symmetric-difference cancellation
gives

\[
 \prod_Dx=\left(\prod_{D_0}x\right)
 P(\ell_1(a))P(\ell_1(b))P(\ell_2(a))P(\ell_2(b)).
\tag{2}
\]

If \(\mathcal O\) is the set of opposite directions, then \(D_0\) is the
disjoint union of the punctured opposite pencil rays, so

\[
 \prod_{D_0}x=x_a\prod_{L\in\mathcal O}P(\ell_L(a)).
\tag{3}
\]

Equations (2)--(3) apply to every all-bad three-spike completion.  The only
geometric difference between the two near-pencil types is whether the
connector direction \(\langle b-a\rangle\) occurs in the opposite pencil
product in (3).  They also apply after an arbitrary relative affine or
signed-projective alignment: first express the transported Boolean
completion in the affine chart and gauge in which the near-pencil quotas are
stated, and then use its resulting line products in (2)--(3).  This is an
identity in that gauge, not a claim that the boundary product itself is
invariant under diagonal switching.

## What controls the phase

The boundary phase is not determined by the three spike incidences or by the
unique-circle mismatch count.  Exact \(p=53\) witnesses use the same
Kiss--Somlai completion and change only its translation relative to a fixed
near pencil:

| type | translation with phase \(+1\) | translation with phase \(-1\) | finite spikes in \(D\) |
|---|---:|---:|---:|
| ordinary | \((0,0)\) | \((0,1)\) | both / both |
| triple | \((0,0)\) | \((3,0)\) | both / both |

In every row both finite spike coordinates belong to \(D\), yet the two
translations give opposite products in (1).  Translation also carries the
unique oriented circle and its mismatch set along with the completion, so
the mismatch number is unchanged.  Thus (1) is an independent global
line-product bit of the relative alignment.

For this particular Kiss--Somlai family, the common-Max test is now even
stronger.  The same construction has circle mismatch \(\mu=0\): subtracting
the augmented hard line leaves a Boolean \(+p\) eigenvector at distance
\(p-2\).  Equations (14)--(22) of the
[two-half near-pencil reduction](NOTE_2026-09-02_TWO_HALF_NEAR_PENCIL_REDUCTION.md)
use the degree-surplus cut bound to exclude that zero-mismatch family in both
two-outlier types throughout \(p\ge53\), regardless of the phase.  The
phase-\(+1\) half remains relevant for other, positive-mismatch completions
of the same signed triple; those completions are not classified or excluded
here.

The exact replay is

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q \
  tests/test_residual_near_pencil_spike_phase.py
```

The implementation is `src/residual_near_pencil_spike_phase.py`.  It checks
both factorizations (2)--(3), both two-outlier types, and the stated
same-incidence/opposite-phase witnesses.  This is a proved necessary gate and
a route correction, not residual closure.
