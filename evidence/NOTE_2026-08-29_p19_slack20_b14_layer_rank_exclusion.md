# p=19 slack-20 b=14 layer-rank exclusion

**Date:** 2026-08-29
**Proposition:** 15.695
**Status:** proved; two of four slack-20 profiles excluded, endpoint open

Two slack-20 profiles have phase-one odd-fibre histogram

```text
9*(b=2) + 1*(b=14).
```

Their exact directional floors consume the whole type budget:

```text
9*18 + 38 = 200 = (19+1)^2/2.
```

Thus the `b=14` direction attains its floor. Its symmetrized quadratic is
`q(t)=1`. The exact positive quadrature has support and weights

```text
t=6: 25/76,   t=8: 25/38,   t=10: 1/76.
```

On each of these even intersection layers the phase-one parity lower bound
is one. Equality of the weighted average therefore forces the original
integer-valued slack `A_d(X)` to equal one at every point of all three
layers, not merely after symmetrization.

Pair monomials span all degree-at-most-two functions on `J(19,10)`: linear
and constant terms are recovered from pair sums because the slice weight is
fixed. A fixed 171-row minor of the pair-vs.-10-subset inclusion matrix,
using 91 rows from `t=10`, 75 from `t=8`, and 5 from `t=6`, has exact rank
171 modulo 101. Hence it also has full rank over the rationals. The
quadratic `A_d-1` vanishes on this determining set, so it vanishes on the
entire slice.

But the `t=5` layer is nonempty, and phase-one parity requires `A_d` to be
even there. This contradicts `A_d=1`. Both `b=14` profiles are impossible.
The p=19 remainder falls from seven profiles to five:

```text
{20:2, 24:1, 28:1, 32:1}.
```

The two slack-20 survivors both have five undetermined directions; their
phase-one profiles are `{2:10}` and `{2:9,16:1}`.

## Context check

Targeted literature searches found the standard Gottlieb/Wilson theory of
full subset-inclusion matrices and later rank-resilience work, but not this
restricted three-intersection-layer `171 x 171` certificate. An OEIS search
for the intermediate rank tuple `91,166,171` returned no matching entry.
These searches are context only; no novelty or priority claim is made.

## Reproduction

- `src/e1_gmin_m4_prop15695.py`
- `evidence/e1_gmin_m4_prop15695.json`
- `tests/test_prop15695.py`

The source stores the 171 ten-subset masks and recomputes the rank by exact
finite-field elimination. No solver timeout or floating-point rank is used.
