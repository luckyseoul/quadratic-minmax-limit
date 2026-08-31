# Isolated-chart close of the critical residual-(ii) layer

**Date:** 2026-08-31
**Proposition:** 15.734
**Status:** proved theorem

Proposition 15.734 excludes every residual-(ii) witness at `k=4p`, for
every prime `p>=13` and every odd-degree boundary size. It does not exclude
even `k>4p`, and at `p=11` it stops at a sharp lift equality.

## 1. The chart is free

The residual witness supplies an odd flip graph `H` with

```text
|H|=4p+1.
```

At most `2|H|=8p+2` of the `p^2+1` projective vertices are incident with an
edge of `H`. For `p>=13`,

```text
p^2+1-(8p+2)=p^2-8p-1>0.
```

Choose an isolated vertex `w`. It cannot lie in the odd-degree boundary
`D=partial H`. Proposition 15.721's signed PSL transport sends `w` to
infinity while preserving `|H|`, the relative flip set, and both residual
separator inequalities. In the transported chart,

```text
infinity notin D,    I=deg_H(infinity)=0,
```

and `D` is all finite. Handshake makes `|D|` even, so in every affine
direction the number `b_d` of odd boundary fibres is even.

This is the step missed by the boundary-shell campaign: it depends only on
the support size of `H`, not on `|D|`, its pair slack, or its geometry.

## 2. The hard type has only three branches

Put

```text
q=(p-1)/2,    m=q+1.
```

For `A_d=(eps_d S_H-3)/2`, Proposition 15.632 gives each quadratic
direction type the exact scaled-mean budget

```text
sum a_d = 2m^2,    a_d=2p E[A_d].
```

Because `(|H|-3)/2=2p-1` is odd and every `b_d` is even, the directions
with `eps_d=c_H` have phase one. Within that type the exact directional mean
formula has one common residue modulo `p+1=2m`; write

```text
a_d=2u+2m k_d,    0<=u<m,    sum_d k_d=m-u.
```

The phase-one floor is at least `p-1`. Therefore every interior residue
`1<=u<=m-2` would require all `m` quotients `k_d>=1`, although their sum is
strictly below `m`.

At `u=0`, every quotient equals one. A `b=2` cell would be a two-unit lift
of its exact baseline, forbidden by Proposition 15.688 because `2<p-3`.
Thus this residue survives only for `p=1 mod 4`, with all `m` directions in
the exact `b=p-1`, mean-`p+1` baseline.

At `u=m-1`, there are `q` mean-`p-1` directions and one mean-`2p`
direction. For `p=1 mod 4` the low directions all have `b=2`. For
`p=3 mod 4` their possible exact baselines have `b=2` and `b=p-1`; equal
means give one common parallel count, while the two coefficient residues
below differ by one, so they cannot mix.

For `p>=17` these are the exact floors of Proposition 15.669. At `p=13`,
the direct Proposition 15.632 LP gives

```text
phase 0:  b=0,2,4,6,8,10,12 -> 0,14,20,26,24,26,12
phase 1:  b=0,2,4,6,8,10,12 -> 26,12,26,24,26,20,14.
```

The two exceptional middle values are too large to enter the argument, so
the same three branches are exhaustive from `p=13` onward. The LP values
alone are not being used as an equality classification. Proposition
15.652's positive quadrature has strictly positive weights at all three
`b=2` intersection values, so equality forces the pointwise XNOR baseline
`(1-x_i-x_j)^2`. After complementing a `b=p-1` set to its missing fibre,
the two positive `b=1` weights similarly force `x_j` for `p=1 mod 4` and
`1-x_j` for `p=3 mod 4`. This supplies the actual baseline polynomials used
below, including at `p=13`.

## 3. Coefficient residues fix the hard parallel count

Let `P` be the common parallel-edge count in a hard baseline direction and
put `z_s=2x_s-1`, so `sum_s z_s=1`. As in Proposition 15.733, if two
multilinear quadratics agree on this slice, their difference is the
multilinearization of

```text
(sum_s z_s-1)(c+sum_s alpha_s z_s).
```

Constant and linear coefficient comparison, followed by integral
polarization, gives

```text
q | I+P-4                         for target 4+z_i z_j,
q | I+P-(4+sigma)                 for target 4+sigma z_j.
```

The three hard baselines and offsets are therefore

| branch | exact baseline | target | offset |
|---|---|---|---:|
| A | `(1-x_i-x_j)^2` | `4+z_i z_j` | 4 |
| B, `p=1 mod 4` | `x_j` | `4+z_j` | 5 |
| C, `p=3 mod 4` | `1-x_j` | `4-z_j` | 3 |

Write `(P-C)/q=rho` for the appropriate offset `C` and put `s=P+rho`.
Since `I=0`, divisibility and `P-C>-q` give `rho>=0`. Exact finite-edge
accounting gives `q(8-s)` opposite edges in A/B and `q(8-s)+1` in C.
Thus `P<=s<=8`. Since `q>=6`, the three congruences force

```text
(P,rho,s)=(4,0,4), (5,0,5), (3,0,3).
```

## 4. The opposite type is too small

For A/B, an opposite direction with parallel count `Q` has

```text
a=(p-1)s+(p+1)Q+9-7p,
sum Q=q(8-s).
```

For C the corresponding identities are

```text
a=(p-1)(s-7)+(p+1)Q,
sum Q=q(8-s)+1.
```

Nonnegativity forces respectively `Q>=3,2,4`. The total surplus over those
uniform minima is respectively `q-3,q-2,q-3`, always smaller than the
`m=q+1` opposite directions. Hence some opposite direction attains the
minimum and has scaled mean

```text
A: 8,    B: 6,    C: 8.
```

That direction has phase zero. Every nonzero even `b` costs at least
`p-1>=12`, so `b=0`. Its parity is even, hence `A_d=2B_d` for a nonzero
nonnegative integer-valued quadratic `B_d`, and

```text
a_d=4p E[B_d] >= p-3
```

by Proposition 15.688. But `p-3>=10`, larger than each of `8,6,8`. This is
the contradiction.

## 5. Exact scope and next frontier

The theorem closes the complete `k=4p` residual-(ii) layer for every
boundary size and every prime `p>=13`. The same isolated chart exists at
`p=11` and reduces both surviving branches to scaled mass `8`; there
`8=p-3`, so the argument reaches equality rather than contradiction. No
classification of that equality case is asserted here.

Even `k>4p` changes the exact type budget and is not a consequence of this
proposition. The global residual-(ii) predicate, multi-level Type I, and the
quadratic-minmax limit remain open.

## Reproduction

- theorem and exact arithmetic: `src/e1_gmin_m4_prop15734.py`;
- regression tests: `tests/test_prop15734.py`;
- deterministic certificate: `evidence/e1_gmin_m4_prop15734.json`.
