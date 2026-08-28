# Pair-deficit equality close in the infinity-plus-`p` shell

Date: 2026-08-28. This is Proposition 15.676. Let a residual candidate have
odd-degree boundary consisting of infinity and `p` finite points. If its
directional pair-deficit inequality is an equality, the finite boundary is
a `p`-arc. The proposition excludes that complete equality branch for both
product signs and every prime `p>=17`. The strict pair-deficit branch, hence
the full infinity-plus-`p` shell, remains open.

## 1. Equality gives exactly two affine conic profiles

For a direction `d`, let `b_d` be the number of fibres containing an odd
number of the `p` finite boundary points. Fibre by fibre,

```text
p-b_d <= 2 * (number of colliding unordered pairs in direction d).
```

After summing over all directions,

```text
sum_d (p-b_d) <= p(p-1).
```

Equality holds exactly when no fibre contains three points. Thus the finite
set is a `p`-arc in `PG(2,p)`. Segre's odd-order `p`-arc theorem puts it on
a nondegenerate conic.

A conic has `p+1` projective points. If the line at infinity is secant, only
`p-1` are affine, too few. Two cases remain.

- If infinity is tangent, all `p` affine conic points occur. The tangent
  direction has one point in every affine fibre, so `b=p`; every other
  direction has one tangent singleton and all other occupied fibres are
  pairs, so `b=1`. The profile is `p*b=1 + one b=p`.

- If infinity is external, all `p+1` conic points are affine and the
  boundary omits one point. The external infinity line has `m=(p+1)/2`
  internal and `m` external points relative to the conic. Internal
  directions give `b=1`. Of the external directions, the tangent at the
  omitted point also gives `b=1`; the other `m-1` directions give `b=3`.
  The profile is `(m+1)*b=1 + (m-1)*b=3`.

The equal split on an external line is also immediate in the normal form
`x^2-nu*y^2=z^2`, with `nu` nonsquare: the affine character sum
`sum_r chi(r^2-nu)=-1`, while the projective infinity direction contributes
`+1`.

## 2. The external-conic profile exceeds the type budgets

Put `P=p+1=2m`. The exact floors at `b=1,3` are

| phase | `b=1` | `b=3` |
|---:|---:|---:|
| zero | `P` | `2P-8` |
| one | `P-2` | `2P-2` |

If one quadratic direction type contains `t` of the `b=3` directions, its
floor sum is

```text
phase zero: mP+t(P-8),
phase one:  mP+(t-1)P.
```

The type budget is exactly `mP`. Phase zero permits no `b=3` direction;
phase one permits at most one per type. Globally the conic profile has
`m-1>=8` such directions, exceeding both two-type capacities. Hence the
external-conic-minus-one-point branch is impossible in both phases.

## 3. The tangent profile fails coefficient arithmetic

In phase zero a baseline `b=1` mean is `P`; in phase one it is `P-2`. The
exact common-residue equation

```text
a_d=2u+P k_d,   sum_d k_d=m-u
```

together with the forbidden two-unit lift leaves a baseline `b=1`
direction in each quadratic type. Relative to its baseline parallel-edge
count, the sum of the parallel-count offsets in one type is zero in phase
zero and one in phase one. Thus, writing `x,y` for the two baseline counts,

```text
E=m(x+y)       in phase zero,
E=m(x+y)+2     in phase one.
```

Explicitly, phase zero forces common residue `u=0`. A type containing the
`b=p` direction can spend its one remaining quotient either there or on one
`b=1` direction, but in either case at least `m-2` `b=1` directions stay at
mean `P` and the aggregate offset is zero. Phase one forces `u=m-1`; each
type has one quotient unit, at least `m-2` `b=1` directions stay at mean
`P-2`, and the aggregate offset is one. The source independently enumerates
all residues and quotient allocations to verify these are exhaustive.

There is one immediate exception to this arithmetic reduction: for
`p=1 mod 4` in phase zero, the type containing the `b=p` direction already
has floor sum

```text
(m-1)P+(2P-2)=mP+(P-2)>mP.
```

It is therefore impossible.

In every remaining row the baseline coefficient identity from Proposition
15.673 is

```text
q | I+P_d-(4+sigma),   q=(p-1)/2,
```

where `sigma=+1` in phase zero and `-1` in phase one. Substituting
`I=4p+1-E` in either phase gives

```text
q|x,  q|y.
```

The edge count gives `x+y<=7`, while `q>=8`, hence `x=y=0`. This leaves
`(E,I)=(0,4p+1)` in phase zero or `(2,4p-1)` in phase one. Both violate the
elementary boundary-support inequality

```text
I <= p+2E.
```

Thus the tangent-conic branch is impossible as well.

## 4. Checks, literature, OEIS, and scope

`src/e1_gmin_m4_prop15676.py` constructs canonical tangent and external
conics for sample primes, computes every directional fibre parity, and
reproduces the two profiles and exact pair-deficit equality. It then checks
both phase ledgers through `p=101`. The tests independently assert all four
residue/phase outcomes.

The conic-containment input is Segre's theorem as restated in Ball--Lavrauw,
[Planar arcs](https://arxiv.org/abs/1705.10940). Ball--Csajbók,
[On sets of points with few odd secants](https://arxiv.org/abs/1711.10876),
and Van de Voorde,
[On sets without tangents and exterior sets of a conic](https://arxiv.org/abs/1201.0484),
are adjacent odd-secant and external-conic literature. Targeted searches
found no source combining these conic profiles with the Paley phase floors,
same-type mean sums, and coefficient congruences. OEIS searches for sample
profile-count lists returned unrelated arrays; these are elementary linear
counts, not a sequence claim.

Reproduction:

```bash
python src/e1_gmin_m4_prop15676.py
python -m pytest -q tests/test_prop15676.py
```

The generated record is `evidence/e1_gmin_m4_prop15676.json`. This closes
only pair-deficit equality in the infinity-plus-`p` shell. Strict deficit,
larger all-finite boundaries, general residual (ii), R1, global QVAR,
Type I, and the limit remain open.
