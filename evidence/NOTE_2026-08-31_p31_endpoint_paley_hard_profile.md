# The `p=31` endpoint has a Paley-hard near-pairing type

**Date:** 2026-08-31
**Status:** proved necessary normal form, not endpoint exclusion.

## Exact hypotheses

Assume the Proposition 15.632 residual affine inequalities for an odd edge
set `H` with `p=31` and `|H|=4p+1=125`.  In an outside chart, assume its
odd-degree boundary `D` is all finite, has `|D|=p+1=32`, and has the first
unexcluded pair slack `R=10`.  Proposition 15.727 then applies: if `x,y`
count the disjoint trisecants and 4-secants, then

`x+2y=10`.

The conclusion below is not a statement about an arbitrary 32-point set.
It uses both exact quadratic-type slack budgets of the residual separator.

## Odd-fibre and phase identities

For a direction `d`, let `b_d` be the number of its 31 affine fibres that
meet `D` oddly.  Since `|D|=32`, every `b_d` is even.  Proposition 15.722's
pair identity gives

`sum_d b_d = 32+4R = 72`.                                      (1)

The parity sign in Proposition 15.632 is

`eps_d*(-1)^((|H|-3)/2)*c_H*(-1)^b_d`.

Here the exponent is 61 and `b_d` is even, so this is `-eps_d*c_H`.
Consequently all sixteen directions with `eps_d=c_H` have phase one, while
the other sixteen have phase zero.  Each type has exact scaled-mean budget

`(p+1)^2/2=512`.                                               (2)

## The phase-one type is rigid

For even `b`, the exact `p=31` phase-one parity floors are

`floor(b)=30` for `b=2,30`, and `floor(b)=62` otherwise.       (3)

Same-type scaled means have one residue modulo 32.  Write

`a_d=2u+32k_d`, `0<=u<16`.

Summing over the type and using (2) gives

`sum_d k_d=16-u`.                                             (4)

If `1<=u<=14`, even a `b=2` or `30` direction needs `k_d>=1`, so
the sixteen directions already require more quotient mass than (4)
provides.  If `u=0`, the only possible quotient-one cells are `b=2,30`,
where the mean would be 32, exactly two above the floor.  Their pointwise
parity baselines are respectively `q_0=(1-x_i-x_j)^2` and, after
complementing the 30-set, `q_0=1-x_i`; both have scaled mean 30.  A
putative slack `A` of scaled mean 32 would make `C=(A-q_0)/2` a nonzero
nonnegative integral quadratic with `4p E[C]=2`, contradicting Proposition
15.688's sharp floor `4p E[C]>=p-3=28`.

Therefore `u=15`.  Equation (4) now has quotient sum one, and (3) forces
the exact mean multiset

`15 copies of 30, one copy of 62`.                            (5)

The opposite type has common residue `2u_0`.  Adding the two instances of
the exact mean formula

`a_d=I+32P_d-eps_d*T-93`

cancels the finite signed sum `T` and gives

`30+2u_0 = 2I+6 (mod 32)`.

Infinity is outside the boundary and so has even degree `I`.  Consequently

`u_0 in {0,2,4,6,8,10,12,14}`.                              (5a)

The fifteen mean-30 directions all have `b in {2,30}`.  Two `b=30`
directions, together with the other thirteen required special directions at
`b=2`, would already contribute

`2*30+13*2=86>72`

to (1).  Hence there is at most one `b=30` direction in this type and

`at least fourteen directions with b=2`.                     (6)

## Interaction with the disjoint rich blocks

Fix a `b=2` direction, and let `r_3,r_4` count its trisecants and
4-secants.  Writing `l_j` for its number of `j`-point fibres, the point and
odd-fibre identities give

`l_0=14+r_3+r_4`, `l_1=2-r_3`,

`l_2=15-r_3-2r_4`, `l_3=r_3`, `l_4=r_4`.                    (7)

In particular, a non-rich `b=2` direction has profile

`(l_0,l_1,l_2,l_3,l_4)=(14,2,15,0,0)`.                       (8)

At block parameter `y`, Proposition 15.727 has `x=10-2y` and only
`x+y=10-y` rich lines.  They occupy at most that many directions.  Combining
this with (6), one Paley type therefore contains at least

`14-(10-y)=4+y`                                              (9)

non-rich directions with profile (8).

## What this advances

Proposition 15.727 left the `p=31,R=10` endpoint as an arbitrary disjoint
3/4-secant normal form.  The new necessary condition says that, in addition,
one Paley half of the direction set is almost entirely made of two-odd-fibre
projections and contains at least `4+y` near-perfect pairing directions.
This is a strict symbolic reduction and uses no arc classification or
configuration solver.

It does **not** exclude the endpoint.  The next geometric implication to
prove is that a 32-point endpoint set with the disjoint block structure
cannot have the same-Paley near-pairing directions forced by (8)--(9), or is
forced into one of the already-closed circle/conic configurations.

## Artifacts

- `src/e1_gmin_m4_prop15728.py`
- `tests/test_prop15728.py`
- `evidence/e1_gmin_m4_prop15728.json`
