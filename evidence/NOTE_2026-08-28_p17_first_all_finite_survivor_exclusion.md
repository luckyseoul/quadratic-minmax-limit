# Exceptional `p=17` first all-finite survivor excluded

Date: 2026-08-28. This is Proposition 15.678. It closes the `p=17`,
`s=14` endpoint deliberately left by Proposition 15.677. Together with
15.675 and 15.677, the first even all-finite boundary size strictly above
`3(p-1)/4` is now excluded for every prime `p>=17`. Later all-finite sizes,
strict pair deficit in infinity-plus-`p`, residual (ii), R1, global QVAR,
Type I, and the limit remain open.

The final step imports Sticker's exhaustive finite classification of arcs in
`PG(2,17)`. That dependency is explicit below; it is not presented as a
locally reproduced classification.

## 1. Exact residue ledger

Here

```text
p=17,  s=14,  m=9,  P=p+1=18,
pair-deficit budget = s(s-1)=182.
```

For one quadratic direction type, the exact means have the form

```text
a_d=2u+18k_d,     sum_d k_d=9-u.                    (1)
```

Enumerating every even `0<=b_d<=14` against the exact symbolic floor, while
excluding a nonzero lift of size two as in Proposition 15.642, gives

```text
phase zero u :   0   2   3   4   5    6    7    8
minimum deficit: 84  82  84  96  98  110  112  112

phase one u :    8 only
minimum deficit: 96.
```

The source regenerates this ledger directly; it is not copied from a solver
output.

### The `u_0=2` exit

Equation (1) has quotient sum seven across nine directions, so some
direction has `k_d=0` and mean four. Every nonzero even `b_d` has phase-zero
floor at least 18, hence this direction has `b_d=0`. Its pointwise slack is
`A_d=2B_d`, where `B_d` is nonzero, nonnegative, integer-valued, and
quadratic. At `p=17`, Proposition 15.642 gives exact nonzero scaled cost at
least six, contradicting its scaled mean four.

### The `u_0=3` exit

The coefficient congruence from Proposition 15.677 has `q=8`,
`1<=j<=7`, and `q | u_0-j`, so `j=3`. Its three surviving coefficient rows
are

```text
ell   I   E   xnor l1 lower bound   transverse capacity
 0   44  25           65                    25
 2   26  43           63                    41
 4    8  61           63                    57.
```

All three are impossible. Residues `u_0>=4` have total minimum deficit at
least `96+96=192>182`. Thus only `u_0=0` remains.

## 2. Pair slack leaves exactly two profiles

If an affine line has occupancy `n`, its contribution to pair-budget slack
is

```text
2*(C(n,2)-floor(n/2)).                                (2)
```

For `n=2r`, (2) is `4r(r-1)`; for `n=2r+1`, it is `4r^2`. Therefore global
pair slack is a nonnegative multiple of four. It vanishes exactly when every
line has occupancy at most two.

At `u_0=0`, phase zero can have deficit 84 or 86 before phase one's minimum
would exhaust the budget. Phase one can have deficit 96 or 98. The nominal
combination `84+96=180` has slack two and is impossible by (2). The two
remaining combinations are

```text
case A:
  phase zero b-profile {0:6, 14:3},       deficit 84
  phase one  b-profile {2:8, 12:1},       deficit 98

case B:
  phase zero b-profile {0:6, 12:1, 14:2}, deficit 86
  phase one  b-profile {2:8, 14:1},       deficit 96.
```

Both total 182, so every affine line has occupancy at most two: the fourteen
finite boundary points form an arc. Writing `t_d=(14-b_d)/2` for the number
of secants in direction `d`, both cases have the same global distribution

```text
t=7: 6 directions
t=6: 8 directions
t=1: 1 direction
t=0: 3 directions.                                  (3)
```

In particular there are exactly three undetermined directions, represented
by three collinear points on the line at infinity.

## 3. The 16-arc classification exit

Heide Sticker's exhaustive classification lists the number of
PGL-inequivalent `k`-arcs in `PG(2,q)`, without requiring completeness. In
the `q=17` column the counts for `k=14,15,16,17,18` are

```text
4, 1, 1, 1, 1.
```

See [Section 5.3, printed page 119 of the
thesis](https://cage.ugent.be/geometry/Theses/57/PhDHeideSticker.pdf). The
classification uses isomorph-free canonical augmentation; the thesis also
reports an orbit-stabilizer/double-count consistency check for every
`(k,2)`-arc computation with `q<=27`.

A nondegenerate conic in `PG(2,17)` is an 18-arc, so deleting two points
produces a 16-arc. Since the classification has exactly one 16-arc class,
every 16-arc is projectively equivalent to this representative and is
contained in a conic.

Let `D_1,D_2,D_3` be the three undetermined infinity points from (3). Adding
any two, say `D_1,D_2`, to the affine 14-arc `S` gives a 16-arc: neither
point lies on an `S`-secant, and their common line at infinity contains no
point of `S`. Hence

```text
S union {D_1,D_2} subset C
```

for some conic `C`. The line at infinity already meets `C` at its two
distinct points `D_1,D_2`, so `D_3` is off `C`. Also `|C\S|=18-14=4`.

An external off-conic point in odd order lies on `(q-1)/2=8` conic secants;
an internal point lies on `(q+1)/2=9`. Deleting four conic points destroys
at most four secants through a fixed point. Thus `D_3` lies on at least four
secants of `S` if external and at least five if internal. This contradicts
that `D_3` is undetermined. Both profiles A and B are impossible.

## 4. Literature, OEIS, and reproduction

The finite classification is the only new imported ingredient. Keri's
peer-reviewed MDS/superregular-matrix classification is an earlier
independent classification in the same range; Sticker reports agreement
with it. The conic secant counts are standard finite-conic incidence counts.

An OEIS search for the distinctive larger block of the `q=17` class-count
row, `17633,21064,6814,629`, returned no result. This is only a context and
duplicate check; Proposition 15.678 is not a sequence claim.

Reproduce the exact local arithmetic record with

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15678.py
PYTHONPATH=src pytest -q tests/test_prop15642.py tests/test_prop15669.py \
  tests/test_prop15675.py tests/test_prop15677.py tests/test_prop15678.py
```

The generated record is `evidence/e1_gmin_m4_prop15678.json`.
