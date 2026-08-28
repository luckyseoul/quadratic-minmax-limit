# Strong integral lift and the `p=29` next all-finite endpoint

Date: 2026-08-28. This is Proposition 15.681. It proves a stronger
all-prime mass bound for nonnegative integer-valued quadratics and uses it,
together with an explicit finite-geometry classification, to close the
`p=29,s=24` second all-finite endpoint. It also removes every
positive-residue row at the corresponding `p=31,37,41` endpoints. The
`p=31,41` residue-zero rows, `p=17,19,23`, later all-finite sizes, the
infinity-present remainder, residual (ii), R1, global QVAR, Type I, and the
limit remain open.

## 1. The paired-cube bound applies before Booleanization

Let `p=2m-1` and let `B` be a nonzero, nonnegative, integer-valued
quadratic on `J(p,m)`. Choose `X` with `B(X)=h>=1`. Leave one uniformly
chosen point of `X` fixed, pair the other `m-1` points bijectively with the
complement, and choose one endpoint of each pair. Averaging first over the
resulting Boolean cube and then over all pairings defines a Markov operator
`T`. With `rho=1/(p+1)`, direct calculation gives

```text
T(1)       = 1,
T(x_i)     = 1/2 + rho*x_i,
T(x_i*x_j) = 1/4 + rho*x_i*x_j,

T(B) = rho*B + (1-rho)*E[B].                         (1)
```

Every paired cube contains `X`, so the restriction of `B` is nonzero and
has degree at most two. The elementary cube-distance bound says that a
nonzero degree-two polynomial is supported on at least one quarter of a
Boolean cube. Since every nonzero value of `B` is a positive integer, the
average of `B` on every cube is at least `1/4`. Put `c=4p E[B]`. Evaluating
(1) at `X` gives

```text
c >= p+1-4h.                                         (2)
```

Proposition 15.642's stabilizer certificate, applied at the same point,
gives

```text
p=4r+3: c >= 4h,
p=4r+1: c >= 4r*h/(r+1).                             (3)
```

For `p=4r+3`, averaging the two inequalities in (2)--(3) gives
`c>=(p+1)/2`. For `p=4r+1`, take weight `(r+1)/(2r+1)` on (3) and the
complementary weight on (2); the `h` coefficient cancels and gives
`c>=(p-1)/2`. Thus

```text
4p E[B] >= (p+1)/2  when p=3 mod 4,
4p E[B] >= (p-1)/2  when p=1 mod 4.                  (4)
```

The scaled floors at `p=29,31,37,41` are respectively

```text
14, 16, 18, 20.                                      (5)
```

This improves the earlier degree-two support floor without assuming that
`B` is Boolean.

## 2. Exact endpoint residues

At the second even all-finite size above `3(p-1)/4`, exact quotient/floor
dynamic programming gives the following phase-zero pair survivors:

```text
p   s    surviving u_0       largest positive scaled mean 2u_0
29  24   0,2,3,4,5                         10
31  26   0,2,3,4,5,6                       12
37  30   2,3,4,5                           10
41  34   0,2,3,4,5,6,7                     14
```

Every positive residue has quotient sum `m-u_0<m`, hence forces a
quotient-zero direction. Its scaled mean `2u_0` is below every positive
phase-zero fibre floor, so its odd-fibre count is zero and its pointwise
slack is `A=2B` for a nonzero `B` of the kind covered by (4). Comparison
with (5) excludes every positive row. The already closed `p=37` endpoint
gets a shorter independent proof. Only residue zero survives at `p=29,31`
and `41`.

## 3. The five exact `p=29` residue-zero profiles

Here `m=15`, `P=30`, `s=24`, and the pair-deficit budget is
`s(s-1)=552`. Phase zero has minimum profile

```text
10*b=0 + 5*b=24,   deficit 240,
```

while phase one has

```text
14*b=2 + 1*b=24,   deficit 308.
```

If a line has occupancy `n`, its contribution to the slack between the
actual pair count and the parity pair floor is

```text
2*(C(n,2)-floor(n/2)) =
    4r(r-1), n=2r,
    4r^2,    n=2r+1.                                 (6)
```

Global pair slack is therefore a nonnegative multiple of four. Enumerating
every exact type profile within the four available deficit units leaves
five phase-labelled profiles but only three global secant distributions.
Writing `t=(24-b)/2`, they are

```text
pair slack 4: {t=12:10, t=11:14, t=0:6},
pair slack 0: {t=12:10, t=11:14, t=2:1, t=0:5},
pair slack 0: {t=12:10, t=11:14, t=1:2, t=0:4}.      (7)
```

By (6), the first shape has exactly one line of occupancy three and every
other line has occupancy at most two. Its six `t=0` directions are
undetermined. Deleting one point from the unique collinear triple leaves a
23-arc and preserves those six directions. The other shapes are 24-arcs
with at least four undetermined directions.

## 4. The finite-classification exit

Coolsaet and Sticker's exhaustive classification of all arcs, not merely
complete arcs, gives the following numbers of projective classes in
`PG(2,29)` (Table 5):

```text
25-arcs: 10 classes,
26-arcs:  5 classes.                                 (8)
```

The source independently generates `PGL(2,29)` from translation,
inversion, and primitive scaling. Its order is exactly
`29*(29^2-1)=24360`. Exact orbit enumeration on subsets of the 30-point
projective line gives

```text
5-subset complements: 10 orbits among C(30,5)=142506,
4-subset complements:  5 orbits among C(30,4)=27405. (9)
```

Every complement gives a 25- or 26-point subset of a nondegenerate conic.
Five arc points determine that conic uniquely, so two such subsets are
projectively equivalent exactly when their complements are in the same
`PGL(2,29)` orbit. The lower bounds (9) therefore exhaust the total class
counts (8): every 25- and every 26-arc in `PG(2,29)` lies on a conic. This
also agrees with Chao--Kaneta's earlier result that the largest
nonclassical arc in `PG(2,29)` has size 24.

For either 24-arc shape in (7), adjoin any two undetermined infinity points
to get a 26-arc. For the one-triple shape, first delete one triple point and
then adjoin any two undetermined infinity points to get a 25-arc. Choose
three distinct undetermined points `D_1,D_2,D_3`. The conics containing the
extensions by `{D_1,D_2}` and `{D_1,D_3}` share at least 23 affine arc
points, so they coincide. That conic would contain the three collinear
points `D_1,D_2,D_3` on the line at infinity, impossible for a
nondegenerate conic. Thus every residue-zero profile is excluded and the
`p=29,s=24` endpoint is closed.

## 5. Literature, OEIS, and reproduction

The `PG(2,29)` arc classification is an explicit external dependency. The
geometry itself is not claimed as new: Chao--Kaneta already recorded
maximum nonclassical size 24, and Coolsaet--Sticker later supplied the full
class counts. The new use is the combination with the Paley phase ledger
and the integral paired-cube lift (4). Targeted searches found no source
combining those ingredients or stating (4) in this form.

OEIS searches found `27405` and `142506` in binomial-coefficient tables, as
expected from `C(30,4)` and `C(30,5)`. Searches for the distinctive endpoint
deficit block `240,214,216,238,240,308`, the arc-class tail
`43,10,5,1,1,1`, and the PGL orbit-size block found no relevant sequence.
No sequence or priority claim is made.

Reproduce with

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15681.py
PYTHONPATH=src pytest -q tests/test_prop15642.py tests/test_prop15669.py \
  tests/test_prop15675.py tests/test_prop15679.py tests/test_prop15680.py \
  tests/test_prop15681.py
```

The generated exact record is `evidence/e1_gmin_m4_prop15681.json`.
