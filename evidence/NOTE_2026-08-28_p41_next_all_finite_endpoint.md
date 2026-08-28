# The `p=41` next all-finite endpoint is impossible

Date: 2026-08-28. This is Proposition 15.683. It closes the `p=41,s=34`
second all-finite endpoint by combining Proposition 15.681's integral lift,
an exact nine-profile residue-zero ledger, and Segre's degree-`2t` tangent
envelope in the polynomial form of Ball--Lavrauw. The same boundary remains
open only at `p=17,19,23`. Later all-finite sizes, the infinity-present
remainder, residual (ii), R1, global QVAR, Type I, and the limit remain open.

## 1. Exact residue-zero ledger

At `p=41`, put

```text
m=21, P=42, s=34, pair-deficit budget=s(s-1)=1122.
```

The exact quotient/floor ledger has phase-one residue `u_1=20` and
phase-zero pair survivors

```text
u_0 = 0,2,3,4,5,6,7.
```

Every positive residue forces a quotient-zero `b=0` direction with scaled
quadratic-lift mass at most 14. Proposition 15.681 gives the stronger floor

```text
4p E[B] >= (p-1)/2 = 20,
```

so only `u_0=0` remains. Its type minima are

```text
phase zero: 14*b=0 + 7*b=34, deficit 476,
phase one:  20*b=2 + 1*b=34, deficit 640.
```

Exact completion-bounded enumeration gives seven phase-zero rows, four
phase-one rows, and nine compatible phase-labelled pairs. As before, a line
of occupancy `n` contributes

```text
2*(C(n,2)-floor(n/2))
```

to pair slack. Thus slack zero is a 34-arc, while slack four means exactly
one 3-secant and every other line has occupancy at most two. Writing
`t=(34-b)/2` for the number of floor-secants in a direction, the nine rows
have only four global shapes:

```text
slack 4 (two phase-labelled rows):
  {t=17:14, t=16:20, t=1:1,                  t=0:7};

slack 0 (seven phase-labelled rows):
  {t=17:14, t=16:20, t=3:1,                  t=0:7},
  {t=17:14, t=16:20, t=2:1, t=1:1,           t=0:6},
  {t=17:14, t=16:20,          t=1:3,          t=0:5}.
```

The fourteen `t=17` directions are perfect matchings and the twenty `t=16`
directions are matchings missing two points. Call the other eight
directions *exceptional*.

## 2. The tangent-envelope component lemma

For an arc `A` of size `q+2-tau` in odd order, Ball--Lavrauw Theorem 11
applies when `|A|>2*tau+2`. It supplies a nonzero homogeneous polynomial
`Phi` of degree `2*tau` in the dual plane. If `P` is an arc point and `P*` is
its dual pencil line, then

```text
Phi restricted to P* = (the tangent polynomial at P)^2.
```

Consequently every tangent dual point is a double zero on its point-pencil
line.

Let `D` be a point at infinity and `D*` the dual line parametrising all
affine lines of direction `D`. If that direction contains `b>2*tau` tangents,
then `Phi|D*` has more roots than its degree, so `D*` divides `Phi`. Write
`Phi=D* Psi`. At each of those tangent points, the corresponding affine
point-pencil meets `D*` transversely. Its required double zero therefore
forces `Psi` to vanish there. If also `b>2*tau-1`, then `D*` divides `Psi` and

```text
(D*)^2 divides Phi.                                   (1)
```

After such direction squares are removed, every tangent in another
direction remains a double zero of the residual polynomial on its
point-pencil. If a degree-`d` residual has `r` such tangents through `P`
with `2r>d`, then its restriction to `P*` vanishes identically, so `P*` is
a line component.

## 3. The seven arc profiles

For a slack-zero row, the eight exceptional directions contain three
secant edges in total. Hence each has at most three secants and at least

```text
34-2*3 = 28
```

tangents. A 34-arc in `PG(2,41)` has `t=41+2-34=9` tangents through each
point, and `34>2*9+2`, so its envelope has degree 18. Since `28>18` and
`28>17`, equation (1) makes all eight exceptional direction lines double
components. Their total degree is 16, leaving a residual conic `Q`.

If an arc point `P` is incident with `e_P` of the three exceptional secant
edges, then exactly `8-e_P` exceptional directions are tangent at `P`.
Since `P` has nine tangents in all, it has

```text
9-(8-e_P) = 1+e_P
```

tangents in the remaining directions. The three exceptional edges are
distinct edges of a simple graph, so they touch at least three distinct arc
points. At each touched point `e_P>=1`, and `Q|P*` has at least two distinct
double zeros. Their total multiplicity is at least four, greater than the
degree two of `Q`; therefore `P*` divides `Q`. This would put at least three
distinct line components in a conic, impossible. All seven arc rows are
excluded.

## 4. The two one-triple profiles

The slack-four shape has seven exceptional directions with no floor-secant
and one with one floor-secant. The unique 3-secant need not be the latter:
it could lie in one of the `t=16` directions. Choose a point of the triple
whose deletion preserves the exceptional floor-secant as an ordinary pair.
If the floor-secant is the triple, every deletion does this. Otherwise its
ordinary pair meets the triple in at most one point, so a suitable triple
point exists.

After deletion there is a 33-arc. Seven exceptional directions contain 33
tangents each; the direction of the surviving exceptional pair contains 31.
This arc has `t=41+2-33=10` tangents through each point and satisfies
`33>2*10+2`, so its envelope has degree 20. Since `31>20` and `31>19`, the
eight exceptional direction squares again contribute degree 16, leaving a
quartic `Q`.

The two endpoints of the surviving exceptional pair see seven exceptional
tangents and hence three tangents in the remaining directions. Their dual
point-pencil lines each carry three distinct double zeros of `Q`, total
multiplicity six greater than four, so both lines divide `Q`. Divide them
out and call the residual conic `R`.

Every other one of the 31 arc points sees all eight exceptional directions
as tangents and therefore has two remaining tangents. Neither tangent dual
point lies on either removed point-pencil line: such incidence would mean
that the tangent line contains a second arc point. Thus both are still
double zeros of `R` on the current point-pencil. Two distinct double zeros
force that line to divide `R`. Taking any three of the 31 points would put
three distinct line components in a conic, impossible. Both near-arc rows
are excluded.

## 5. Literature, OEIS, and reproduction

The external theorem is:

- S. Ball and M. Lavrauw, *Planar arcs*, J. Combin. Theory Ser. A **160**
  (2018), 261--287, Theorem 11, doi:10.1016/j.jcta.2018.06.015;
  arXiv:1705.10940v4.

Hirschfeld's survey restates the same classical Segre envelope as Theorem
5.2, including the twice-counted tangent intersections. Searches for an
existing application combining eight high-tangent direction pencils with
the `p=41` profiles found none. OEIS searches for the exact profile blocks
`14,7,20,1`, `18,16,2`, and `34,41,561` returned unrelated sequences; no
sequence claim is made.

Reproduce with

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15683.py
PYTHONPATH=src pytest -q tests/test_prop15681.py tests/test_prop15682.py tests/test_prop15683.py
```

The generated record is `evidence/e1_gmin_m4_prop15683.json`.
