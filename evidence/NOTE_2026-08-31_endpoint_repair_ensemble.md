# Endpoint equality forces a complementary-arc repair ensemble

**Date:** 2026-08-31
**Status:** proved structural reduction, not endpoint exclusion.

**Later status:** Proposition 15.734 closes every `k=4p` endpoint boundary
for `p>=13`. The repair ensemble remains a valid conditional structure, but
is superseded as a live residual gate.

## Hypothesis inherited from Proposition 15.727

Let `p>=17` be prime, put

`R=floor((p-1)/3)`, `p=3R+c`, `c in {1,2}`,

and suppose endpoint equality survives.  Proposition 15.727 gives an affine
set `D` of `p+1` points and a minimum deletion

`D=A disjoint-union T`, `|A|=k=p+1-R`, `|T|=R`,

such that `A` is an arc.  Every line containing at least three points of `D`
is one of a pairwise `D`-disjoint family of `x` trisecants and `y` 4-secants,
each containing exactly two points of `A`, with

`x+2y=R`.                                                   (1)

Call these lines the rich blocks, and let `S` be the points of `D` outside
them.  Exact point counting gives

`|S|=p+1-(3x+4y)=c+1+2y`.                                 (2)

## Every maximum repair, not just one repair

On every rich block `L`, choose a two-point subset `Q_L`, and put

`A_Q=S union (union_L Q_L)`, `T_Q=D\A_Q`.                  (3)

Any trisecant of `A_Q` would be a rich line of `D`, hence one of the blocks;
but (3) retains exactly two points on every block.  Thus `A_Q` is an arc.
The same argument applies to `T_Q`: it retains one point on a trisecant and
two on a 4-secant, so it is also an arc.  Their sizes are

`|A_Q|=2(x+y)+|S|=p+1-R=k`, `|T_Q|=x+2y=R`.               (4)

Conversely, any `k`-point arc contained in `D` must delete at least one point
from every trisecant and at least two from every 4-secant.  These disjoint
demands already total `x+2y=R`; equality forces it to retain all of `S` and
exactly two points on every block.  Therefore (3) lists **all** maximum
`D`-subarcs, and their number is exactly

`3^x 6^y`.                                                 (5)

For every `z in T_Q`, its own block contains two points of `A_Q`, so it gives
an `A_Q`-secant through `z`.  A second such secant would be another rich
block through `z`, contradicting block disjointness.  Hence every repair has
an `R`-point complementary arc consisting entirely of index-one points.  The
unique secants form a matching on `A_Q`; their fibres have size one on each
trisecant and size two on each 4-secant.

In particular, all `R` sets `A_Q union {z}`, `z in T_Q`, are affine
`(k+1,3)`-arcs with exactly one trisecant.  This simultaneous family is
strictly more information than selecting one set in Proposition 15.729.

## Exact two-colour projective line census

For a fixed repair write

`n_ij = #{projective lines L: |L intersect A_Q|=i,
                              |L intersect T_Q|=j}`.

Both colours are arcs, and a line of type `(1,2)` would be a forbidden rich
line not having two `A_Q` points.  The complete census is

| `i\j` | `0` | `1` | `2` |
|---:|---:|---:|---:|
| `0` | `p(p-1)/2-R-y` | `2R+2y` | `binom(R,2)-y` |
| `1` | `k+2R` | `R(k-2)` | `0` |
| `2` | `binom(k,2)-R+y` | `x=R-2y` | `y` |

The `(2,1)` and `(2,2)` cells are the rich blocks.  Counting `A_Q` pairs,
`T_Q` pairs, and cross-colour pairs gives

`n_20+n_21+n_22=binom(k,2)`,

`n_02+n_22=binom(R,2)`,

`n_11+2n_21+4n_22=kR`.                                  (6)

An `A_Q` point has `R+1` tangents and a `T_Q` point has `k+1` tangents, so

`n_10+n_11=k(R+1)`,

`n_01+n_11+n_21=R(k+1)`.                                (7)

Equations (1), (6), and (7), followed by the total number
`p^2+p+1` of projective lines, give the displayed table.  Aggregating by
`i+j` recovers the Proposition 15.727 occupancy census exactly.

## Point signatures and co-tangent multiplicity

The local line signatures are also forced.  Superscripts below are
multiplicities of lines through the indicated point.

- Each of the `c+1+2y` singleton `A_Q` points has
  `(2,0)^(k-1), (1,1)^R, (1,0)^1`.
- Each of the `2x` `A_Q` points on a trisecant block has
  `(2,1)^1, (2,0)^(k-2), (1,1)^(R-1), (1,0)^2`.
- Each of the `2y` `A_Q` points on a 4-secant block has
  `(2,2)^1, (2,0)^(k-2), (1,1)^(R-2), (1,0)^3`.
- Each of the `x` `T_Q` points on a trisecant block has
  `(2,1)^1, (1,1)^(k-2), (0,2)^(R-1), (0,1)^2`.
- Each of the `2y` `T_Q` points on a 4-secant block has
  `(2,2)^1, (1,1)^(k-2), (0,2)^(R-2), (0,1)^3`.

Let a rich block contain the `A_Q` pair `{a,b}`.  The arc
`A_Q\{a}` has all points in

`{a} union (T_Q intersect L)`

as individually valid extension points on the tangent `L` through `b`, and
no two of them can be adjoined together.  Repeating with `a,b` exchanged
shows that a fixed repair has

- `2x` size-`p-R` arc bases with two pairwise-incompatible co-tangent
  extensions; and
- `2y` size-`p-R` arc bases with **three** pairwise-incompatible co-tangent
  extensions.

The second clause strengthens the two-extension consequence of 15.729
whenever a 4-secant is present.

## Direction-level refinement

Fix an affine direction `d`.  Let `sigma_d` and `tau_d` be the numbers of
`A_Q`- and `T_Q`-secants in that direction, let `r_3(d),r_4(d)` count its
rich trisecants and 4-secants, and let `m_d` count its `(1,1)` lines.  Among
the `p` affine lines of that direction the two-colour cells are

`n_20=sigma_d-r_3-r_4`, `n_21=r_3`, `n_22=r_4`,

`n_10=k-2sigma_d-m_d`, `n_11=m_d`, `n_12=0`,

`n_02=tau_d-r_4`, `n_01=R-2tau_d-r_3-m_d`,

`n_00=sigma_d+tau_d+r_3+r_4+m_d-1`.                      (8)

All nine quantities are nonnegative.  If `b_d` is the number of odd
`D`-fibres in direction `d`, then

`b_d=p+1-2(sigma_d+tau_d+m_d)`.                           (9)

Consequently the ordinary `D`-fibre profile is

`l_4=r_4`, `l_3=r_3`,

`l_2=(p+1-b_d)/2-r_3-2r_4`, `l_1=b_d-r_3`,

`l_0=(p-1-b_d)/2+r_3+r_4`.                               (10)

Summing (9) over the `p+1` directions recovers

`sum_d b_d=p+1+4R`.                                      (11)

Unlike a one-object arc classification, (8)--(10) retain the directional
data needed by the residual separator.

## Historical next implication

This proposition does not exclude endpoint equality.  It replaces the
incorrect Bartoli--Storme ceiling route recorded after 15.729 with a stronger
internally proved target: a surviving endpoint would support an entire
product family of complementary arcs, exact local swaps, and the census
(8) in every direction.

Each `A_Q` is one point (`c=2`) or two points (`c=1`) below the sufficient
size hypothesis in Ball--Lavrauw's explicit tangent interpolation formula.
Proposition 15.731 shows directly that the envelopes nevertheless exist:
they are unique in the first residue and form line-product pencils in the
second. One-block swaps then have quadratic or cubic transition data. The
exact open implication is exclusion of the common `D` under (8) and the
phase-specific lift constraints. Proposition 15.732 proves that the bare
linear cycle circulation of the transitions is exact and identifies a
nonzero rich-direction first jet, but no phase bridge. Proposition 15.734
later bypasses this implication with an isolated-chart coefficient
contradiction. No universal size
ceiling for odd-order unique-trisecant 3-arcs is available from [70].

## Artifacts

- `src/e1_gmin_m4_prop15730.py`
- `tests/test_prop15730.py`
- `evidence/e1_gmin_m4_prop15730.json`
