# First-slack endpoint rigidity and the first four prime closes

**Date:** 2026-08-30
**Status:** proved structural reduction for every prime `p>=17`, plus
classification-assisted exclusion of the endpoint at `p=17,19,23,29`. The
`p+1` shell, residual (ii), Type I, and the limit remain open.

## Result

Let `D` be the `p+1` affine points in an outside chart and put

`R=floor((p-1)/3)`.

This is the first positive slack not excluded by Proposition 15.726. Choose
a minimum-cardinality set `T` such that `A=D\T` is an arc, and write
`t=|T|`. Then:

- `t=R`;
- every `z in T` lies on exactly one `A`-secant;
- every rich line of `D` is a trisecant or 4-secant containing exactly two
  points of `A`;
- the rich lines are pairwise disjoint as subsets of `D`.

If `x` and `y` count the trisecants and 4-secants, respectively, then

`x+2y=R`.

In particular, the arc `A` has at least `R` outside points of secant index
one. In the standard notation, `c_1(A)>=R`.

The archived exhaustive arc classifications then give

| `p` | endpoint `R` | `|A|` | required `c_1(A)` | classified maximum |
|---:|---:|---:|---:|---:|
| 17 | 5 | 13 | 5 | 4 |
| 19 | 6 | 14 | 6 | 4 |
| 23 | 7 | 17 | 7 | 1 |
| 29 | 9 | 21 | 9 | 0 |

Thus the first possible positive slack is at least `6,7,8,10` at
`p=17,19,23,29`, respectively. The first prime whose endpoint is not
excluded by this proposition is `p=31`, where `R=10`.

## Universal endpoint proof

Write `p=3R+c`, where `c in {1,2}`. Proposition 15.726 gives the linewise
upper bound

`I=sum_(z in T) s_A(z) <= R`.

Suppose `t<R`. The Ball--Lavrauw tangent envelope still applies because at
the worst value `t=R-1`,

`|A|-(2(t+1)+2)=p-3t-3=c>=1`.

Exactly as in Proposition 15.726 it gives

`I>=F(t)=t(p-1-3t)/2`.

The function `F` is concave. On `1<=t<=R-1`, its endpoint margins above
`R` are

| residue | `F(1)-R` | `F(R-1)-R` |
|---|---:|---:|
| `p=3R+1` | `(R-3)/2` | `(R-3)/2` |
| `p=3R+2` | `(R-2)/2` | `R-2` |

They are positive for every prime `p>=17`. Hence `t<R` is impossible and
`t=R`.

Minimum-cardinality makes `T` inclusion-minimal, so every integer
`s_A(z)>=1`. Combining `|T|=R` with `I<=R` forces

`I=R` and `s_A(z)=1` for all `z in T`.

On a line with `a<=2` retained arc points and `u` deleted points, its
contribution to `I` is `u` when `a=2` and zero otherwise. Equality with the
slack contribution forces either total occupancy at most two, or

`(a,u)=(2,1),(2,2)`.

These are exactly a trisecant and a 4-secant. If two rich lines shared a
point of `D`, choose that shared point in both deletion demands `n_l-2`, and
choose arbitrary required deletions on every other rich line. Their union
would repair `D` to an arc with at most

`sum_l(n_l-2)-1=R-1`

deletions, contradicting the minimum `t=R`. Thus the rich lines are
`D`-disjoint and `x+2y=R`.

The complete line census is

`N_4=y`, `N_3=R-2y`, `N_2=p(p+1)/2-3R`,

`N_1=p+1+3R+2y`, `N_0=p(p-1)/2-R-y`.

Removing every point outside the rich blocks and one point from each
4-secant leaves a regular trisecant core of size `3(R-y)`: every point lies
on exactly one trisecant and has `c+3+3y` tangents. This is the next
all-prime geometric target; no general semiarc theorem is claimed here.

## Classification-assisted closes

### `p=17`

The endpoint arc has size 13 and needs five index-one outside points.
Sticker's eight complete 13-arc classes have index-one counts

`0,0,0,0,0,0,2,3`.

If the 13-arc is incomplete, extend it to a 14-arc. If that extension is
complete, Sticker's unique complete-14 class applies; auditing all fourteen
13-subarcs gives index-one counts at most four. If the 14-arc is incomplete,
extend once more to a 15-arc. The unique 15-arc class is conic-contained. A
missing conic point has secant index zero, while an off-conic point retains
at least `8-5=3` secants after five conic points are omitted. Hence this
branch has no index-one point.

### `p=19`

The endpoint arc has size 14 and needs six index-one outside points.
Al-Zangana's exhaustive table of all 83 projective 14-arc classes, including
13 incomplete and 70 complete classes, has `c_1<=4` in every class.

### `p=23`

The endpoint arc has size 17 and needs seven index-one outside points.
The five complete 17-arc classes have index-one counts `0,0,1,0,0`. If the
arc is incomplete, extend it to size 18. The Coolsaet--Sticker complete-arc
spectrum has no complete size from 18 through 23 and a unique size-24 class,
the conic. A missing conic point has secant index zero, while an off-conic
point retains at least `11-7=4` secants. Thus this branch also has no
index-one point.

### `p=29`

The endpoint arc has size 21 and needs nine index-one outside points.
Coolsaet--Sticker's full complete-arc classification has exactly two
complete 21-arc classes. Direct audits of their published representatives
give outside secant-index histograms

`{4:18,5:75,6:190,7:312,8:189,9:63,10:3}`

and

`{3:3,4:21,5:66,6:187,7:294,8:243,9:27,10:9}`.

Thus both have `c_1=0`. If the 21-arc is incomplete, extend it to a complete
arc. Segre's odd-order bound caps its size at 30, and the published spectrum
has no complete sizes 22, 23, or 25 through 29, so the extension has size 24
or 30. The unique complete 24-arc is the Klein quartic
`x^3y+y^3z+z^3x=0`; its exact outside secant-index histogram is

`{6:28,8:126,9:504,10:84,11:84,12:21}`.

A 21-subarc deletes three Klein points. Missing Klein points have index zero,
and every point outside the Klein arc retains at least `6-3=3` secants. The
unique complete 30-arc is a conic. A 21-subarc omits nine conic points, so a
missing conic point has index zero and an off-conic point retains at least
`14-9=5` secants. Hence every incomplete branch also has `c_1=0`.

These are classification-assisted theorems, not new finite solver runs. The
`p=29` coordinate and Klein-curve checks are small exact enumerations used to
audit the published representatives, not a search for configurations.

## Artifacts

- `src/e1_gmin_m4_prop15727.py`
- `tests/test_prop15727.py`
- `evidence/e1_gmin_m4_prop15727.json`

The imported classification certificates remain in Propositions 15.685,
15.693, and 15.701--15.703.
