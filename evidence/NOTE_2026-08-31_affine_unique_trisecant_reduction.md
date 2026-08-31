# Endpoint blocks force a large affine unique-trisecant set

**Date:** 2026-08-31
**Status:** proved structural reduction, not endpoint exclusion.

## Hypothesis inherited from Proposition 15.727

Let `p>=17` be prime and write

`R=floor((p-1)/3)`, `p=3R+c`, `c in {1,2}`.

At endpoint equality, `D` consists of `p+1` affine points.  Its rich lines
are pairwise disjoint as subsets of `D`: there are `x` trisecants and `y`
4-secants, and

`x+2y=R`.                                                   (1)

These are all lines containing at least three points of `D`.

## Retain three points on exactly one rich block

Choose any one rich block `L`.  Keep three points of `L`, keep two points
on every other rich block, and keep every point outside the rich blocks.
Call the retained set `U`.

If `L` is a trisecant, the number deleted is

`(x-1)+2y=R-1`.                                             (2)

If `L` is a 4-secant, delete one point of `L`; the number deleted is

`1+x+2(y-1)=R-1`.                                          (3)

Thus in either case

`|U|=(p+1)-(R-1)=p+2-R`.                                   (4)

The set `U` is affine because `U` is a subset of the already-affine set
`D`.  Any line containing at least three points of `U` also contains at
least three points of `D`, so it is one of the rich blocks above.  The
chosen block retains three points and every other rich block retains two.
Therefore

`U is an affine (p+2-R,3)-arc with exactly one trisecant.`   (5)

No blocking-set theorem, configuration search, or classification is used.

## Two extension points on one tangent

Write the unique trisecant of `U` as `L={P,Q,Z}` and set

`B=U\{P,Q}`.

Then `B` is an affine arc of size

`|B|=p-R`.                                                  (6)

Both `B union {P}` and `B union {Q}` are arcs.  Indeed, if adding `P`
created a trisecant, that line would already be a second trisecant of `U`;
the same argument applies to `Q`.  Also `L` meets `B` only at `Z`, so `L`
is a tangent of `B` at `Z`.  Hence `P,Q` are two distinct affine extension
points of `B` on the same tangent.  Any of the three points of `L` may be
chosen as the one retained in `B`.

The two residue classes are

| residue | unique-trisecant set `U` | co-tangent arc `B` |
|---|---:|---:|
| `p=3R+1` | `2R+3` | `2R+1` |
| `p=3R+2` | `2R+4` | `2R+2` |

## Exact remaining implication

This proposition does not rule out either object.  A geometric endpoint
close can now target either of two equivalent necessary consequences:

- exclude the above large affine 3-arcs with exactly one trisecant, while
  retaining the fact that they arise from the common disjoint-block
  completion `D`; or
- exclude the above large affine arcs having two distinct extension points
  on one tangent, again with the compatible disjoint-block completion.

Dropping the completion data would be a weaker problem and is not justified
as an equivalent reformulation.

## Corrected tangent route

The reduction lands immediately next to two published polynomial thresholds.
Remove one point from the unique trisecant and call the resulting arc `A`.
Then

`|A|=p+1-R`, and in the standard notation `|A|=p+2-t` one has `t=R+1`.

For odd order, Ball--Lavrauw [44, Theorem 11] give an explicit
degree-`2t` tangent interpolation under the sufficient hypothesis
`|A|>=2t+2=2R+4`.  The endpoint arc misses that stated hypothesis by exactly
two points when `p=3R+1` and exactly one point when `p=3R+2`:

| residue | actual `|A|` | standard threshold | deficit |
|---|---:|---:|---:|
| `p=3R+1` | `2R+2` | `2R+4` | 2 |
| `p=3R+2` | `2R+3` | `2R+4` | 1 |

The previously recorded parallel ``unique-trisecant ceiling'' is
**RETRACTED**.  Bartoli--Storme [70, Theorem 1 in the public manuscript;
Corollary 2.7 in the published organization] use

`d <= 2(p+2)/3+2`

as the upper endpoint under their other hypotheses, including
`d>3+2sqrt(p)` and existence of the configuration, for which the arrangement
arising from a unique-trisecant `(d,3)`-arc is the second-smallest
hyperplane arrangement. It is not an upper bound on the size of a
unique-trisecant 3-arc. Their construction
section explicitly says that, for odd order, existence in a substantial
larger interval was unknown.  Consequently the former claims that `U`
``attains the ceiling'' for `p=3R+2`, and that an equality/one-defect
stability theorem for that ceiling should be proved, have no valid source.
They are not used in Proposition 15.729 itself.

The numerical deficit does **not** mean that tangent-envelope conditions are
missing. Proposition 15.730 first strengthens one chosen `U` to all
`3^x 6^y` maximum repairs. Proposition 15.731 then applies Segre's tangent
lemma and an elementary dual-line gluing argument directly. Every repair has
a degree-`2t` tangent envelope: it is unique for `p=3R+2`, while for
`p=3R+1` it is an affine line-product pencil after fixing tangent
normalization. Adjacent repair envelopes satisfy an exact transition law
with a quadratic quotient in the first residue and a cubic class in the
second.

The exact open implication is exclusion of the common completion under the
15.728/15.730 direction, phase, and lift constraints. Deriving a nontrivial
cycle identity from the low-degree transitions is the proposed next attack,
not a condition established here. It must not assume either a nonexistent
unique-trisecant size ceiling or a failure of tangent-envelope existence.

## Artifacts

- `src/e1_gmin_m4_prop15729.py`
- `tests/test_prop15729.py`
- `evidence/e1_gmin_m4_prop15729.json`
