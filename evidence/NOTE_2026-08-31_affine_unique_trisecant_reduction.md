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

## Specific next attack: the one/two-defect tangent interpolation

The reduction lands immediately next to two published polynomial thresholds.
Remove one point from the unique trisecant and call the resulting arc `A`.
Then

`|A|=p+1-R`, and in the standard notation `|A|=p+2-t` one has `t=R+1`.

For odd order, the Segre--Blokhuis--Bruen--Thas tangent-form interpolation
reviewed by Ball--Lavrauw [44] has its usual uniqueness range at
`|A|>=2t+2=2R+4`.  The endpoint arc misses that range by exactly two points
when `p=3R+1` and exactly one point when `p=3R+2`:

| residue | actual `|A|` | standard threshold | deficit |
|---|---:|---:|---:|
| `p=3R+1` | `2R+2` | `2R+4` | 2 |
| `p=3R+2` | `2R+3` | `2R+4` | 1 |

There is a parallel extremality statement for `U`.  Bartoli--Storme [70,
Corollary 2.7] bound an odd-order 3-arc with a unique trisecant by
`floor(2(p+2)/3+2)`.  Thus `U` is one below that ceiling for `p=3R+1` and
attains it for `p=3R+2`:

| residue | actual `|U|` | unique-trisecant ceiling | deficit |
|---|---:|---:|---:|
| `p=3R+1` | `2R+3` | `2R+4` | 1 |
| `p=3R+2` | `2R+4` | `2R+4` | 0 |

The recommended next lemma is therefore not a broad classification.  It is
an equality/one-defect stability theorem for that ceiling, or equivalently
an endpoint-specific repair of the missing one or two tangent-interpolation
conditions.  The extra hypotheses available here are precisely what a
general `(k,3)`-arc lacks: `U` is affine, it comes with all `R` index-one
points of the same arc repair, and every one of those secants belongs to the
common pairwise-disjoint 15.727 completion `D`.  A proof that these data force
the standard tangent form, or classify the equality/one-defect cases into a
non-affine or conic completion, would close endpoint equality for every
prime at once.

This paragraph records a route, not a theorem.  No such stability statement
is claimed here.

## Artifacts

- `src/e1_gmin_m4_prop15729.py`
- `tests/test_prop15729.py`
- `evidence/e1_gmin_m4_prop15729.json`
