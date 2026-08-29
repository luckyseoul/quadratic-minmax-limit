# p=17 complete-14-arc low-slack exclusion

**Date:** 2026-08-29
**Proposition:** 15.702
**Status:** proved conditional on Sticker's exhaustive complete-arc
classification; reduces but does not close the endpoint

Proposition 15.701 leaves 932 exact profiles at `p=17,s=16`. The unresolved
low-slack portion is 67 slack-eight rows with no undetermined direction and
112 slack-twelve rows, split by the number `t` of undetermined directions as
`{0:33,1:79,2:43}`.

Sticker's complete-arc table records exactly one complete 14-arc class in
`PG(2,17)`, with automorphism group `D8`. The checked representative in
`src/e1_gmin_m4_prop15702.py` has line occupancy census

```text
{0:146, 1:70, 2:91}
```

and outside-point secant-index census

```text
{2:4, 3:4, 4:76, 5:128, 6:75, 7:6}.
```

Thus the unique complete class has no outside point of secant index zero or
one.

At slack eight, a realization not already excluded by the conic-core argument
must repair by exactly two deletions to a complete 14-arc. Each deleted point
then lies on at least two core secants, so the line-slack incidence bound gives
`slack >= 4(2+2)=16`, contradiction. This removes the remaining 67 rows and
closes the entire slack-eight block.

At slack twelve, a surviving repair must use exactly three deletions. For the
79 rows with one undetermined direction, adjoin its infinity point to the
repaired 13-arc. If the resulting 14-arc is incomplete, it extends to a
conic-contained 15-arc and Proposition 15.701 applies. If it is complete, the
undetermined direction ensures that a line from any deleted point to the
adjoined infinity point contains no repaired affine point. Every secant through
a deleted point therefore uses two repaired affine points. The minimum index
two forces `slack >= 4*3*2=24`, again a contradiction.

Another `67+79=146` profiles are excluded. The exact p17 remainder falls from
932 to 786:

- two tangent-conic profiles of slack zero;
- 33 slack-twelve profiles with no undetermined direction;
- 751 profiles of slack at least sixteen.

The p17 endpoint and all top-level gates remain open.

External input: H. Sticker, *Classification of Arcs in Small Desarguesian
Projective Planes*, PhD thesis, Ghent University (2012), Section 5.1, printed
page 102 (PDF page 111).

## Search audit

Targeted searches recovered Sticker's unique complete-class count and `D8`
stabilizer, but found no prior listing of the outside secant-index histogram
above or this `67+79` profile exclusion. OEIS searches found no matching count
sequence. No search result is used in place of the local exhaustive census.

Reproduction:

- `src/e1_gmin_m4_prop15702.py`
- `evidence/e1_gmin_m4_prop15702.json`
- `tests/test_prop15702.py`
