# p=17 slack-twelve complete-13-arc exclusion

**Date:** 2026-08-29
**Proposition:** 15.703
**Status:** proved conditional on Sticker's published complete-class count;
closes the slack-twelve block but not the endpoint

Proposition 15.702 leaves 33 slack-twelve profiles, all with no undetermined
direction. Every realization must repair by exactly three deletions to a
13-arc `A`.

Sticker's table has eight complete 13-arc classes in `PG(2,17)`, with
stabilizer orders

```text
1,2,2,2,2,3,4,6.
```

The exact normalized PGL generator fixes a quadrangle, finds a complete arc,
enumerates and blocks every normalized image of its full projective orbit,
and repeats. It found eight pairwise invariant-distinct representatives with
exactly that stabilizer fingerprint. Their outside index-one point counts are

```text
0,0,0,0,0,0,2,3.
```

Slack-twelve equality requires all three deleted points to have secant index
one relative to `A`. Only one complete class even supplies a triple, and that
triple reconstructs slack 16 rather than 12.

If `A` is incomplete, extend it to a 14-arc. An incomplete 14-arc extends to a
conic-contained 15-arc and is already impossible. For the unique complete
14-arc, deleting each possible point gives index-one counts

```text
{0:4, 1:8, 4:2}.
```

There are eight candidate index-one triples, and every one reconstructs slack
20. If the original deleted set includes the 14-arc extension point, the
other two deleted points have complete-14-arc secant index at least two and
force slack at least 16 directly.

Thus all 33 slack-twelve rows are impossible. The exact p17 remainder falls
from 786 to 753:

- two tangent-conic slack-zero profiles;
- 751 profiles of slack at least sixteen.

The p17 endpoint and all top-level gates remain open.

External input: H. Sticker, *Classification of Arcs in Small Desarguesian
Projective Planes*, PhD thesis, Ghent University (2012), Section 5.1, printed
page 102 (PDF page 111).

## Search audit

Targeted searches found the published eight-class count but no prior listing
of the eight outside secant-index histograms, their index-one count vector
`0,0,0,0,0,0,2,3`, or this slack-twelve exclusion. OEIS searches found no
matching sequence. The web results are not used in place of the local class
and triple census.

Reproduction:

- `scripts/p17_complete_arc_class_generator.py`
- `evidence/e1_gmin_m4_prop15703_complete13_classes.json`
- `src/e1_gmin_m4_prop15703.py`
- `evidence/e1_gmin_m4_prop15703.json`
- `tests/test_prop15703.py`
