# p=17 slack-sixteen free-direction exclusion

**Date:** 2026-08-29
**Proposition:** 15.704
**Status:** proved conditional on Sticker's published complete-arc class
counts; excludes 99 profiles but does not close the endpoint

The 112 slack-sixteen profiles left by Proposition 15.703 have undetermined-
direction histogram

```text
{0:13, 1:47, 2:47, 3:5}.
```

Pair slack sixteen permits repair by at most four deletions. With two
undetermined directions, adjoining their infinity points gives a 14-arc in
the depth-four branch. If complete, the unique complete-14 class has outside
secant index at least two, so the four deleted points force slack at least
32. If incomplete, it extends to a conic and the positive off-conic slack
floor is at least 20. Three directions are handled by two overlapping pairs:
both pair arcs must extend to the same conic, which cannot contain the three
collinear infinity points.

For one direction, depths below four similarly reach either a conic or a
complete 14-arc whose three deleted points force slack at least 24. At depth
four, adjoining the infinity point gives a 13-arc. If complete, slack sixteen
would require four outside points of secant index one, whereas the eight
complete-13 classes have index-one counts

```text
0,0,0,0,0,0,2,3.
```

If the 13-arc is incomplete, extend it by one point. The incomplete-14 branch
reaches a conic. In the complete-14 branch, including the extension point
among the deleted points forces slack at least 24. Otherwise, deleting each
of the 14 complete-arc points gives index-one-count histogram
`{0:4,1:8,4:2}` and two candidate index-one quadruples. Across all 26 choices
of infinity point their reconstructed slack histogram is

```text
{16:2, 28:16, 32:8}.
```

The two apparent slack-sixteen placements are not undetermined directions.
After imposing the exact condition that no boundary chord pass through the
chosen infinity point, eight placements remain and all have slack 32.

Therefore all 99 profiles with at least one undetermined direction are
impossible. The exact p17 remainder falls from 753 to 654:

- two tangent-conic slack-zero profiles;
- thirteen slack-sixteen profiles with no undetermined direction;
- 639 profiles of slack at least twenty.

The p17 endpoint and every top-level gate remain open.

## Search audit

Targeted literature searches found the published arc classifications but no
prior statement of this complete-14-minus-one infinity-placement census or
the 99-profile exclusion. OEIS searches found only incidental appearances of
the isolated counts in unrelated sequences; those matches are not evidence.

Reproduction:

- `src/e1_gmin_m4_prop15704.py`
- `evidence/e1_gmin_m4_prop15704.json`
- `tests/test_prop15704.py`
