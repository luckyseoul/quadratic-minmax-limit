# p=17 low-positive-slack conic-core reduction

**Date:** 2026-08-29
**Proposition:** 15.701
**Status:** proved conditional on Sticker's published `PG(2,17)` arc
classification; reduces but does not close the endpoint

Proposition 15.700 leaves 1,330 exact phase-labelled profiles at the second
all-finite boundary `p=17,s=16`: two tangent-conic profiles of pair slack zero
and every positive-slack profile.

Sticker's exhaustive classification has one PGL class of 15-arcs in
`PG(2,17)`. Deleting three points from a nondegenerate conic supplies that
class, so every 15-arc is conic-contained.

Pair slack `4r` permits deleting at most `r` points to obtain an arc. An
undetermined direction (`b=16`) supplies a point at infinity that can be
adjoined while preserving the arc property. The exact ledger is:

```text
slack  profiles  t histogram             t needed  excluded  remain
    4       227  {0:113,1:102,2:12}             0       227       0
    8       195  {0:67,1:104,2:24}              1       128      67
   12       155  {0:33,1:79,2:43}               2        43     112
```

At slack four the repaired arc already has size at least 15. At slack eight,
one infinity point raises the worst repaired size from 14 to 15. At slack
twelve, two infinity points raise the worst repaired size from 13 to 15.
Every listed qualifying profile therefore has a conic core.

Let `h` original boundary points lie off that conic. Since the boundary has
16 points and the conic has 18, the retained conic subset omits `2+h` conic
points. An off-conic point has at least eight full conic secants, so at least
`6-h` still have both endpoints retained. If one such secant line contains
`a` off-conic points, its occupancy is at least `2+a` and its line-slack is at
least `4a`. Therefore

```text
pair slack >= 4 h (6-h).
```

Repair gives `1<=h<=r<=3`; the floors for `h=1,2,3` are `20,32,36`. All are
strictly above slack `4,8,12`. Thus 398 profiles are impossible, and the exact
remainder falls from 1,330 to 932:

- two slack-zero tangent-conic profiles;
- 67 slack-eight profiles with no undetermined direction;
- 112 slack-twelve profiles with at most one undetermined direction;
- all 751 profiles of slack at least sixteen.

The p17 endpoint and all top-level gates remain open. Bounded coefficient-lift
solver runs on representative slack-zero allocation branches returned
`UNKNOWN`; they are not used as evidence.

External input: H. Sticker, *Classification of Arcs in Small Desarguesian
Projective Planes*, PhD thesis, Ghent University (2012), Section 5.3, printed
page 119 (PDF page 129).

## Search audit

Targeted searches for the exact `PG(2,17)` pair-slack reduction, the profile
counts `227,195,155`, and the excluded split `227,128,43` found no prior
statement of this result. OEIS searches found no matching count sequence.
The isolated values `398`, `932`, and `1330` all happen to occur in the
Secondary Wythoff Array (A137707), but its Fibonacci/Wythoff construction has
no identified relation to this finite-geometric profile census and is not
used as evidence.

Reproduction:

- `src/e1_gmin_m4_prop15701.py`
- `evidence/e1_gmin_m4_prop15701.json`
- `tests/test_prop15701.py`
