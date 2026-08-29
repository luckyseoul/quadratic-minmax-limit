# p=19 slack-16 complete-14-arc exclusion

**Date:** 2026-08-29
**Proposition:** 15.693
**Status:** proved conditional on the exhaustive published classification;
reduces the p=19 endpoint from fourteen profiles to seven

All seven slack-16 profiles left by Proposition 15.689 have three or four
undetermined directions. Repair deletes at most four boundary points.

With at most three deletions, adjoining two undetermined infinity points
gives an arc of size at least 15. The `PG(2,19)` complete-arc spectrum puts
it on the conic, where a third undetermined infinity point is impossible.

Suppose all four deletions are needed. Write `S=A union D`, where `A` is a
12-arc and `|D|=4`, and adjoin two undetermined infinity points `U1,U2`:

```text
K = A union {U1,U2},   |K|=14.
```

If `K` is incomplete, it extends through the classified size gap to the
20-point conic. A third undetermined infinity point lies off that conic and
has exactly one `K`-secant, the line at infinity. But an off-conic point has
at least nine conic secants, and the six omitted conic points destroy at
most six of them. At least three `K`-secants remain, a contradiction. Hence
`K` is complete.

For `x in D`, no secant of `K` through `x` can use `U1` or `U2`, because
their directions are undetermined for the original boundary. Every such
secant is therefore an `A`-secant and charges four units of pair slack.
Completeness gives at least one secant through each deleted point, while
total slack 16 forces exactly one through each. Thus all four deleted
points have secant index one outside `K`.

Every unused undetermined infinity point also has secant index exactly one:
its only `K`-secant is the line at infinity through `U1,U2`. Therefore `K`
has at least

```text
4 + (t-2) >= 5
```

outside points of secant index one.

Al-Zangana's exhaustive `PG(2,19)` classification has 83 projective
14-arc classes, 70 complete, and reports `c1<=4` for every 14-arc, where
`c1` counts outside points of secant index one. Sticker's independent
classification gives the same complete-class count. This contradiction
excludes all seven slack-16 profiles.

The exact p19 remainder is now

```text
{20:4, 24:1, 28:1, 32:1}  (seven profiles total).
```

The same argument shows that every slack-20 profile must use exactly five
repair deletions and every slack-24 profile at least five. It does not yet
exclude those profiles or close residual (ii).

Primary classification input: E. B. Al-Zangana, *The Geometry of the Plane
of Order Nineteen and its Application to Error-Correcting Codes*, PhD
thesis, University of Sussex, 2011, Chapter 4, Section 4.22, p. 105.

Reproduction:

- `src/e1_gmin_m4_prop15693.py`
- `evidence/e1_gmin_m4_prop15693.json`
- `tests/test_prop15693.py`
