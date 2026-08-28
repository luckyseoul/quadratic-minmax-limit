# The unique `p=23` slack-12 profile is impossible

Date: 2026-08-28. This is Proposition 15.685. It removes one of the 203
profiles left by Proposition 15.684, leaving exactly 202. It does not close
the `p=23,s=20` endpoint.

The proof uses the slack-to-arc repair lemma from 15.684 and the five
classified complete 17-arcs in `PG(2,23)`.

## 1. The exact profile

The unique remaining pair-slack-12 row is

```text
phase 0: 8*b0 + 4*b18, deficit 168
phase 1: 11*b2 + b18,  deficit 200
global floor-secants: {t1:5, t9:11, t10:8}
pair slack: 380-168-200 = 12.
```

For a projective line containing `n` boundary points, recall

```text
delta(n) = 2*(C(n,2)-floor(n/2)).
```

Proposition 15.684 proves that pair slack `4r` permits deleting at most
`r` points to obtain an arc. Thus a hypothetical 20-point realization
`S` of this profile has an arc `A=S\D` with `|D|<=3`.

## 2. Why the repaired arc must be a complete 17-arc

If `|D|<=2`, then `A` has at least 18 points. Coolsaet--Sticker's
complete-arc classification forces `A` onto a conic. The conic-core lemma
of 15.684 says that a positive-slack 20-set differing from that conic core
in at most four points has pair slack at least 24. This contradicts slack
12.

It follows that `|D|=3` and `A` is a 17-arc. If `A` were incomplete,
one point could be adjoined to make an 18-arc, again forcing `A` onto a
conic and giving the same contradiction. Therefore `A` must be complete.

## 3. Slack forces three points of secant multiplicity one

For a point `x` outside `A`, write

```text
mu_A(x) = number of secants of A through x.
```

Consider one secant of `A` containing `r` deleted points. Its final
occupancy in `S` is `2+r`, and

```text
r   delta(2+r)   4r
1       4         4
2       8         8
3      16        12
```

Hence, after summing over the secants of `A`,

```text
slack(S) >= 4 * sum_{x in D} mu_A(x).                (1)
```

Completeness means every outside point lies on at least one secant, so
`mu_A(x)>=1` for all three points of `D`. Since the left side of (1) is
12, equality is forced:

```text
mu_A(x)=1 for each x in D.                            (2)
```

A realization therefore requires a complete 17-arc with at least three
outside points of secant multiplicity one.

## 4. The five complete 17-arc classes

Coolsaet--Sticker's exhaustive classification states that `PG(2,23)` has
exactly five projective classes of complete 17-arcs. The machine record
contains five explicit homogeneous-coordinate representatives. The verifier
checks, for each representative:

- 17 distinct projective points;
- 136 distinct secant lines, so no three arc points are collinear;
- all 536 outside projective points lie on at least one secant, so the arc
  is complete;
- all outside-point secant multiplicities and their total incidence count
  `136*(23-1)=2992`.

The resulting histograms `{multiplicity: outside-point count}` are

```text
class 1: {2:2, 3:6, 4:68, 5:172, 6:190, 7:86, 8:12}
class 2: {2:1, 3:15, 4:59, 5:159, 6:208, 7:86, 8:8}
class 3: {1:1, 3:6, 4:69, 5:171, 6:196, 7:78, 8:15}
class 4: {3:14, 4:58, 5:170, 6:206, 7:72, 8:16}
class 5: {2:1, 3:8, 4:63, 5:185, 6:176, 7:91, 8:12}.
```

This histogram is invariant under projective equivalence. The five
histograms are distinct, so the five verified representatives are
inequivalent. Since the external classification has exactly five classes,
they exhaust it.

Their counts of multiplicity-one outside points are

```text
0, 0, 1, 0, 0.
```

No complete 17-arc has the three points required by (2). The slack-12
profile is impossible.

## 5. Exact remainder

Proposition 15.684 left

```text
{12:1, 16:1, 20:68, 24:49, 28:35, 32:21, 36:13,
 40:7, 44:4, 48:1, 52:1, 56:1, 60:1}.
```

After Proposition 15.685 the exact arithmetic remainder is

```text
{16:1, 20:68, 24:49, 28:35, 32:21, 36:13,
 40:7, 44:4, 48:1, 52:1, 56:1, 60:1},
```

which sums to 202 profiles.

## 6. Literature, OEIS, and reproduction

The external classification input is:

- K. Coolsaet and H. Sticker, *A full classification of the complete k-arcs
  of PG(2,23) and PG(2,25)*, J. Combin. Des. **17** (2009), 459--477,
  doi:10.1002/jcd.20211.
- H. Sticker, *Classification of Arcs in Small Desarguesian Projective
  Planes*, Ghent PhD thesis, 2012, Section 5.1. Its table gives
  `N_17=5` and notes that `PGL=PΓL` for prime order.

The notion of an outside point's secant index is standard in the arc
literature. Targeted searches found no published table containing these five
exact `PG(2,23)` multiplicity histograms. Exact OEIS searches for the
high-count blocks `68,172,190,86` and `69,171,196,78` returned no
relevant sequence. No novelty or OEIS-submission claim is made.

Reproduce with

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15685.py
PYTHONPATH=src pytest -q tests/test_prop15685.py
```

The generated machine-readable record is
`evidence/e1_gmin_m4_prop15685.json`.
