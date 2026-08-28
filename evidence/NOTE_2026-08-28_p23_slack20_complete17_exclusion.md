# All 68 `p=23` slack-20 profiles are impossible

Date: 2026-08-28. This is Proposition 15.687. It reduces the exact
`p=23,s=20` remainder from 201 profiles to 133. Every remaining profile
has pair slack at least 24.

## 1. Exact profile split

The 68 exact slack-20 profiles have the following numbers of undetermined
directions:

```text
t0=2:  2 profiles
t0=3: 36 profiles
t0=4: 30 profiles.
```

The repair lemma deletes at most five points from a hypothetical 20-set to
obtain an arc.

## 2. Conic-core count through five deleted points

Proposition 15.684 used the off-conic count through four deleted points. The
same proof remains valid for five. If a 20-set has `h` points off a conic,
an off-conic point has at least 11 full conic secants. The `20-h` retained
conic points omit `4+h` conic points, so at least `7-h` secants retain
both endpoints. Summing the line-slack charge gives

```text
slack >= 4h(7-h).
```

For `1<=h<=5`, the values are

```text
h:       1  2  3  4  5
floor:  24 40 48 48 40.
```

Thus a positive-slack conic-core set cannot have slack 20.

## 3. The 66 profiles with at least three undetermined directions

Choose three undetermined infinity points `U1,U2,U3`, but adjoin them only
in the overlapping pairs `{U1,U2}` and `{U1,U3}`.

If repair uses at most four deletions, each pair gives an arc of size at
least 18 and hence lies on a conic. The two conics share the repaired arc,
which has at least 16 points, so they coincide. That conic would contain
the three collinear infinity points, impossible.

If all five deletions are needed, each pair gives a 17-arc. If either pair
arc were complete, the slack argument in the next section would force all
five deleted points to have secant multiplicity one, contradicting the
classified maximum one. Thus both pair arcs are incomplete and extend to
18-arcs on conics. Those conics again share the repaired 15-arc, coincide,
and contain all three collinear infinity points.

This excludes all

```text
36+30=66
```

profiles with three or four undetermined directions.

## 4. The two profiles with two undetermined directions

If repair uses at most four deletions, adjoining both infinity points again
gives an arc of size at least 18 and the conic contradiction.

In the hard branch, repair uses five deletions:

```text
S=A union D,  |A|=15, |D|=5.
```

Adjoin the two undetermined infinity points to obtain a 17-arc `K`. If
`K` were incomplete, it would extend to an 18-arc and give the same
conic-core contradiction. Thus `K` is complete.

No secant of `K` through a deleted point can use either infinity point:
that would contradict the corresponding direction being undetermined for
`S`. Hence every deleted point has at least one secant already in `A`.
The exact secant-line inequality gives

```text
20 = slack(S) >= 4 * sum_{d in D} mu_A(d) >= 4*5.
```

All five deleted points must therefore have secant multiplicity one outside
`K`. Proposition 15.685 exhausts the five complete-17-arc classes and
shows that their maximum number of such points is one. This excludes the
last two profiles.

## 5. Exact remainder

The remaining histogram is

```text
{24:49, 28:35, 32:21, 36:13, 40:7,
 44:4, 48:1, 52:1, 56:1, 60:1},
```

which sums to 133. The endpoint remains open.

## 6. Context and reproduction

The external classification remains Coolsaet--Sticker,
doi:10.1002/jcd.20211. Targeted literature and OEIS searches found no
matching statement or sequence for the exact obstruction or remainder.
No novelty or sequence-submission claim is made.

Reproduce with

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15687.py
PYTHONPATH=src pytest -q tests/test_prop15687.py
```

The machine-readable record is
`evidence/e1_gmin_m4_prop15687.json`.
