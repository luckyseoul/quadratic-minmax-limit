# The unique `p=23` slack-16 profile is impossible

Date: 2026-08-28. This is Proposition 15.686. It removes the second
exceptional low-slack profile left by Proposition 15.684. Exactly 201
arithmetic profiles remain at `p=23,s=20`, all of pair slack at least 20.

## 1. The exact row

The unique remaining slack-16 profile is

```text
phase 0: 7*b0 + b2 + 3*b18 + b20, deficit 164
phase 1: 11*b2 + b18,                 deficit 200
global floor-secants: {t0:1, t1:4, t9:12, t10:7}
pair slack: 380-164-200 = 16.
```

In particular, there is exactly one undetermined direction. Write `U` for
its point on the line at infinity.

## 2. Repair and the complete 17-arc

The repair lemma of 15.684 deletes at most four points from a hypothetical
20-set `S` to obtain an arc `A`. If at most three points are deleted,
`A union {U}` is an arc of size at least 18. The complete-arc
classification puts it on a conic, and the conic-core lemma contradicts
positive slack `16<24`.

The only branch left has

```text
S = A union D,   |A|=16,   |D|=4.
```

Now

```text
K = A union {U}
```

is a 17-arc. If `K` were incomplete, it could be extended to an 18-arc
and would give the same conic-core contradiction. Therefore `K` is a
complete 17-arc.

## 3. Four multiplicity-one points are forced

For `d in D`, the line `Ud` contains no second point of `S), because
`U` is undetermined for `S`. Consequently no secant of `K` through
`d` uses `U): every such secant is already a secant of `A). Thus

```text
mu_A(d) = mu_K(d) >= 1,
```

where the inequality is completeness of `K`.

As in Proposition 15.685, a secant of `A` containing `r` deleted points
contributes at least `4r` pair slack. Hence

```text
16 = slack(S) >= 4 * sum_{d in D} mu_A(d) >= 4*4.
```

Equality is forced throughout, so all four deleted points satisfy

```text
mu_K(d)=1.
```

Proposition 15.685 verified explicit representatives for all five
Coolsaet--Sticker classes of complete 17-arcs in `PG(2,23)`. Their numbers
of outside points of secant multiplicity one are

```text
0, 0, 1, 0, 0.
```

No class has the required four points. This excludes the profile.

## 4. Exact remainder

The remainder is now

```text
{20:68, 24:49, 28:35, 32:21, 36:13, 40:7,
 44:4, 48:1, 52:1, 56:1, 60:1},
```

which sums to 201 profiles. The `p=23` endpoint remains open.

## 5. Context and reproduction

The external classification is K. Coolsaet and H. Sticker, *A full
classification of the complete k-arcs of PG(2,23) and PG(2,25)*,
J. Combin. Des. **17** (2009), 459--477,
doi:10.1002/jcd.20211.

Targeted literature and OEIS searches found no prior statement matching this
four-point repair obstruction or the exact remaining histogram. No novelty
or sequence-submission claim is made.

Reproduce with

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15686.py
PYTHONPATH=src pytest -q tests/test_prop15686.py
```

The machine-readable record is
`evidence/e1_gmin_m4_prop15686.json`.
