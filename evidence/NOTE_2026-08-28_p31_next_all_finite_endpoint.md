# The `p=31` next all-finite endpoint is impossible

Date: 2026-08-28. This is Proposition 15.682. It closes the `p=31,s=26`
second all-finite endpoint by combining Proposition 15.681's integral
paired-cube lift, an exact fourteen-profile ledger, and Coolsaet's exhaustive
classification of complete arcs in `PG(2,31)`. The same boundary remains
open at `p=17,19,23,41`; `p=41` has only residue zero left. Later
all-finite sizes, the infinity-present remainder, residual (ii), R1, global
QVAR, Type I, and the limit remain open.

## 1. Only residue zero remains

At `p=31`, put

```text
m=16, P=32, s=26, pair-deficit budget=s(s-1)=650.
```

The exact quotient/floor ledger has phase-one residue `u_1=15` and
phase-zero pair survivors

```text
u_0 = 0,2,3,4,5,6.
```

Every positive residue forces a quotient-zero `b=0` direction with scaled
quadratic-lift mass at most 12. Proposition 15.681 gives the stronger floor

```text
4p E[B] >= (p+1)/2 = 16,
```

so all positive rows are impossible. Only `u_0=0` remains.

## 2. Exact residue-zero profiles

The two type minima are

```text
phase zero: 10*b=0 + 1*b=2 + 5*b=26, deficit 284,
phase one:  15*b=2 + 1*b=26,          deficit 360.
```

Their nominal total is 644, leaving six pair-deficit units. The exact
enumeration includes a second near-minimal phase-zero quotient allocation,
so merely lowering the six high fibre counts is not exhaustive. Enumerating
every profile through deficits 290 and 366 gives fourteen phase-labelled
profiles.

For an affine line of occupancy `n`, the pair-slack contribution is

```text
2*(C(n,2)-floor(n/2)),
```

a nonnegative multiple of four. Thus the total slack can only be four or
zero. Writing `t=(26-b)/2`, the fourteen profiles have seven distinct
global distributions:

```text
slack 4:
  {t=13:10, t=12:16, t=1:1, t=0:5},
  {t=13:11, t=12:15,        t=0:6};

slack 0:
  {t=13:10, t=12:16, t=3:1,            t=0:5},
  {t=13:10, t=12:16, t=2:1, t=1:1,     t=0:4},
  {t=13:10, t=12:16,        t=1:3,     t=0:3},
  {t=13:11, t=12:15, t=2:1,            t=0:5},
  {t=13:11, t=12:15,        t=1:2,     t=0:4}.
```

There are three phase-labelled slack-four profiles and eleven equality
profiles. Slack zero means that the 26 finite points form an arc; every
such profile has at least three undetermined directions. Slack four means
that exactly one affine line has occupancy three and every other line has
occupancy at most two; those profiles have at least five undetermined
directions. Deleting any point of the unique triple gives a 25-arc and
preserves all undetermined directions.

## 3. Every 27-/28-arc in `PG(2,31)` lies on a conic

Coolsaet's full isomorph-free classification of complete arcs in
`PG(2,31)` has no complete arc of any size from 23 through 31. Every arc in
a finite projective plane can be extended greedily until it is complete.
Therefore a 27- or 28-arc must extend to a complete 32-arc. By Segre's
odd-order theorem, every 32-arc in `PG(2,31)` is a nondegenerate conic.
Consequently every 27- and every 28-arc is conic-contained.

This implication is also summarized by the standard parameter
`m'(2,31)=22`: 22 is the largest size of a nonconic complete arc.

## 4. Three undetermined infinity points contradict a conic

In a slack-zero profile, adjoin any two undetermined infinity points to the
26-arc, producing a 28-arc and hence a conic subset. In a slack-four
profile, first delete one triple point and then adjoin any two undetermined
infinity points, producing a 27-arc and hence a conic subset.

Choose three distinct undetermined points `D_1,D_2,D_3` on the line at
infinity. The conics containing the extensions by `{D_1,D_2}` and
`{D_1,D_3}` share at least 25 affine arc points and therefore coincide.
That nondegenerate conic would contain the three collinear points
`D_1,D_2,D_3`, impossible. Every residue-zero profile is excluded, so the
`p=31,s=26` endpoint is closed.

## 5. Literature, OEIS, and reproduction

The complete-arc classification is an explicit external dependency:

- K. Coolsaet, *The Complete Arcs of PG(2,31)*, J. Combin. Des. **23**
  (2015), 522--533, doi:10.1002/jcd.21410.

The conic classification is Segre's classical theorem, restated in
Ball--Lavrauw's *Planar arcs*. No searched source combines these facts with
the Paley phase ledger or Proposition 15.681's integral lift.

The OEIS context check was directly useful here. A000509 records the second
largest complete-arc size and explicitly states `m'(31)=22` because no
complete `n`-arc exists for `23<=n<=31`; this pointed to the exact
classification exit above. Searches for the endpoint deficit block
`284,286,288,290,360,362,364,366` and the seven secant profiles found no
additional relevant sequence. No new sequence claim is made.

Reproduce with

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15682.py
PYTHONPATH=src pytest -q tests/test_prop15681.py tests/test_prop15682.py
```

The generated record is `evidence/e1_gmin_m4_prop15682.json`.
