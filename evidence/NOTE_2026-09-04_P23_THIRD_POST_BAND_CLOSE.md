# The third p23 post-band endpoint: all-boundary proof

Proposition 15.771. At `p=23,t=11,k=114`, the residual-(ii) isolated-chart
branch is empty for every boundary size. This is an **endpoint theorem**
using previously proved local theorems and a fixed finite coefficient
certificate. It is not a global residual-(ii) or original-limit theorem.

The checkpoint `f8af445a0569cd03032f71ed8e18b30822a63721` intentionally left
this endpoint `REVIEW_PENDING`. The three then-missing bridges are supplied
below and in the two linked standalone proofs. The mesh checks are
independent verification of these bridges, not replacements for them.

## 1. The exhaustive residue split

Here `q=11`, `m=12`, and `|H|=k+1=115`. There are at least
`23^2+1-2*115=300` isolated vertices, so the proved signed transport gives
an all-finite chart with `I=0` and every directional boundary size even.
Write `h` for the hard-direction sign, `T` for the common signed edge sum,
and `P_L` for the hard parallel counts. The exact row identities are

```text
a_L = 2u+24 k_L,        sum_hard k_L = 23-u,       0<=u<=11,
hT  = 24 P_L-69-a_L.                                      (1)
```

The phase-one floors are 22 for `b=2,22`, and 46 for all other even
boundary sizes. At `b=2,22`, subtract the respective pointwise parity
minima `(1-x_i-x_j)^2` and `1-x_j`. Half the difference is a nonnegative
integral quadratic. Every positive such lift has scaled mass at least
20 by Proposition 15.688.

A quotient-zero row is therefore impossible unless `u=11`, where it must
be an exact `b=2` or `b=22` baseline. In the absence of quotient zero,
there are at least `24-(23-u)=u+1` quotient-one rows. Their scaled mean is
`24+2u`. For `u<=8`, the only possible boundaries would leave positive
lift mass `2+2u<20`, a contradiction. The complete residue list is thus

| u | Low hard rows | Low scaled mean | Task |
|---|---:|---:|---|
| 9 | at least 10 | 42 | sharp mass-20 lift, two-unit carry |
| 10 | at least 11 | 44 | forbidden mass-22 lift |
| 11, some quotient zero | at least one exact baseline | 22 | common baseline ledger |
| 11, no quotient zero | all 12 | 46 | full equality classification |

## 2. The carried residues u=9 and u=10

For `u=9`, reuse the sharp equality families and integral coefficient
comparison of Propositions 15.769--15.770. Their offsets are `2,4,3,5`.
Equation (1) makes `P_L-k_L` common. Relative to a low row's parallel
count `P`, the hard parallel total is `12P+2`; the slice-kernel congruence
and the total-edge bound fix `P` to the family's offset. Therefore

```text
hT=24P-111,
sum_opposite Q=113-12P,
a_opposite(Q)=24(P+Q)-180.
```

The putative `Q=8-P` row has mass12 and is impossible by the phase-zero
floors and sharp lift bound. Thus every opposite count is at least `9-P`.
Above those twelve minima the total surplus is five, so at least seven
rows have scaled mass36. Their nonzero boundaries are excluded by the
same pointwise subtraction: only `b=2,22` survive the floors, and their
lift mass12 is below20. At `b=0`, reuse the proved mean-half equality
classification `F4/F5`. The coefficient congruence leaves only
`P=4,Q=5,F5`. The ten low hard directions are ten projective roots of
the common quartic and octic forms, so both forms vanish identically.
The fixed exact 33,649-five-set certificate from Proposition 15.769 has
no simultaneous quartic/octic zero for `F5`, a contradiction. That
certificate is reused, not enlarged or rerun as a new search.

For `u=10`, every low row must have `b=2,22`; subtracting the exact parity
minimum leaves a nonnegative integral quadratic of mass22. The local
`p-1` exclusion in the already proved second-p23 endpoint package
excludes it. At least eleven low rows exist, so one suffices.

## 3. Every hard equality at u=11

Suppose first that there is no quotient-zero row. There are twelve
nonnegative integer quotients summing to twelve, so every quotient is
one. Every hard row has `E[A]=1`, with scaled mean46.

* At `b=0`, parity makes `A` odd and at least one everywhere; its mean
  forces `A=1`, giving signed-target offset5.
* At `b=2,22`, subtract the exact parity baseline and divide by two.
  The lift has mass24. The sharp stabilizer and paired-cube bounds force
  height6 if the height is at least two. Every maximizing cube would
  then have mean1/2, contradicting the dimension-free maximum3 theorem.
  The lift is Boolean. The corrected Johnson influence bound is
  `15708/2645<6`, reducing to at most five slice coordinates and hence
  the fixed four-bit Boolean quadratic catalog. Its 30 relevant tables
  are the selected pair, oriented pair, and compact triangle families.
  On the original slice their `4L` offset increments are `3,1,1`.
  Adding baseline offsets4 or3 gives offsets7,5 or6,4.
* At `b=4,20`, positive contact quadrature makes `A=1` on all the
  required even intersection layers. The general-slice reduction in
  [the small-boundary proof](NOTE_2026-09-04_P23_SMALL_BOUNDARY_EQUALITY_PROOF.md)
  removes all outside-coordinate dependence. The complete types are
  `4000/2200` at `b=4` (offset5), and `000;4`, `200;2`, `220;0`,
  `400;0` on the three-point complement at `b=20` (offsets8,6,4,4).
* At `b=6,8,...,18`, the positive quadrature forces `A-1=0` on every
  even intersection. Through every twelve-set there is a cross-boundary
  swap cube of dimension `min(b,23-b)>=5`. Restriction of degree-two
  functions to its even half is injective by the exact Walsh character
  Gram identity. It forces `A=1` at an odd vertex, violating parity.
  [The covering-cube proof](NOTE_2026-09-04_P23_MIDDLE_BOUNDARY_CUBE_PROOF.md)
  includes the injection capacities and extreme intersection layers.

This accounts for every even boundary, and every surviving family has an
explicit integral signed-target representative with offset in `4,5,6,7,8`.

## 4. The all-one common-row contradiction

Equation (1) with `a_L=46` gives a common parallel count `P` and

```text
hT=24P-115,
sum_opposite Q=115-12P,
a_opposite(Q)=24(P+Q)-184.                                (2)
```

In every row the established isolated-chart coefficient comparison gives
`P=offset (mod11)`. Since `0<=12P<=115`, we have `0<=P<=9`, and each
offset lies between4 and8. Hence `P=offset` exactly; distinct offsets
cannot mix. The `Q=8-P` value in (2) has scaled mean8 and is excluded in
Section6. Smaller counts give negative means and are impossible too.
Thus `Q>=9-P` in every opposite direction. The surplus above twelve such
counts is

```text
(115-12P)-12(9-P)=7.
```

At least five opposite directions therefore have `Q=9-P` and scaled
mean32. Section6 excludes every such row, finishing this branch.

## 5. The quotient-zero common-row contradiction

An exact quotient-zero row has boundary2 or22 and offset4 or3. For every
hard row, (1) now gives one common integer `c=P_L-k_L`, with

```text
hT=24c-91,     sum_hard P_L=12c+12.
```

The baseline row has `k_L=0`, so `c>=0`; the edge bound gives `c<=8`.
Its coefficient congruence fixes `c=4` or3 respectively. The opposite
ledger is

```text
sum_opposite Q=103-12c,
a_opposite(Q)=24(c+Q)-160.                                (3)
```

In (3), `Q=7-c` has forbidden mean8, and lower values give negative
means. Therefore `Q>=8-c`. The total surplus is again seven:
`(103-12c)-12(8-c)=7`. At least five opposite rows have mean32, again
contradicting Section6. These two cases exhaust all quotient profiles.

## 6. The phase-zero bridge, including pointwise nonnegativity

For an opposite row, write `a=46 E[A]`, with `A` a nonnegative integral
quadratic of parity `|X intersect B|`. The exact phase-zero floors are

```text
b       0   2   4   6   8  10  12  14  16  18  20  22
floor   0  24  40  46  46  46  46  46  46  46  40  24.
```

At mean8, only `b=0` survives. Then `A=2L` with `L` a nonzero
nonnegative integral quadratic and `4p E[L]=8<20`, impossible.

At mean32, only `b=0,2,22` survive. At `b=2`, the pointwise parity
minimum is the Boolean quadratic `(x_i-x_j)^2`. At `b=22`, it is `x_j`
for the omitted coordinate, since `12-x_j` has the parity of `x_j`.
Both minima have mean12/23 and scaled mean24. The difference between
`A` and its minimum is nonnegative and even at EVERY slice point.
Thus half that difference is genuinely a nonnegative integral quadratic,
of mass `32-24=8<20`, a contradiction. This cannot be justified by a
floor subtraction alone; the pointwise Boolean parity minima are essential.

The remaining `b=0` row again has `A=2L`, now with
`4p E[L]=32=p+9`. The local theorem of Proposition15.752 excludes it at
`p=23`. That theorem has no shell-band or boundary hypothesis, so its use
here does not attempt to extend the old band by assertion.

All branches are empty. The next unclosed p23 layer is `t=12,k=116`.
The generic frontiers, residual(ii) globally, E1, `L=1/2`, and original
convergence remain open.
