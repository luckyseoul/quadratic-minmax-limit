# Proposition 15.770: the next two post-band residual layers close

Date: 2026-09-04

Scope: the layer immediately after Proposition 15.768 in the
`p=1 (mod 4)` class and the layer immediately after Proposition 15.769 in
the `p=3 (mod 4)` class. This is an infinite-family theorem, but it is not a
global closure of residual (ii).

## 1. Changed premise and common notation

Put

```text
q=(p-1)/2,   m=q+1,   a_L=2u+(p+1)k_L,
sum_L k_L=m+t-u.
```

The changed premise is exact: Propositions 15.768 and 15.769 now classify
and exclude the equality branches at their first post-15.752 layers. At the
next layer, an all-low equality profile becomes `m-1` low rows and one row
whose quotient is larger by one. If a low row has parallel count `P`, the
common difference-row identity is

```text
hT=(p+1)P-3p-a_L.                                      (1)
```

The high row has parallel count `P+1`. Consequently its extra hard edge and
the two extra total edges cancel in the signed difference: `hT` is unchanged,
while the opposite parallel-count total rises by one. This one-row carry is
the precise new implication used below.

No prime, graph, orbit, Johnson-cell, or residual-candidate census is used.
The sole finite enumeration is Proposition 15.751's already pinned catalog
of Boolean quadratics on four cube coordinates.

## 2. The `p=1 (mod 4)` next layer

Here

```text
t=q-2=m-3,   k=5p-5,   |H|=5p-4.
```

The phase-one residue ledger leaves exactly four arithmetic branches:

```text
u=0:     exact b=p-1 complement literal;
u=m-4:   exact b=p-3 complement triple;
u=m-3:   b=2 XNOR plus a sharp lift of mass p-3;
u=m-1:   exact b=2 XNOR endpoint.
```

All other candidate rows either have a forced quotient-zero mean below
every phase-one floor, are nonzero integral lifts from the pointwise
`b=2,p-1` parity baselines below Proposition 15.688's sharp floor `p-3`,
or are the `b=p-3` complement-triple cell at excess two. That last cell
requires the separate [punctured-gap theorem](NOTE_2026-09-04_COMPLEMENT_TRIPLE_PUNCTURED_GAP.md),
not ordinary nonnegative lift subtraction: `(r-2)^2=4` at `r=0`, while
the pointwise parity minimum there is zero. The theorem's neighboring-slice
bound, positive contact quadrature, and small-side kernel reduction give
`delta=4*sum(nonnegative integer pair values)`, excluding `delta=2`.
This supplies the previously implicit justification without changing the
four surviving branches or the endpoint conclusion.

The two old endpoint branches retain their standard forced opposite cells:

```text
hard branch             P     hT       excluded row   forced local mass
complement literal      5     p+4      Q=2, mass 6    Q=3, mass p+7
XNOR endpoint           4     5        Q=3, mass 8    Q=4, mass p+9.
```

The local masses `p+7` and `p+9` are excluded by Propositions 15.751 and
15.752.

### 2.1 Carried complement triple

Proposition 15.768's pointwise hard baseline is

```text
A=(2-|X intersect C|)^2,   |C|=3,
```

with coefficient offset two. At the new layer there are `m-1` low rows with
`P=2` and one high row with `P=3`. Formula (1) still gives

```text
hT=8-3p.
```

The hard and opposite edge totals are `p+2` and `4p-6`. An opposite row has

```text
a(Q)=(p+1)Q-6p+8.
```

The `Q=6` row has mass 14, below both the nonzero phase-zero boundary floor
and the sharp lift floor. After raising all opposite rows to `Q=7`, the
surplus is

```text
(4p-6)-7m=m-10,
```

so at least ten directions have `Q=7`. Such a row has local mass `p+15`,
excluded by Proposition 15.768.

### 2.2 New sharp XNOR branch

At `u=m-3`, all hard quotients equal one and the scaled mean is `2p-4`.
The `b=2` parity baseline has mean `p-1`, so

```text
B=(A-A_0)/2,   4p E[B]=p-3.
```

The sharp lift theorem forces `B` to be Boolean. The corrected Johnson
influence bound leaves at most five slice coordinates, and cube influence
leaves at most four active coordinates. The fixed four-bit catalog has
exactly ten tables at this density: six omitted-pair tables and four
all-equal-triple tables. Their signed-target increments have offsets `-1`
and `+1`. Adding the XNOR baseline offset four leaves exactly

```text
P=3 or P=5.
```

For either value, (1) gives

```text
hT=(p+1)P-5p+4.
```

The opposite row `Q=8-P` has mass 12 and is below both relevant floors.
The next row `Q=9-P` has mass `p+13`; the opposite surplus is `m-9`, so at
least nine such rows occur. Every nonzero-boundary alternative is a lift of
mass 12 below `p-3`. A surviving row would therefore give a nonzero
nonnegative integral quadratic `C` with

```text
4p E[C]=p+13.                                           (2)
```

### 2.3 The local `p+13` theorem includes `p=29`

For `p=1 (mod 4), p>=37`, the standard half-mean argument excludes height at
least two: it gives

```text
H >= (p-11)/4 > 3,
```

while stabilizer averaging bounds the paired-cube average by
`(p+13)/(2(p-1))<3/4`; hence some paired cube has mean `1/2` and maximum at
most three.

The first prime `p=29` is not lost. Here the half-mean step gives `H>=5`.
Thus every paired cube through a maximizer has mean at least `3/4`, giving
`H>=12`. The exact stabilizer bound is also `H<=12`, so `H=12`, and the
paired-cube average is exactly `3/4`. Proposition 15.768's dimension-free
three-quarter theorem bounds the maximum of every such cube by six, a
contradiction.

At height one, the Johnson bound at the endpoint is

```text
L <= 391608/54665 < 8.
```

The same four-bit catalog misses density `21/58`. The uniform polynomial
certifying `L<8` for every `p>=29` is

```text
p^4-25p^3+229p^2-559p+338
=x^4+91x^3+3100x^2+47204x+274272,   x=p-29.
```

Thus (2) is impossible for every claimed prime, including `p=29`. Therefore

```text
p=1 (mod 4), p>=29:  t=q-2 and k=5p-5 are empty.
```

## 3. The `p=3 (mod 4)` next layer

Here

```text
t=q-1=m-2,   k=5p-3,   |H|=5p-2.
```

The phase-one residue ledger leaves three arithmetic residues:

```text
u=m-3:   the two sharp p-3 baselines from Proposition 15.769;
u=m-2:   the same baselines with a lift of mass p-1;
u=m-1:   the two old exact endpoint baselines.
```

The old XNOR and complement-literal endpoints again force mass `p+9` after
excluding a mass-eight row.

For each of Proposition 15.769's four sharp families, the coefficient offset
is `P in {2,3,4,5}`. The new quotient ledger has `m-1` low rows of parallel
count `P` and one high row of count `P+1`. Hence

```text
hard edges=mP+1,
hT=(p+1)P-5p+4,
```

exactly the same signed sum as in Proposition 15.769. The row `Q=8-P` has
mass 12; after excluding it, at least eight rows have `Q=9-P` and mass
`p+13`. Proposition 15.769's local theorem excludes them for every
`p=3 (mod 4), p>=31`.

It remains to exclude the new `u=m-2` residue. Every hard quotient is one,
and relative to either parity baseline the difference is a nonzero
nonnegative integral quadratic of scaled mass `p-1`. The height-at-least-two
floor is `p+1`, so this quadratic is Boolean. Its corrected Johnson bound is

```text
L <= 2(p-1)^2(p-2)(3p+1)/(p^2(p+1)(p-3)) < 6.
```

The positivity gap after clearing denominators is

```text
5p^3-20p^2+p+2
=5x^3+445x^2+13176x+129768,   x=p-31.
```

Thus at most five slice coordinates remain and cube influence leaves at most
four active coordinates. The fixed catalog has no table of density
`(p-1)/(4p)`, which lies strictly between its adjacent densities
`(p-3)/(4p)` and `(p+1)/(4p)`. Therefore

```text
p=3 (mod 4), p>=31:  t=q-1 and k=5p-3 are empty.
```

### 3.1 The exceptional `p=23` endpoint

The mass-`p-1` argument above already works at `p=23`: its lift has mass
22 below the height-at-least-two floor 24, the Johnson bound is
`5929/1058<6`, and density `11/46` is absent from the fixed four-bit catalog.

The carried sharp branch needs the exceptional equality theorem from the
preceding layer.  It has eleven low hard rows of one common family and one
high row.  The same opposite `F4/F5` coefficient test leaves only
`P=4,Q=5,F5`.  All eleven low hard coefficient graphs are
triangle-minus-full-star, so they give eleven distinct roots of the common
degree-four and degree-eight forms.  Since `11>8`, both forms vanish
identically.  The opposite `F5/K5` row would then satisfy both moment
identities, contradicting the already certified zero intersection among all
33,649 five-sets.  Thus

```text
p=23:  t=10 and k=112 are empty for every boundary size.
```

See `NOTE_2026-09-04_P23_SECOND_POST_BAND_MOMENT_CLOSE.md` for the exact
residue, edge, sign, and moment ledger.

## 4. Scope

Proposition 15.770 closes the two displayed infinite one-layer families and
the exceptional `p=23,t=10,k=112` endpoint. It makes no claim for later
layers, for the remaining exceptional small-prime regimes, for global
residual (ii), for E1, or for the original quadratic min-max limit.

## Replay

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_prop15770.py tests/test_p23_second_post_band_moment_close.py

PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python src/e1_gmin_m4_prop15770.py
```
