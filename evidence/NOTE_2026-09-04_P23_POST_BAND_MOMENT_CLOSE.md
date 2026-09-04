# The p=23 first post-band endpoint closes by equality rigidity and two moments

Date: 2026-09-04

Scope: the exceptional `p=23,t=9,k=110` first layer beyond Proposition
15.752's `p=3 mod 4` band.  This is a proved endpoint theorem with a fixed
finite coefficient certificate.  It is not a graph census, a later-layer
theorem, or a global closure of residual (ii).

## 1. Exact isolated-chart ledger

Put

```text
p=23, q=11, m=12, t=9, k=110, |H|=111.
```

There are

```text
p^2+1-2|H|=308
```

guaranteed isolated vertices, so the signed transport used in Proposition
15.752 gives an all-finite chart with `I=0`.  The phase-one means and quotient
sum are

```text
a_L=2u+24 k_L,             sum_L k_L=21-u.
```

The same exact floor ledger leaves only `u=9` and `u=11`.  At `u=11`, the
two old hard branches force an opposite cell of mass `p+9=32`; Proposition
15.752 excludes both.

At `u=9`, all twelve hard quotients equal one and every hard cell has mean
`2p-4=42`.  Relative to its phase-one parity baseline, its difference is a
sharp nonnegative integral lift of scaled mass `p-3=20`.  Proposition 15.688
makes that lift Boolean.  The corrected transposition bound is

```text
L <= 6(p-1)(p-2)/p^2 = 2772/529 < 6.
```

Thus at most five slice coordinates remain; cube influence reduces to four
active coordinates.  The fixed 65,536-table catalog has exactly ten tables
at density `5/23`: six selected-pair tables and four all-equal-triple tables.
Complementing back to `J(23,12)` gives an omitted pair or an all-equal
triple.

Combining these two lifts with the XNOR and complementary-literal baselines
gives offsets `2,3,4,5`.  Equal hard means and the common unspecialized row
sum force one common hard parallel count `P`.  The coefficient congruence

```text
11 divides P-offset
```

and `P<=floor(111/12)=9` force `P=offset`; different families cannot mix.
For each family

```text
hT=24P-111,
a(Q)=24Q+hT-69,
sum opposite Q=111-12P.
```

The row `Q=8-P` has mass twelve, below both the nonzero phase-zero floor 24
and the sharp lift floor 20.  Hence every opposite row has `Q>=9-P`.
After assigning that count to all twelve directions, the surplus is exactly
three.  At least nine directions therefore have `Q=9-P` and mean 36.  Every
nonzero-boundary alternative is a phase-zero baseline of mass 24 plus a
nonzero lift of mass twelve, again below the sharp lift floor.  Consequently
each of those exact low rows has `b=0` and

```text
A=2C,        C>=0 integral quadratic on J(23,12),
92 E[C]=36.
```

## 2. The mass-36 cell has height exactly three

Let `H=max C`, attained at `X`.  Height one is impossible by the corrected
Johnson/cube reduction: the exact largest-zero-class complement bound is

```text
19404/2645 < 8,
```

so at most seven slice coordinates remain and then at most four cube
coordinates remain.  The fixed four-bit density list is

```text
0, 5/23, 6/23, 11/23, 12/23, 17/23, 18/23, 1,
```

which misses the required Boolean density `9/23`.

For `H>=2`, every paired cube through `X` is a nonzero integral quadratic.
Its mean is quarter-integral.  Mean `1/4` would attain the degree-two support
floor and hence be Boolean, impossible because the cube contains the value
`H>=2`.  Thus every paired cube has mean at least `1/2`.  The exact paired
operator is

```text
T C(X)=(C(X)+23 E[C])/24=(H+9)/24,
```

so initially `H>=3`.  The `p=3 mod 4` stabilizer inequality gives `H<=9`.
If `H>=4`, no paired cube may have mean `1/2`, since Proposition 15.751's
half-mean theorem bounds its maximum by three.  Every cube would then have
mean at least `3/4`, forcing `H>=9`.  Hence `H=9`, their average is exactly
`3/4`, and every paired cube has mean `3/4`.  Proposition 15.768's
dimension-free three-quarter theorem bounds its maximum by six, contradicting
the value nine at `X`.  Therefore

```text
H=3,
```

and the paired-cube average is exactly `1/2`.  Since every paired cube has
mean at least `1/2`, every one of them has mean exactly `1/2`.

## 3. All half-mean height-three cube quadratics

The following equality classification is dimension-free.

> If `g` is a nonnegative integer-valued polynomial of degree at most two on
> a Boolean cube, `E[g]=1/2`, and `g(0)=max g=3`, then, up to permutation and
> dummy coordinates,
>
> ```text
> g=F_r(s):=3-2s+binom(s,2),       r in {4,5},
> ```
>
> where `s` is the sum of the `r` active coordinates.

Here is a complete proof.  For any coordinate, write `h=g|x_i=0` and
`k=g|x_i=1`.  Their means sum to one and lie in `(1/4)Z`.  Since `h(0)=3`,
its mean cannot be `1/4`: support-floor equality would make every nonzero
value one.  Its mean cannot be one either.  In that case `k=0`, the affine
difference `k-h=-h` would make `h` affine, and writing
`h=3+sum c_j x_j` would give `sum c_j=-4` from the mean while origin
maximality and nonnegativity at the all-one point give `sum c_j>=-3`.
Thus the origin facet has mean `1/2` or `3/4`.

For dimensions zero, one, and two the total cube mass is respectively
incompatible with integrality or strictly below the value three at the
origin, so these are the induction bases.  Induct on dimension.  Suppose
some origin facet has mean `1/2`.  By
induction it is `F_4` or `F_5` on active coordinates `R`, with possible
dummy coordinates `D`.  The opposite facet also has mean `1/2`, and their
difference is integral affine:

```text
ell=c+sum_(r in R)u_r x_r+sum_(z in D)v_z x_z.
```

For an `F_4` facet, `h=0` on every active two-set.  Average `ell` over the
dummies and write `c'=c+(1/2)sum_D v_z`, `U=sum_R u_r`.  The six numbers

```text
q_rs=c'+u_r+u_s
```

are nonnegative.  Their sum is `6c'+3U=0`, because `E[ell]=c'+U/2=0`.
Thus every `q_rs=0`.  Pointwise `ell` on an active two-set is the
nonnegative opposite-facet value and has dummy average zero, so it vanishes
for every dummy assignment.  Hence all `v_z=0` and
`c+u_r+u_s=0` for every pair.  All `u_r=u`, `c=-2u`, and the empty active
vertex gives `-3<=c<=0`.  Integrality leaves only

```text
(c,u)=(0,0)          or          (c,u)=(-2,1).
```

These are respectively a dummy extension and
`F_4(s)+(-2+s)x_i=F_5(s+x_i)`.

For an `F_5` facet, the empty and full active vertices both have value
three.  Their dummy-averaged differences are `c'` and `c'+U=-c'`; both are
nonpositive because the opposite-facet values are at most three.  Therefore
`c'=0`.  At the empty active vertex the pointwise difference is nonpositive
with dummy average zero, so it vanishes identically: `c=0` and every
`v_z=0`.  Now `U=0`.  On every active two-set `F_5=0`, so
`u_r+u_s>=0`.  Summing the ten inequalities gives `4U=0`; every pair sum is
zero and all `u_r=0`.  Only a dummy extension remains.

It remains to handle the case where every origin facet has mean `3/4`.
Every opposite facet then has mean `1/4`, so support-floor equality makes it
Boolean.  Write

```text
g=3+sum_i a_i x_i+sum_{i<j} b_ij x_i x_j,
A=sum_i a_i, B=sum_{i<j}b_ij.
```

Booleanity at the opposite-facet origins gives `a_i in {-3,-2}`.  The
global mean and each facet-mean difference give

```text
2A+B=-10,
2a_i+sum_(j!=i)b_ij=-1.
```

Summing the second equations yields `A+B=-d/2`, hence
`A=-10+d/2`.  In particular `d` is even.  Since `A<=-2d`, one has `d<=4`,
while the total mass `2^(d-1)` is at least the value three at the origin, so
`d>=3`.  Thus `d=4`, `A=-8`, and every `a_i=-2`.  Each row of `b` has sum
three.  Nonnegativity on a two-set gives `b_ij>=1`; each row has three such
entries, so every `b_ij=1`.  This is `F_4`.  The induction and its base are
complete.

The two forms have layer values

```text
F4: 3,1,0,0,1,
F5: 3,1,0,0,1,3,
```

and both have cube mean `1/2`.

## 4. Compatibility between all pairings globalizes the slice form

Fix the maximizing twelve-set `X`.  For `x in X` and `y outside X`, define

```text
D_xy=1  iff  C(X-x+y)=1.
```

Every cross edge occurs in a paired cube.  The equality classification says
that its one-flip value is either one (active) or three (dummy), so `D` is a
well-defined binary `12 by 11` matrix.

Take distinct `x,x' in X` and `y,y' outside X`.  The same double swap

```text
X-x-x'+y+y'
```

appears once in a paired cube containing `(x,y),(x',y')` and once in a
paired cube containing `(x,y'),(x',y)`.  Both `F_4` and `F_5` take the
distinct values `3,1,0` when zero, one, or two of these coordinates are
active.  Equality of the double-swap value therefore gives

```text
D_xy+D_x'y' = D_xy'+D_x'y.
```

All additive `2 by 2` minors vanish.  Hence `D_xy=r_x+c_y`; because `D` is
binary, either all rows are identical (column-only) or every row is constant
(row-only).

Every near-perfect matching uses all eleven columns and all but one of the
twelve rows.  Its active-edge count is four or five.  In the column-only
case this leaves exactly four or five active columns.  In the row-only case,
omitting arbitrary rows leaves both `r` and `r-1` active edges, forcing
exactly five active rows.

Every `Z in J(23,12)` lies in a paired cube through `X`: retain any point of
`Z intersect X`, match the swapped points, and extend to a bijection.  Thus
the column cases globalize to `F_4(|Z intersect R|)` or
`F_5(|Z intersect R|)`.  The row case gives
`F_5(5-|Z intersect R|)=F_5(|Z intersect R|)`.  These are all possible
mass-36 slice quadratics; there are no unclassified nonsymmetric equality
forms.

For `|R|=r`, direct conversion with `z_i=2x_i-1` gives

```text
3+4C = 15-4r+binom(r,2)
       +(r-5) sum_(i in R) z_i
       +sum_({i,j} subset R) z_i z_j.
```

The coefficient offset (constant plus the linear coefficients) is one for
`F_4` and five for `F_5`.  In the four hard-family ledgers the exact low
opposite counts are respectively

```text
(P,Q)=(2,7),(3,6),(4,5),(5,4).
```

The opposite congruence is `11 divides Q-offset`.  Its only solution is

```text
P=4, Q=5, C=F5.
```

This is precisely the complementary-literal plus all-equal-triple hard
family.  The four-coordinate offset-one example is real, but it is not
coefficient-compatible with any forced opposite count.

## 5. Twelve hard roots force two global moment identities

In the surviving hard family, if `j` is the literal coordinate and `T` is
the triple, its signed target is

```text
5-z_j+sum_({a,b} subset T) z_a z_b.
```

On the hard slice `sum_i z_i=1`, the identity

```text
sum_(i!=j) z_i z_j=z_j-1
```

rewrites this as parallel constant four plus the coefficient graph

```text
triangle(T) - full_star(j).
```

Here the isolated-chart normalization has `I=0`, and the coefficient
congruence has forced the parallel count to equal the target offset `P=4`.
Thus the slice-kernel scalar is zero: the displayed triangle-minus-star is
the actual normalized coefficient graph, not merely its class modulo the
slice ideal.  On the opposite row the same facts are `I=0` and `Q=5` equal
to the `F5` offset, so its actual normalized coefficient graph is `K5`.

For even `d<22`, every full-star moment vanishes over `F_23`:

```text
sum_(a in F_23^*) a^d=0.
```

For a triangle with difference coordinates `x,y,x-y`, put
`S_d=x^d+y^d+(x-y)^d`.  Direct polynomial identities give

```text
2S4-S2^2=0,
24S8-32S2*S6+5S2^4=0.
```

Let the genuine global homogeneous moment forms be

```text
M_d(L)=sum_({u,v} in H) chi(u-v)(L(u)-L(v))^d.
```

All twelve hard directions have one common sign `h`, so they are roots of

```text
G4(L)=2h M4(L)-M2(L)^2,
G8(L)=24h M8(L)-32M2(L)M6(L)+5M2(L)^4.
```

These are homogeneous binary forms of degrees four and eight.  Twelve
distinct projective roots force both forms to vanish identically over
`F_23`.

On an opposite `F_5` row the sign is `-h` and its normalized coefficient
graph is the complete graph on a five-set `R subset F_23`.  Writing

```text
S_d(R)=sum_({i,j} subset R)(i-j)^d,
```

the two identically zero global forms require

```text
-2S4-S2^2=0,                                      (G4)
-24S8-32S2*S6+5S2^4=0.                           (G8)
```

## 6. Exact five-set certificate

The authoritative deterministic replay checks all

```text
binom(23,5)=33,649
```

five-sets with integer modular arithmetic.  Its results are

| condition | five-sets | SHA-256 of lexicographic set list |
|---|---:|---|
| `G4=0` | 1,518 | `82460f67f3414a1f461b24605c108861d215f970063c0d0af82772de21240c1a` |
| `G8=0` | 2,024 | `733bc62c7ad8d0d7083388480d307ad7298d56b4f9e1fcd12562848350c8d6c7` |
| both | **0** | empty |

An independent orbit replay partitions all five-sets into 69 orbits under
`AGL(1,23)`: 64 have size 506 and five have size 253.  The representative
list with orbit sizes has SHA-256

```text
34eeb59b625d24907758658f78c0f966291728a72cebb0426a3d4a883fb2022a.
```

The three `G4`-zero representatives are

```text
{0,1,2,3,12}, {0,1,2,4,15}, {0,1,2,7,17},
```

and the four `G8`-zero representatives are

```text
{0,1,2,3,10}, {0,1,2,4,17}, {0,1,2,4,18}, {0,1,2,7,10}.
```

All seven have full orbit size 506 and the two lists are disjoint, giving
the same counts `3*506=1518`, `4*506=2024`, and joint count zero.

The independent mesh implementation is
`scripts/p23_k5_moment_gpu.py`, SHA-256
`4afa5d397ccf38dc6e61f8b006c7a697c533ff067b9c7af944441f35c986ee1d`.
Current replays on NUKA gfx1201 OpenCL, Jellyfin Arc A380 OpenCL, and
Soulkiller V100 CUDA each returned the strict-five-set vector

```text
[33649,1518,2024,0].
```

Those accelerator runs are independent implementation checks; the exact
CPU/all-orbit replay is the authoritative certificate.

Since at least nine opposite rows would require a five-set satisfying both
identities and not even one exists, the last branch is impossible.  Therefore

```text
p=23,t=9,k=110 is empty for every boundary size.
```

Later layers, residual (ii) globally, `E(1)`, `L=1/2`, and the original
MathOverflow limit remain open.

## Replay

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_p23_post_band_moment_close.py
```

Machine-readable evidence:
`evidence/e1_gmin_m4_p23_post_band_moment_close.json`.
