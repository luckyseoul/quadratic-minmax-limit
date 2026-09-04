# Proposition 15.768: first `p=1 mod 4` layer beyond 15.752

**Status:** proved infinite-family theorem with the fixed four-bit catalog of
Proposition 15.751.  No prime, graph, orbit, slice, or residual-cell census is
used.

## Statement

For every prime `p>=29` with `p=1 mod 4`, residual (ii) is empty at

```text
q=(p-1)/2,  m=q+1,  t=q-3=(p-7)/2,
k=4p+2t=5p-7,  |E(H)|=5p-6.
```

This is exactly the first layer beyond Proposition 15.752's
`p=1 mod 4` interval, which ends at `t=q-4`.  The theorem does not treat the
first uncovered `p=3 mod 4` layer, any later layer, primes below `29`, or
residual (ii) globally.

## 1. Exact new hard residue

An isolated projective vertex exists because

```text
p^2+1-2(5p-6)=p^2-10p+13>0.
```

Signed transport sends it to infinity, so `I=0`, the boundary is all finite,
and every directional odd-fibre count is even.  The phase-one type has budget

```text
2m(m+t)=m(2p-6).
```

Write

```text
a_L=2u+(p+1)k_L,       sum_L k_L=m+t-u.
```

The same exact floors and sharp `p-3` lift obstruction used by Proposition
15.752 give the following exhaustive residue list.

- At `u=0`, the old exact `b=p-1` literal branch survives.
- For `1<=u<=t-1`, every forced `k_L=1` cell is a positive lift strictly
  below `p-3`.
- At `u=t`, `sum k_L=m`, so every `k_L=1`; the new exact floor is
  `b=p-3`, `a_L=2p-6`.
- For `t+1<=u<=m-2`, a `k_L=0` cell is forced below the phase-one floor.
- At `u=m-1`, the old exact `b=2` XNOR branch survives.

Thus the only branches are the old A/B branches and the new complement-triple
branch.

## 2. The complement-triple equality is pointwise

Let `C` be the three-point complement of the `p-3` odd fibres and put
`r=|X intersect C|`.  Since `m` is odd, phase one on the large set becomes
phase zero in `r`.  Proposition 15.652's exact quadrature has candidate

```text
A_0=(2-r)^2,
```

contact layers `r=1,2,3`, and strictly positive weights.  Equality therefore
forces the original quadratic `A`, not only its stabilizer average, to agree
with `A_0` at every point in those three layers.

The omitted `r=0` layer follows without a finite rank computation.  Fix such
an `m`-set, choose three of its outside points, and swap them independently
with the three points of `C`.  The difference `A-A_0` restricts to a
degree-at-most-two function on this three-cube and vanishes at all seven
nonzero vertices.  Its third finite difference is zero, so its value at the
origin also vanishes.  Hence equality is pointwise everywhere.

In signed coordinates,

```text
3+2A_0
 =5-sum_(i in C) z_i+sum_({i,j} subset C) z_i*z_j.       (1)
```

The coefficient offset in (1) is `5+(-1-1-1)=5-3=2`.  More precisely, the
slice-ideal coefficient comparison is

```text
I+P_L-offset=(p-1)c,       with 2c integral.
```

Consequently the slice-kernel comparison gives

```text
q divides I+P_L-2,       hence P_L=2 mod q.               (2)
```

## 3. The common-row normalization

For every normalized direction row, the local coefficient sum and the one
common graph give, respectively,

```text
sum_a q_L(a)=p(P_L-3)-a_L,
sum_a q_L(a)=hT-P_L                                      (3)
```

on the hard type.  Equating the two sides of (3) is essential: this is not a
collection of independent local cells.  Since all new hard means equal
`2p-6`, (3) forces every hard parallel count to be the same.  Moreover,

```text
P_L<=floor((5p-6)/m)<10<q+2.
```

Together with (2), this gives

```text
P_L=2,
hard edges=mP_L=p+1,
opposite edges=4p-7,
hT=(p+1)-(4p-7)=8-3p.                                  (4)
```

For an opposite direction of parallel count `Q`, the second sign in the
common-row identity gives

```text
a(Q)=(p+1)Q-6p+8,       sum Q=4p-7.                    (5)
```

Now `a(5)=13-p<0` and `a(6)=14`.  The latter is below both the least nonzero
phase-zero odd-fibre floor and the sharp `p-3` empty-fibre lift floor.  Hence
every `Q>=7`.  But

```text
sum(Q-7)=(p-21)/2<m,
```

so some direction has `Q=7` and

```text
a(7)=p+15.                                              (6)
```

At the mean in (6), the only nonzero even-`b` floors that fit are `b=2` and
`b=p-1`.  Their excesses over the pointwise equality baselines are `14` and
`16`; both lie strictly between zero and `p-3`.  Proposition 15.688 excludes
them.  Thus `b=0`, `A=2B`, and the new branch requires

```text
B:J(p,m)->Z_{>=0},  B nonzero,  4p E[B]=p+15.          (7)
```

## 4. A sharp dimension-free three-quarter-mean cube theorem

We need one endpoint strengthening of Proposition 15.751's half-mean cube
theorem:

> If `g` is a nonnegative integer-valued multilinear polynomial of degree at
> most two on a Boolean cube and `E g=3/4`, then `max g<=6`.

Take a counterexample of minimum dimension, move a maximizing point to the
origin, and write its value as `M>=7`.  Dimensions at most three are already
impossible because the total mass is below seven, so every coordinate facet
has dimension at least three.  Integer cube values give integral multilinear
coefficients by finite differences; summing a degree-at-most-two polynomial
on such a facet therefore shows that its mean lies in `(1/4)Z`.
For every coordinate, the means of the facets through and opposite the origin
sum to `3/2`.  The possible positive through-origin means are therefore

```text
1/4, 1/2, 3/4, 1, 5/4, 3/2.
```

The first three are impossible: support equality bounds the maximum at mean
`1/4` by one, Proposition 15.751 bounds it at mean `1/2` by three, and mean
`3/4` contradicts minimum dimension.  At mean `3/2`, the opposite facet has
mean zero and vanishes.  Thus `g=(1-x_i)h`; degree at most two makes `h`
affine, and `h(x)+h(1-x)=3` bounds its maximum by three.  The only remaining
through-origin means are `1` and `5/4`.  Their opposite means are `1/2` and
`1/4`, so every opposite facet has maximum at most three.  Every nonorigin
vertex lies on one of those facets and consequently has value at most three.

If the dimension is at least five, fix any five-coordinate subcube through
the origin and average over coordinate permutations.  Its layer averages
`q(s)` form a quadratic in `s`.  Exact Lagrange interpolation at `s=1,3,5`
gives

```text
M=q(0)=15q(1)/8-5q(3)/4+3q(5)/8 <= 27/4 < 7.          (C)
```

Here `q(1),q(5)<=3` and `q(3)>=0`.  Integrality gives `M<=6`.  In dimensions
at most three the entire mass is below seven.  In dimension four, the
vanishing fourth alternating difference equates the even- and odd-parity
masses; their total is `12`, so each is six and again `M<=6`.

The theorem is sharp.  On six variables,

```text
g=6-3s+binom(s,2),
```

has layer values `6,3,1,0,0,1,3`, total mass `48`, mean `3/4`, and maximum
six.

## 5. Uniform local `p+15` theorem

Let `H=max B`.  If `H>=2`, every paired cube through a maximizing point has
quarter-integral mean at least `1/2`.  The paired-cube identity gives

```text
H >= (p-13)/4.                                          (8)
```

For `p=1 mod 4`, the sharp stabilizer bound gives

```text
H <= (p+15)(p+3)/(4(p-1)),
T B(X) <= (p+15)/(2(p-1)) < 3/4.                        (9)
```

The last inequality holds for `p>=37`.  Therefore one paired cube has mean
exactly `1/2`; Proposition 15.751's dimension-free half-mean theorem gives
`H<=3`, contradicting `(p-13)/4>3`.

The endpoint `p=29` closes by a two-stage exact bootstrap.  First the
paired-cube lower bound (8) gives `H>=4`, so a cube containing the maximizer cannot
have mean `1/2`; the half-mean theorem would give `H<=3`.  Hence every paired
cube has mean at least `3/4`.  Since

```text
E B=11/29,  rho=1/30,
T B(X)=rho H+(1-rho)E B=(H+11)/30,
```

this forces `H>=23/2`, hence `H>=12`.  The stabilizer bound is `H<=88/7`,
hence integrally `H<=12`.  Thus `H=12` and `T B(X)=23/30<1`.  Quarter
integrality now forces one paired cube to have mean exactly `3/4`, while the
theorem of Section 4 bounds its maximum by six.  It contains the maximizing
point of height twelve, a contradiction.

It remains that `B` is Boolean with density

```text
mu=(p+15)/(4p).
```

The corrected Johnson transposition-influence calculation gives, for the
complement size `L` of the largest zero-influence coordinate class,

```text
L <= 2(p-1)(p-2)(p+15)(3p-15) /
       (p^2(p+1)(p-3)) < 8.                             (10)
```

After writing `p=x+29`, the cleared numerator for the strict inequality in
(10), up to a positive factor two, is

```text
x^4+87x^3+2820x^2+40880x+228912>0.
```

Thus `L<=7<(p-1)/2`.  This weaker-than-15.752 bound is sufficient: every
pattern on those seven coordinates extends to the complementary middle
slice.  Symmetrization gives a Boolean degree-at-most-two cube function, and
cube influence reduces it to at most four active coordinates.  Proposition
15.751's fixed four-bit catalog has density list

```text
0, 1, (p-3)/(4p), (p+1)/(4p), (p-1)/(2p), (p+1)/(2p),
(3p-1)/(4p), 3(p+1)/(4p).
```

For `p>=29`, `(p+15)/(4p)` lies strictly between `(p+1)/(4p)` and
`(p-1)/(2p)`, so it is absent.  This proves the local theorem excluding (7).

## 6. The old branches still close

At `t=q-3`, the old XNOR branch still forces an opposite `Q=4` cell of mass
`p+9`; the surplus after raising every opposite count to four is
`(p-15)/2<m`.  Proposition 15.752 excludes that cell.

The old literal branch still forces an opposite `Q=3` cell of mass `p+7`;
the remaining surplus is `(p-13)/2<m`.  Proposition 15.751 excludes that
cell.  Section 5 excludes the sole new branch, so the three-branch residue
ledger is empty.

Therefore

```text
boxed: k=5p-7 is impossible for every prime p=1 mod 4, p>=29.
```

## Replay

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python src/e1_gmin_m4_prop15768.py
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 tests/test_prop15768.py
```
