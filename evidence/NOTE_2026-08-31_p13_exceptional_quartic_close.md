# The corrected p=13 exceptional branch and the generic moment frontier

**Status.**  Proposition 15.739 is a proved branch theorem: the exceptional
`p=13,t=3,u=3` row at `k=58`, `|H|=59` is empty.  It does **not** close the
generic `p=13` row, the generic `p>=17` fourth shell, residual (ii),
multi-level Type I, or the quadratic min-max limit.

## 1. The coefficient correction

Let `C` be the three-point complement in an exceptional hard direction and
put `r=|X intersect C|`.  The exact phase-one baseline is

```text
A=(2-r)^2.
```

With `z_i=2x_i-1`, direct expansion gives

```text
epsilon*S_H
 = 3+2A
 = 5-sum_(i in C) z_i+sum_({i,j} subset C) z_i*z_j.
```

The coefficient offset is therefore `5-3=2`.  The previous value five in
the first-three-shell note was a normalization error.  If `P` is the common
hard parallel count, polarization of the slice kernel gives

```text
6 | P-2.
```

Seven hard directions use `7P` of the 59 finite edges, so `P<=8`.  Hence
only `P=2,8` remain.

The common directional mean formula is

```text
20=14P-h*T-39,
```

where `h` is the hard direction sign.  An opposite direction of parallel
count `Q` therefore has

```text
a(Q)=14Q+h*T-39=14(P+Q-7).
```

Exact finite-edge accounting gives

| hard `P` | opposite `sum Q` | opposite mean |
|---:|---:|---:|
| 2 | 45 | `14(Q-5)` |
| 8 | 3 | `14(Q+1)` |

For `P=2`, the nonnegative minimum `Q=5,a=0` would have target
`epsilon*S_H=3`, offset three.  It is incompatible with `Q=5` modulo six.
Thus every opposite count is at least six; the total 45 forces at least four
`Q=6,a=14` directions.  For `P=8`, the total three forces at least four
`Q=0,a=14` directions.

## 2. Importing the exact mass-fourteen cell catalog

The phase-zero floors at `b=0,2,...,12` are

```text
0, 14, 20, 26, 24, 26, 12.
```

At mean fourteen, the exact `b=2` target is `4-z_i*z_j`, with offset four,
so it is incompatible with both `Q=0` and `Q=6`.  A `b=12` cell would be a
two-unit lift of its pointwise parity baseline; Proposition 15.688 requires
any nonzero integral lift to cost at least ten.  Every other nonzero `b` has
floor above fourteen.  Therefore

```text
b=0,       A=2B,       4p*E[B]=14.
```

Proposition 15.738 first uses the two inequalities in Proposition 15.688 to
reduce `max B` to one or four.  Its exact `Q=0` and `Q=6`
coefficient/cut/row-parity/l1 models exclude height four.  It then certifies
rank 78 for the pair-evaluation space on `J(13,7)`, constructs a rank-1638
third-difference annihilator, and proves by an anchored exact CP-SAT no-good
model that the following 1,092 Boolean supports exhaust support 462:

| family | count | coefficient offset |
|---|---:|---:|
| `x_i*x_j` | 78 | 6 |
| `x_i*(1-x_j)` | 156 | 4 |
| `z_i=z_j=-z_k` | 858 | 4 |

Modulo six, only `B=x_i*x_j` survives at `Q=0,6`.

## 3. The sign-safe quartic

For even `d`, keep the global homogeneous binary form

```text
M_d(L)=sum_({u,v} in H) chi(u-v)*(L(u)-L(v))^d.
```

Do not multiply this form by the direction-dependent sign before applying a
root count.  Instead fix `h`, the common sign on all seven hard directions,
and use the local normalized matrix `W=hK` only to compute its values.

For either hard gauge `P=2,8`, the even normalized moments equal those of
the complement triangle.  After normalizing its labels to `{0,1,r}` and
putting `q_0=r^2-r+1`,

```text
S_2=2q_0,       S_4=2q_0^2,       2S_4=S_2^2.
```

Consequently the genuine homogeneous quartic

```text
G(L)=2h*M_4(L)-M_2(L)^2
```

vanishes in all seven distinct hard projective directions.  A nonzero
binary quartic has at most four projective roots, so `G` is identically
zero.

The opposite selected-pair cell has sign `-h`.  At either `Q=0` or `Q=6`,
its normalized moments are gauge-independent:

```text
S_2=(i-j)^2,       S_4=(i-j)^4.
```

Evaluation in the global quartic gives

```text
G=-3(i-j)^4 != 0 mod 13,
```

contradicting `G=0`.  This closes the exceptional row.

## 4. What remains in the generic branch

For generic branch B at `p=4a+1>=17`, the hard quotient identity has only
three excess units.  Hence at least

```text
N=m-3=(p-5)/2
```

hard directions are exact stars.  For every orientation-independent even
degree `d<N`, the star power sum is zero and the degree-`d` binary form has
more roots than its degree.  Thus

```text
M_d=0,       d=2,4,...,(p-9)/2.
```

For an opposite normalized matrix this adds the congruences

```text
sum W_st*(s-t)^d=0 mod p
```

to the exact conditions

```text
sum W=-(p+7),
even row sums,
sum |W|<=4p+4,
cut_W(X)<=-(p+7)/2 for every middle-slice X.
```

Conditional cut averages force every entry into `{-1,0,1,2,3}`.  At
`p=17`, they additionally give even row degrees in `[-16,16]` and the exact
pair, triple, and four-set inequalities recorded by Proposition 15.739.
The source recomputes all of these from exact hypergeometric conditional
averages rather than treating the displayed inequalities as trusted input.

There is also a two-intersection stabilizer consequence.  Fix any nine-set
`X` and average `W` over its stabilizer, writing `a,b,c` for the coefficients
on the 36 inside, 28 outside, and 72 crossing edges.  The total sum and the
cut constraints averaged over nine-sets meeting `X` in four and five points
are

```text
36a+28b+72c=-24,
20a+15b+37c<=-12,
20a+16b+36c<=-12.
```

The exact combination `-9*(four-intersection inequality)
-(45/4)*(five-intersection inequality)+(45/4)*(total equality)` gives
`72c>=-27`.  Every cut is even because every row degree is even, so

```text
-26 <= cut_W(X) <= -12.
```

Consequently `B(X)=-6-cut_W(X)/2` takes values only in `{0,...,7}`, has
mean `6/17`, and has total mass 8,580 on `J(17,9)`.
The complete-domain solver route was terminated after a bounded mesh audit.
Every CP-SAT run materialized all 24,310 cuts and both moments; none returned
an incumbent or an infeasibility certificate:

| encoding | host | limit / solve time | status | exported model SHA-256 |
|---|---|---:|---|---|
| direct cuts, pair/triple/four averages | NUKA | 1200.402 s | `UNKNOWN` | `a492c3eb1b4266d365b63eabf9968436ee08cac36861c9efe0582f7df14c77f4` |
| half cuts, `l1` objective | ORIN | 1200.301 s | `UNKNOWN` | `d5248bc94d32afccf9daa6cf8b4d22baec5b6bf3d14a7cb10f2a78712645d1ed` |
| split signs, conditioned through size five | NUKA | 700.012 s | `UNKNOWN` | `b0b079a25ef6720017ba6acf10283f104b44a8393eac1f668d7106283d4c6969` |
| split signs, conditioned through size six | NUKA | 1200.116 s | `UNKNOWN` | `80deae7c562c788c7708dfe6fd89e0cbf1d5e7e035a5ca325f309000237ba789` |
| split signs, conditioned through size five | ORIN | 1200.094 s | `UNKNOWN` | `20887f2c7c2c3452614b6dab8aa47209f31aa253467323502cf002001d558392` |

A lazy exact-cut model returned `UNKNOWN` before finding its first incumbent.
Five independent Z3 variants never returned from `check()`: three reached
external guards and two were explicitly stopped.  Their result is
`NO_SOLVER_STATUS`, not Z3 `unknown` and not evidence of infeasibility.  The
temporary run archive and full hash manifest are at
`/tmp/qml-t3-p17-20260831T055835Z/` for this session.  This backend/seed/
timeout campaign is closed.  A future p17 attack must prove new structure for
the `{0,...,7}`-valued quadratic, not rerun the same exact model.

At `p=13`, the generic branch forces only `M_2=0`.  There is an exact local
elevated hard cell showing why the one-direction extension stops.  For
`P=6`, take

```text
W=Adj(K_5 on {0,1,2,3,5})+1_{{0,11}}.
```

It has `sum W=l1(W)=11`, odd rows `{0,11}`, and every seven-cut is between
zero and seven.  Therefore `A(X)=7-cut_W(X)` is nonnegative, has scaled mean
28 and `b=2`.  Its moments are

```text
S_2=0,       S_4=5 mod 13,
M_2=h*S_2=0, M_4=h*S_4=5h.
```

This is a counterexample to forcing an elevated hard direction back to the
`b=12` star using its floor and the quadratic moment.  It is not a common
flip graph or a residual witness.

Proposition 15.740 subsequently excludes the generic hard partitions with
five or six exact stars by cyclic translation averages.  The only surviving
p13 partition is `1^4 2^3`.  Its next live implication is precise: exclude
the simultaneous completion of four exact stars, three elevated hard cells,
and seven opposite cells by one common 59-edge graph.  Another independent
direction floor, local catalog, or mean-halving pass cannot do it.
