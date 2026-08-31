# The first three residual shells and the `t=3` barrier

**Status.** Proposition 15.735 is a proved theorem: for every prime
`p>=13`, every boundary size, and `t in {1,2}`, residual (ii) has no witness
at

```text
k=4p+2t.
```

Together with Proposition 15.734 (`t=0`), this closes exactly
`k=4p,4p+2,4p+4` in this range. The following `t=3` audit is an **open
reduction**, not part of Proposition 15.735 and not a proof at `k=4p+6`.

## 1. Isolated chart and the common residue

Let the odd flip graph have

```text
h=|H|=4p+2t+1,       q=(p-1)/2,       m=q+1=(p+1)/2.
```

For `p>=13` and `0<=t<=3`,

```text
p^2+1-2h=p^2-8p-4t-1>0.
```

Thus an isolated, hence nonboundary, vertex can be sent to infinity by the
signed PSL transport of 15.721. In the transported chart `I=0`, the boundary
is all finite and even, and every directional odd-fibre count `b_d` is even.
This uses no boundary-size hypothesis.

The phase-one type changes sign when `t` changes parity, but still contains
exactly `m` directions. Proposition 15.632 gives either type the budget

```text
sum a_d=(p+1)(h-3p)/2=2m(m+t).
```

For the phase-one type, write its common residue as

```text
a_d=2u+(p+1)k_d,       0<=u<m,       sum k_d=m+t-u.       (1)
```

For `p>=17` and `t<=3`, the exact floors and the 15.688 lift cost classify
(1) as follows.

* At `u=0`, only `p=1 mod 4`, `b=p-1`, `a=p+1` survives.
* For `1<=u<=t`, every direction requires `k_d>=1`. The only cells below
  `p+1+2u` are lifts of the `b=2` or `b=p-1` equality baselines. Their
  positive lift is at most eight, strictly below `p-3`, so 15.688 excludes
  them.
* For `t<u<=m-2`, every direction again requires `k_d>=1`, but (1) has
  `sum k_d<m`.
* At `u=m-1`, the `b=2`, `a=p-1` baseline survives; when `p=3 mod 4`, the
  `b=p-1`, `a=p-1` baseline also survives. They cannot mix: equal means
  force equal parallel counts, while their coefficient offsets are `4` and
  `3` modulo `q`.

At `p=13`, the exact phase-one floors for `b=0,2,...,12` are

```text
26, 12, 26, 24, 26, 20, 14.
```

They give the same classification for `t<=2`. The entry `b=10,a=20`
creates one additional branch only at `t=3`; see Section 5.

## 2. The three endpoint branches

The generic hard-type possibilities are therefore

| branch | congruence | low baseline | residue `u` | offset `C` | hard finite edges |
|---|---|---|---:|---:|---:|
| A | all `p` | `b=2,a=p-1` | `m-1` | 4 | `mP+t+1` |
| B | `p=1 mod 4` | `b=p-1,a=p+1` | 0 | 5 | `mP+t` |
| C | `p=3 mod 4` | `b=p-1,a=p-1` | `m-1` | 3 | `mP+t+1` |

If `P` is the common low parallel count, coefficient comparison gives

```text
q | P-C,       rho=(P-C)/q,       s=P+rho.
```

The edge bound forces

```text
(P,rho,s)=(4,0,4), (5,0,5), (3,0,3)                 (2)
```

in A, B, C respectively. This also holds in branch B at `t=3`: although
the crude bound permits `P<=8`, the next value congruent to `5 mod q` is
still unavailable for every applicable `p>=13`.

For A and B, an opposite direction with parallel count `Q` has

```text
a(Q)=(p-1)s+(p+1)Q+9-7p,
sum Q=q(8-s)+t.                                      (3)
```

For C it has

```text
a(Q)=(p-1)s+(p+1)Q+7-7p,
sum Q=q(8-s)+t+1.                                    (4)
```

Nonnegativity gives `Q_min=3,2,4`; at those minima the scaled means are
`8,6,8`. The surplus above `m Q_min` is

| `t` | A | B | C |
|---:|---:|---:|---:|
| 1 | `q-2` | `q-1` | `q-2` |
| 2 | `q-1` | `q` | `q-1` |
| 3 | `q` | `q+1=m` | `q` |

## 3. Proposition 15.735: the proved `t<=2` close

For `t=1,2`, every entry in the last table is strictly below `m`. Some
opposite direction therefore attains `Q_min` and has scaled mean `8`, `6`,
or `8`.

That direction has phase zero. If `b_d` is nonempty, its exact floor is at
least `p-1>=12`. If `b_d=0`, parity gives `A_d=2B_d` for a nonzero
nonnegative integral quadratic, and Proposition 15.688 gives

```text
a_d=4p E[B_d]>=p-3>=10.
```

Both alternatives contradict the forced means. This proves Proposition
15.735 for all boundary sizes. Combined with 15.734, the proved uniform band
is

```text
k in {4p,4p+2,4p+4},       p>=13.
```

## 4. `t=3`: the first generic arithmetic survivor

At `t=3`, A and C still have surplus `q<m` and are excluded. Branch B has
surplus exactly `m`. If no forbidden `Q=2` direction occurs, equality forces
every opposite direction to have

```text
Q=3,       a=p+7.                                    (5)
```

The hard side consists of `m` `b=p-1,a=p+1` baselines with three quotient
elevations. Its finite-edge total is `5m+3`; the opposite total is `3m`.
This is the exact point where the preceding pigeonhole/halving argument
stops.

For `p>=17`, the phase-zero floors and 15.688 force every direction in (5)
to have `b=0`, hence `A=2B` with

```text
4p E[B]=p+7.                                          (6)
```

This mean is locally attainable. On `J(17,9)`, fix a four-set `R`, put
`r=|X intersect R|`, and take

```text
B(X)=binom(r,2)-2r+3.
```

Its values for `r=0,...,4` are `3,1,0,0,1`, and `E[B]=6/17`, so
`4p E[B]=24=p+7`.

**This is only a nonnegative integral quadratic satisfying the local mean
condition. It is not asserted to satisfy the full directional coefficient
normalization and is not a residual graph.** It is a counterexample to any
attempt to finish branch B using only the one-direction lift floor or another
mean-halving step.

## 5. The exceptional `p=13,t=3` hard branch

At `p=13`, `m=7`, and `u=3`, equation (1) permits `k_d=1` in all seven hard
directions. The exact phase-one cell `b=10,a=20` realizes the whole budget:

```text
7*20=140=2m(m+3).
```

It has the baseline `(2-r)^2` on the three-point complement, coefficient
target `P=5`, and hard finite-edge total `35`. The opposite ledger has

```text
sum Q=24,       a(Q)=14(Q-2).
```

The available floors do not exclude this allocation. This is an additional
arithmetic survivor, not a constructed residual graph.

For reference, the `p=13` phase-zero floors at `b=0,2,...,12` are

```text
0, 14, 20, 26, 24, 26, 12.
```

Thus the generic branch-B mean `a=20` can also use the exact `b=4` cell;
unlike `p>=17`, `b=0` is not forced.

## 6. A coefficient-normalized necessary-condition witness

For a `b=0,Q=3` direction in generic branch B, the shared slice
coefficients imply an integral signed edge matrix `W=(epsilon K)` on the
`p` coordinates such that

```text
sum_(s<t) W_st=-(p+7),
every row sum is even,
sum_(s<t) |W_st| <= h-Q=4p+4,
B(X)=-(p+7)/4 - (1/2) sum_(s in X,t notin X) W_st >=0
```

for every `m`-set `X`. Equivalently, every balanced cut has `W`-weight at
most `-(p+7)/2`.

This stronger single-direction condition still does not uniformly exclude
`t=3`. At `p=13`, take

```text
W=-Adj(K_{2,10} disjoint-union K_1).
```

Its row sums are even, `sum W=-20`, and `sum |W|=20<=56`. If a seven-set
contains `alpha` of the two centres, `beta` of the ten leaves, and `gamma`
of the isolated vertex, then

```text
alpha+beta+gamma=7,
cut=10 alpha+2 beta-2 alpha beta.
```

For `alpha=0,1,2`, the possible cut sizes are respectively `12 or 14`,
`10`, and `10 or 12`. Hence

```text
B(X)=-5+cut(X)/2
```

is integral and nonnegative, with `E[B]=5/13` and `4p E[B]=20=p+7`.

**This is an exact witness only for the displayed single-direction
coefficient and cut necessary conditions. It does not construct the other
directions, a common flip graph, either separator inequality, or a residual
witness.**

## 7. Exact scope and next gate

Proposition 15.735 proves `t=1,2`; it makes no claim at `t=3`. The audit
above leaves at `k=4p+6`:

1. generic branch B for `p=1 mod 4`, with every opposite direction at
   `Q=3,a=p+7`; and
2. the extra `p=13,u=3,b=10,a=20` hard branch.

The local and coefficient-normalized examples show that another independent
direction floor, moment bound, or averaging/halving argument cannot close
the whole shell. A further proof must couple several directional matrices,
or couple the hard and opposite types, through the common flip graph and
line-incidence reconstruction.
