# The six-root quartic close at `p=13,t=4,u=3`

**Status.** Proposition 15.744 is a **proved branch theorem**, with a new
exhaustive finite certificate for its only changed local premise. It closes
exactly

```text
p=13,  t=4,  k=60,  hard residue u=3.
```

It does not close the full `p=13,k=60` row. Residues `u=0,4,6` remain.

## 1. The exact hard ledger

The isolated chart has `I=0`, `|H|=61`, seven directions of each sign, and
phase-one type budget

```text
2m(m+t)=2*7*11=154.
```

At residue `u=3`, write `a_L=6+14k_L`. Then

```text
sum k_L=8.
```

The exact phase-one even-`b` floors are

```text
b:       0   2   4   6   8  10  12
floor:  26  12  26  24  26  20  14.
```

At the low mean 20, the `b=2` and `b=12` cells would be lifts of scaled
mass eight and six, below Proposition 15.688's sharp lift floor ten. The
phase-one `b=10` problem becomes the phase-zero `b=3` problem on its
three-point complement.  The live Proposition 15.652 quadrature has
coefficients `(1,-4,4)`, contact nodes `(1,2,3)`, and strictly positive
weights `(15/26,3/13,5/26)`.  At equality, positivity forces the original
integer quadratic—not only its stabilizer average—to equal the parity
minimum on every point in the three contact layers.

The remaining extension to the omitted `r=0` layer is certified separately.
The 78 pair monomials span all degree-at-most-two functions on `J(13,7)` by
the fixed-weight identities `sum_(i<j)x_i*x_j=21` and
`sum_(j!=i)x_i*x_j=6*x_i`.  Their `1596 by 78` evaluation matrix on the
`r=1,2,3` layers has rank 78 modulo 101, hence also over the rationals.  Its
SHA-256 is

```text
996269be45189565eaf8717f97f71f2e2f22ad33c8116da5fd2e154ec8eaf695.
```

Thus vanishing on the contact layers determines a slice quadratic globally,
and the exact `b=10` cell is pointwise, not merely an averaged floor value.
The only low cell is therefore the exact `b=10` complement-triple quadratic

```text
A=(2-r)^2,
epsilon*S_H=5-sum_(i in C) z_i+sum_({i,j} subset C) z_i*z_j.
```

Its coefficient offset is `5-3=2`. Since every quotient is at least one,
the quotient profile is forced to

```text
1^6 2.
```

Thus there are six exact complement triples of mean 20 and one elevated
hard row of mean 34.

The `b=10` cell two units above this floor needs a different argument.  The
ordinary sharp-lift theorem cannot be used: for `B=(A-(2-r)^2)/2`, parity and
nonnegativity give `B>=0` on `r=1,2,3` but only `B>=-2` on `r=0`.  Proposition
15.744 therefore imposes Proposition 15.738's 1,638 independent
third-difference identities, those exact lower bounds, and

```text
sum_X B(X)=66,          equivalently 4p*E[B]=2.
```

The safe coordinate upper bound is 306 because the omitted layer has 120
points.  The deterministic 1,716-variable, 1,639-constraint one-worker model
is `INFEASIBLE`; its text-proto SHA-256 is

```text
b0d1956f0a173f7c4ce94d7f588af92311d42f5017dd72d814d361e964b6bcd4.
```

This is the actual exclusion of the `u=4,b=10` floor-plus-two candidate.
Proposition 15.688 is used only for the ordinary globally nonnegative lifts
above the `b=2` and `b=12` baselines.

The same live floor and lift sieve applied to every residue gives

| `u` | `sum k` | forced low rows | low mean | outcome |
|---:|---:|---:|---:|---|
| 0 | 11 | at least 3 with `k=1` | 14 | exact `b=12` survives |
| 1 | 10 | at least 4 with `k=1` | 16 | excluded below lift floor |
| 2 | 9 | at least 5 with `k=1` | 18 | excluded below lift floor |
| 3 | 8 | at least 6 with `k=1` | 20 | exact `b=10` survives |
| 4 | 7 | all 7 with `k=1` | 22 | `b=10` punctured lift dies; sharp `b=2` mass 10 survives |
| 5 | 6 | at least 1 with `k=0` | 10 | below every phase-one floor |
| 6 | 5 | at least 2 with `k=0` | 12 | exact `b=2` survives |

Thus the exact pre-15.744 survivor set is `{0,3,4,6}`. Removing the branch
proved here leaves `{0,4,6}`; that scope statement is derived, not assumed.

## 2. The two parallel-count ledgers

Let `h` be the hard sign and `T` the signed global edge total. For an exact
hard row with parallel count `P`, local normalization gives

```text
20=14P-hT-39,       hT=14P-59.                    (1)
```

The offset-two congruence is `6 | P-2`. The elevated row has parallel count
`R`, and its mean equation together with (1) gives

```text
34=14R-hT-39,       R=P+1.                        (2)
```

The six exact rows and the elevated row consume `7P+1` edges. Since
`7P+1<=61`, the congruence leaves only `P=2,8`. Exact edge and mean
accounting is then

| `P` | `hT` | `R` | hard edges | `sum Q` opposite | opposite mean |
|---:|---:|---:|---:|---:|---|
| 2 | -31 | 3 | 15 | 46 | `a(Q)=14(Q-5)` |
| 8 | 53 | 9 | 57 | 4 | `a(Q)=14(Q+1)` |

In the first ledger, `Q=5,a=0` would have the constant signed target 3,
whose coefficient offset is three. It would require `6 | Q-3`, which fails
at `Q=5`. Hence every opposite row has `Q>=6`; because their sum is 46,
at least one has `Q=6,a=14`.

In the second ledger, seven nonnegative `Q` values sum to four, so at least
one (in fact at least three) has `Q=0,a=14`.

## 3. The changed `H=61` mass-14 certificate

Proposition 15.738's mass-14 Boolean support classification is independent
of the residual edge count, but its preliminary height-four exclusions are
not: those models used `|H|=59` and

```text
sum |W_st| <= 59-Q.
```

They cannot simply be cited at this layer. Proposition 15.744 rebuilds both
models with the necessary relaxed bounds

```text
Q=0: sum |W_st| <= 61,
Q=6: sum |W_st| <= 55.
```

Each model has 1,716 variables `B(X)` in `[0,4]`, 78 integral coefficients
`W_st`, exact mass 462, even coefficient-row sums, the coefficient total

```text
sum W=13Q-53,
```

and all 1,716 identities

```text
4B(X)=Q-3+sum W-2 cut_W(X).
```

The safe height-four orbit anchor is `B(first lexicographic 7-set)=4`.
Exact one-worker CP-SAT gives:

| `Q` | relaxed `l1` | status | model SHA-256 |
|---:|---:|---|---|
| 0 | 61 | `INFEASIBLE` | `70313e414ca6da2cf6694c11bdd7c7ee8ee985ca05bd30802aa2b6b96353d3d3` |
| 6 | 55 | `INFEASIBLE` | `a94796122b2c1a115b1efec4094031726f118e05fd6825f2e805a406b4f2b9dd` |

The live height dichotomy from Proposition 15.738 is `[1,4]`, so the forced
cell is Boolean. Its edge-count-independent exhaustive support-462 catalog
has three families with offsets 6, 4, and 4. At `Q=0` and `Q=6`, only the
offset-six family survives:

```text
B=x_i*x_j.
```

Both slice gauges have normalized moments `(i-j)^2` and `(i-j)^4`.

## 4. Six roots finish the branch

For even `d`, use the genuine global homogeneous binary form

```text
M_d(L)=sum_({u,v} in H) chi(u-v)(L(u)-L(v))^d.
```

In either exact hard gauge `P=2,8`, a complement triple has

```text
2S_4=S_2^2.
```

Consequently

```text
G(L)=2hM_4(L)-M_2(L)^2
```

vanishes in the six distinct exact hard projective directions. A nonzero
binary quartic has at most four projective roots, so `G` is identically
zero. The selected-pair opposite cell has sign `-h` and instead gives

```text
G=-3(i-j)^4 != 0  in F_13.
```

The executable certificate checks both choices of `h` and both slice
gauges. This contradiction excludes `p=13,t=4,u=3`.

## 5. Exact remaining gate

The other phase-one residues at `p=13,t=4` are not claimed here. The exact
remaining residue set is

```text
u in {0,4,6}.
```

Canonical artifacts for this branch are:

* `src/e1_gmin_m4_prop15744.py`
* `tests/test_prop15744.py`
* `evidence/e1_gmin_m4_prop15744.json`
* `evidence/NOTE_2026-08-31_p13_t4_u3_quartic_close.md`
