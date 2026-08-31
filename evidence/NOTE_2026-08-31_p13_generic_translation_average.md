# Generic `p=13,t=3`: translation averages leave only four exact stars

**Status.** Proposition 15.740 is a proved branch split with an exhaustive
finite certificate. It excludes the two generic `p=13,t=3` hard-quotient
partitions having five or six exact hard stars. It leaves the partition

```text
(1,1,1,1,2,2,2)
```

open. Thus it does not close `p=13,k=58`, residual (ii), multi-level Type I,
or the quadratic min-max limit.

## 1. Duplicate gate and exact branch split

The binary affine-Radon isomorphism is already Proposition 15.692. Over
`F_2`, if `A` is the affine line-point incidence matrix, then

```text
A^T A=I+J.
```

On even point words this is an isomorphism onto the direct sum of the even
directional profile blocks, with inverse `x=A^T r`. In the present branch all
seven opposite profiles are zero, so the boundary is reconstructed from the
seven hard even profiles, an 84-dimensional direct summand. An exact hard
profile is `1+delta_j`. This application imports 15.692; it is not a new
Radon theorem.

In the generic p=13 branch, write the seven hard means as `a_L=14k_L`. The
ledger is

```text
k_L>=1,        sum_L k_L=10.
```

The three excess units have exactly three unordered distributions:

| hard quotient partition | exact `k=1` stars |
|---|---:|
| `(1,1,1,1,1,1,4)` | 6 |
| `(1,1,1,1,1,2,3)` | 5 |
| `(1,1,1,1,2,2,2)` | 4 |

For even `d<12`, every exact star is a zero of the global homogeneous form

```text
M_d(L)=sum_({u,v} in H) chi(u-v)*(L(u)-L(v))^d.
```

The power sum over a star is zero in `F_13`. More than `d` projective roots
force a binary degree-`d` form to vanish identically. Hence five or six
exact stars force both `M_2=M_4=0`; four exact stars force `M_2=0` but do not
force `M_4=0`.

## 2. The opposite signed-matrix relaxation

Every opposite direction has `Q=3`, scaled mean 20, and `b=0`. Put
`W=epsilon_L K^L`. Exact coefficient comparison gives

```text
sum_(i<j) W_ij=-20,
sum_(i<j) |W_ij|<=56,
every row sum is even,
B(X)=-5-cut_W(X)/2>=0                  (|X|=7).
```

Thus every balanced cut has weight at most `-10`. Conditioning a uniform
seven-set on containing a fixed pair gives

```text
E[B | i,j in X]=(20+12W_ij)/44.
```

Nonnegativity and integrality imply `W_ij>=-1`.

For `a=1,...,6`, let `E_a` be the thirteen unordered pairs of cyclic
difference `+-a` and put

```text
n_a=sum_({i,j} in E_a) W_ij.
```

The six aggregates satisfy

```text
-13<=n_a<=18,
sum_a n_a=-20,
sum_a |n_a|<=56,
sum_a a^2 n_a=0 mod 13,
sum_a a^4 n_a=0 mod 13.
```

The lower bound uses the thirteen entries in a distance class and
`W_ij>=-1`. The upper bound uses the total positive mass: from sum `-20` and
`l1<=56`, it is at most `(56-20)/2=18`.

This is a necessary aggregate relaxation. A feasible aggregate vector would
not by itself be a directional cell, a common 59-edge graph, or a residual
witness. Infeasibility, however, excludes every genuine opposite cell.

## 3. Translation-average inequalities

For a seven-set `X subset F_13`, define

```text
c_a(X)=|X triangle (X+a)|.
```

For a fixed edge of difference class `+-a`, exactly `c_a(X)` of the thirteen
translates `X+t` separate its endpoints. Summing their thirteen cut
inequalities therefore gives the exact necessary inequality

```text
sum_(a=1)^6 c_a(X)n_a<=-130.                 (1)
```

The 1,716 seven-sets give 74 distinct vectors `c(X)`. Every entry is even,
each vector sums to 42, and their catalog hash is

```text
bfec2077a81acf1a6719caf93b066313445c55b4e2951c189d357731b437a265.
```

## 4. Exact nine-vector certificate

Row reduction of the sum, quadratic-moment, and quartic-moment congruences
over `F_13` gives, with `n_4,n_5,n_6` free,

```text
n_1+10n_4+n_5+10n_6 = 9 mod 13,
n_2+ 6n_4+3n_5+ 2n_6 =12 mod 13,
n_3+11n_4+10n_5+2n_6 =11 mod 13.
```

Enumerating `n_a in [-13,18]`, then imposing the exact sum and l1 bound,
leaves 32,313 vectors. An independent meet-in-the-middle enumeration gives
the same count. Their deterministic catalog hash is

```text
a06fcc7e9b9d793babb49f480f261d02683be0bc4df5bcc706fe785a23076b6d.
```

From the lexicographically sorted 74-vector catalog, repeatedly choose the
first vector eliminating the most remaining aggregates. The resulting nine
vectors and lexicographically first representative seven-sets are:

| `c(X)` | representative `X` |
|---|---|
| `(2,4,6,8,10,12)` | `(0,1,2,3,4,5,6)` |
| `(8,8,6,8,8,4)` | `(0,1,2,4,7,8,10)` |
| `(8,8,8,6,4,8)` | `(0,1,2,5,6,8,10)` |
| `(12,2,10,4,8,6)` | `(0,1,3,5,7,9,11)` |
| `(4,8,10,8,6,6)` | `(0,1,2,3,4,8,9)` |
| `(8,8,4,8,6,8)` | `(0,1,2,4,5,7,10)` |
| `(6,12,8,2,4,10)` | `(0,1,2,5,6,9,10)` |
| `(10,6,4,12,2,8)` | `(0,1,3,5,6,8,11)` |
| `(8,8,6,6,10,4)` | `(0,1,2,4,7,8,11)` |

They eliminate respectively

```text
14222, 9967, 6087, 1395, 417, 168, 45, 8, 4
```

rows, leaving

```text
18091, 8124, 2037, 642, 225, 57, 12, 4, 0.
```

Thus the nine instances of (1) are already infeasible. Their hash is

```text
2932ffd7f49846c390df8076991141cc496226055dc6c0c1daa18399b28116bd.
```

As an independent exact check, a 14-variable, 19-constraint CP-SAT model
containing only the six `n_a`, their absolute values, two moment quotients,
and the nine inequalities returns `INFEASIBLE` with one worker and seed zero.
Its text-proto hash is

```text
3aca686a9aff08c5118c2a95280d0bd17d17a410365fc5117f428be82e879826.
```

The pure integer enumeration is the primary exhaustive certificate; the
solver is a separately encoded validation.

## 5. Exact conclusion and next gate

If the hard quotient partition has five or six exact stars, every opposite
cell would enter the infeasible aggregate model above. Therefore both

```text
(1,1,1,1,1,1,4),        (1,1,1,1,1,2,3)
```

are impossible. The generic p=13 branch is reduced exactly to

```text
(1,1,1,1,2,2,2).
```

The remaining theorem must couple four exact `P=5` stars, three elevated
`P=6` hard cells, and seven `Q=3,b=0` opposite cells through the common
59-edge graph. The quartic has exactly four forced hard roots and can be a
nonzero scalar multiple of their product, so the present `M_4=0` certificate
does not apply to that last partition.
