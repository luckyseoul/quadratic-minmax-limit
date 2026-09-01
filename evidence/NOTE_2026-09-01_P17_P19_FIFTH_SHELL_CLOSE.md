# Proposition 15.753: the `p=17,19` fifth-shell endpoints close

Result class: **exhaustive finite aggregate certificate and proved endpoint
theorem**.  This is a deliberately narrow finite-prime aggregate-row census.
It does not enumerate graphs, coefficient cells, or common graph
realizations.

## Statement

In residual (ii), the fifth shell `t=4`, hence `k=4p+8`, is empty at both
exceptional endpoints not covered by Proposition 15.752:

```text
p=17: k=76, |E(H)|=77;
p=19: k=84, |E(H)|=85.
```

Thus Proposition 15.752's `p>=23` fifth-shell theorem extends through
`p=17,19`.  This does **not** close residual (ii), `E(1)`, or the original
quadratic min-max limit problem.

The executable certificate is
`src/e1_gmin_m4_prop15753.py`; its exact output is
`evidence/e1_gmin_m4_prop15753.json`, and fail-when-wrong tests are in
`tests/test_prop15753.py`.

## 1. Exact branch exhaustiveness

Put

```text
q=(p-1)/2,  m=q+1,  t=4,
a_d=2u+(p+1)k_d,  sum_d k_d=m+t-u.
```

For `0<=u<=t`, first `k_d=0` would give mean `2u`, below the least phase-one
floor (`16` at `p=17`, `18` at `p=19`), so `k_d>=1`.  Compare the `k_d=1`
mean `p+1+2u` with every entry of the exact phase-one floor table.  The
largest possible lift excess is `2t+2=10`, strictly below the sharp nonzero
integral-lift floors `p-3=14,16`.  Unless equality holds, `k_d=1` is therefore
also impossible and every direction has `k_d>=2`.  Since
`sum k_d=m+t-u<2m`, such a row is empty.  At `p=17`, the only exact low-row
survivor is the `b=p-1` literal at `u=0`; at `p=19` there is no exact low-row
survivor.  For `t<u<=m-2`, `k_d=0` still has mean `2u` below the least
phase-one floor, while

```text
sum_d k_d=m+t-u<m
```

contradicts `k_d>=1` in all `m` directions.  At `u=m-1`, the exact
phase-one floor candidates are `b=2` at `p=17` and `b in {2,p-1}` at
`p=19`.  At `p=19` the two endpoint types cannot mix: their common parallel
count would have to obey both coefficient congruences

```text
P=4 (mod 9),   P=3 (mod 9).
```

Consequently the exhaustive branches are

```text
p=17: A = b=2 XNOR, or B = b=16 literal;
p=19: A = b=2 XNOR, or C = b=18 complement literal.
```

The checked ledger `hard_residue_ledger(p)` derives these lists from the
floor tables and coefficient offsets.  They are not assumptions imported
from the `t<=2` API of Proposition 15.735.

## 2. Common row identity

For a projective direction `L`, let `P_L` be its parallel edge count and
let `q_L(a)`, `1<=a<=(p-1)/2`, be the normalized nonzero part of its common
difference-Radon row.  The coefficient-offset congruence is first solved in
the full integral range for `P`; it independently gives the hard/opposite
edge split and therefore `hT=|E_h|-|E_{-h}|`.  Only then are local rows glued
to the common graph.  Their nonzero-bin sums are

```text
hard:     sum_a q_L(a)= hT-P_L,
opposite: sum_a q_L(a)=-hT-Q_L.
```

The unspecialized local hard sum is

```text
p(P_L-3)-c-(p+1)k_L,
```

where `c=0` only for the `p=17` literal branch and `c=p-1` for the three
other branches.  Equating local and common sums gives, with no alternative
integer `P_L`,

| branch | edge split `(hard,opposite)` | `hT` | exact row | forced `P_L` |
|:---|:---:|---:|:---|:---|
| `p=17` A | `(41,36)` | 5 | XNOR `e_d`, energy 1 | `4+k_L` |
| `p=17` B | `(49,28)` | 21 | literal `(2)^8`, energy 32 | `4+k_L` |
| `p=19` A | `(45,40)` | 5 | XNOR `e_d`, energy 1 | `4+k_L` |
| `p=19` C | `(35,50)` | -15 | complement literal `(-2)^9`, energy 36 | `3+k_L` |

The exact-row sums (`1,16,1,-18`) agree with `hT-P` after the independent
edge split; no exact star is allowed to choose its own normalization.  The
exact identity used in every energy comparison is

```text
sum_{L,a>0} q_L(a)^2
 = p|E(H)| + 2(hT)^2 - 2 sum_L P_L^2 + 2p C,

C = sum_delta binom(m_delta,2) >= 0.
```

If some hard directions are exact stars, their known row energies are
subtracted from the right side to obtain a lower bound for the remaining
nonexact rows.

For every row, all translation-averaged middle-slice cut inequalities are
imposed.  Their exact catalogs are:

| `p` | bins | distinct cut rows | row sum | SHA-256 |
|---:|---:|---:|---:|:---|
| 17 | 8 | 698 | 72 | `a8ac7349cb601db5163ef1526949587c766914d774fe26858fe93eac1d940708` |
| 19 | 9 | 2338 | 90 | `5f07e9ced107e6dc1551b806043a92147c00d80eb009b70d0cbfd3ce9631c5b7` |

Each finite row maximum has an explicit integral maximizer.  Its upper
bound is proved by adding `sum q_a^2 >= maximum+1` and replaying the resulting
CP-SAT model to exact `INFEASIBLE` with one worker.  An empty row system is
replayed without an energy constraint.  The model-proto hashes, solver
version, worker count, status, and maximizers are pinned in the JSON.

## 3. `p=17`, branch B (literal)

Here there are nine hard and nine opposite directions,

```text
hard edges=49, opposite edges=28, hT=21,
P_L=4+k_L,  sum hard excess=4.
```

At least five exact literal stars force `M_2=M_4=0`.  The exact row results
are

| row | result for `sum q_a^2` |
|:---|---:|
| hard excess 1 | infeasible |
| hard excess 2 | 70 |
| hard excess 4 | 218 |
| opposite `Q=3` | 72 |
| opposite `Q=4` | 101 |

Every partition containing excess one is empty.  The two remaining
Parseval ledgers are

| excess partition | hard `P` multiset | opposite `Q` multiset | lower | upper | gap |
|:---|:---|:---|---:|---:|---:|
| `(2,2)` | `5^7,7^2` | `3^8,4` | `1245+34C` | 817 | 428 |
| `(4)` | `5^8,9` | `3^8,4` | `1197+34C` | 895 | 302 |

Both are impossible.

## 4. `p=17`, branch A (XNOR)

Here

```text
hard edges=41, opposite edges=36, hT=5,
P_L=4+k_L,  sum hard excess=5.
```

An opposite `Q=3` cell has scaled mean `8`.  If `b` is nonzero this is below
the least nonzero phase-zero floor `16`; if `b=0` it is a nonzero integral
lift below the sharp floor `14`.  Hence all nine opposite directions have
`Q=4`.
The raw translated-cut maxima are

```text
hard excess e=1,2,3,4: 28,81,200,289;
opposite Q=4:            23.
```

They give the following exact ledgers for every partition except `(5)`:

| partition | lower | upper | gap |
|:---|---:|---:|---:|
| `(1,1,1,1,1)` | `689+34C` | 347 | 342 |
| `(2,1,1,1)` | `684+34C` | 372 | 312 |
| `(2,2,1)` | `679+34C` | 397 | 282 |
| `(3,1,1)` | `675+34C` | 463 | 212 |
| `(3,2)` | `670+34C` | 488 | 182 |
| `(4,1)` | `662+34C` | 524 | 138 |

For partition `(5)`, eight of the nine hard directions are exact XNOR rows.
Let

```text
S_j(L)=epsilon_L M_j(L),
epsilon_L=h on hard directions,
epsilon_L=-h on opposite directions.
```

The homogeneous binary quartic

```text
G(L)=h M_4(L)-M_2(L)^2
```

vanishes at the eight distinct exact hard directions, more than its degree,
so `G` is identically zero over `F_17`.  The normalization sign is essential:

```text
hard:     S_4=S_2^2  (mod 17),
opposite: S_4=-S_2^2 (mod 17).
```

With those relations the exact maxima are

```text
hard excess 5: 384,
opposite Q=4:   11.
```

The opposite maximizer is

```text
(-1,-1,-2,-1,-1,-1,-1,-1),
```

with `(S_2,S_4)=(8,4)` modulo `17`, so `4=-8^2 (mod 17)`.  The formerly
tempting plus-sign row has energy `15` and `(S_2,S_4)=(15,4)`; it satisfies
`S_4=S_2^2` but violates the required opposite relation.  A dedicated
regression test rejects it.

The final ledger is therefore

```text
lower = 645+34C,
upper = 384+9*11 = 483,
gap   = 162.
```

Thus branch A is empty, and `p=17,k=76` is closed.

## 5. `p=19`, branch C (complement literal)

Here

```text
hard edges=35, opposite edges=50, hT=-15,
P_L=3+k_L,  sum hard excess=5.
```

At least five exact complement literals force `M_2=M_4=0`.  An opposite
`Q=4` cell again has scaled mean `8`.  This is below the least nonzero-`b`
phase-zero floor `20`, while the `b=0` option is a nonzero integral lift
below the sharp floor `16`.  Thus all ten opposite directions have `Q=5`.
Its normalized row
would have to satisfy

```text
sum q_a=10,  sum |q_a|<=80,
c_X dot q<=114 for all 2338 translated cuts,
M_2=M_4=0 (mod 19).
```

The exact one-worker model is infeasible.  Branch C is empty.

## 6. `p=19`, branch A (XNOR)

Here

```text
hard edges=45, opposite edges=40, hT=5,
P_L=4+k_L,  sum hard excess=5.
```

An opposite `Q=3` cell has scaled mean `8`; the nonzero-`b` floor is `20`
and the `b=0` integral-lift floor is `16`.  Hence all ten opposite directions
have `Q=4`.  The raw maxima are

```text
hard excess e=1,2,3,4,5: 36,97,194,325,494;
opposite Q=4:             23.
```

The exact Parseval ledgers are

| partition | lower | upper | gap |
|:---|---:|---:|---:|
| `(1,1,1,1,1)` | `930+38C` | 410 | 520 |
| `(2,1,1,1)` | `925+38C` | 435 | 490 |
| `(2,2,1)` | `920+38C` | 460 | 460 |
| `(3,1,1)` | `916+38C` | 496 | 420 |
| `(3,2)` | `911+38C` | 521 | 390 |
| `(4,1)` | `903+38C` | 591 | 312 |
| `(5)` | `886+38C` | 724 | 162 |

All are impossible.  Branch A is empty, and `p=19,k=84` is closed.

## Verification and scope

Run:

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15753.py
PYTHONPATH=src pytest -q tests/test_prop15753.py
```

The first command regenerates the JSON atomically after all nineteen exact
one-worker row-model replays.  The evidence manifest pins the two cut
catalogs, all row specifications, every CP model, and the source/test/note
files.  The tests independently replay every explicit maximizer, all branch
and Parseval arithmetic, the sign-sensitive opposite quartic model, both
empty hinge models, the live artifact hashes, atomic output, and the actual
global residual predicate.

The proposition sets only

```text
p17_k76_closed=True,
p19_k84_closed=True.
```

It keeps

```text
residual_ii_k_ge_4p_ND_closed=False,
E1_closed=False,
quadratic_minmax_limit_closed=False.
```
