# The full translated-cut energy certificate at `p=17,t=3`

**Status.** Proposition 15.743 is an **exhaustive finite certificate**. It
closes the generic residual-(ii) branch at

```text
p=17,  t=3,  k=4p+6=74.
```

It does not close later `p=17` layers, the same layer at every prime, global
residual (ii), multi-level Type I, or the quadratic min-max limit.

## 1. The two-source normalization is the global hinge

Propositions 15.735 and 15.739 leave nine hard directions with

```text
k_L >= 1,       sum_L k_L=12,
```

and nine opposite directions, each having parallel count `Q=3`. At least
six hard directions are exact literal cells. For an exact cell `k_L=1`,
local coefficient comparison and the common difference-Radon sum give

```text
sum q_L=17(P_L-3)-18=17P_L-69=hT-P_L,
hence hT=18P_L-69.
```

Because `hT` belongs to the one common graph, all exact cells therefore
have one common parallel count `P`. Six such directions among 75 edges give
`6P<=75`, hence `P<=12`. Exact literal coefficient comparison independently
gives `P=5 mod 8`, so `P=5`. Only after fixing this normalization do we get

```text
hT=18*5-69=21,       q_L=(2)^8.                   (1)
```

This agrees independently with the opposite ledger:

```text
opposite edges = 9*Q=9*3=27,
hard edges     = 75-27=48,
hT             = 48-27=21.
```

For one hard cell with parallel count `P_L` and quotient `k_L`, local
coefficient comparison gives

```text
sum_(a=1)^8 q_L(a)=17(P_L-3)-18k_L.               (2)
```

But the same row is the nonzero-bin part of the common difference-Radon
transform. Its hard-sign-normalized total is independently

```text
sum_(a=1)^8 q_L(a)=hT-P_L=21-P_L.                 (3)
```

Equating (2) and (3) forces

```text
P_L=4+k_L.                                        (4)
```

This is the genuinely cross-direction input. Locally feasible coefficient
cells with `P_L != 4+k_L` fail (3) and therefore cannot be projections of one
common graph.

Write `e_L=k_L-1`. The hard excess has total three, so its only unordered
partitions are

```text
1+1+1,       2+1,       3.                        (5)
```

For a nonexact hard row of excess `e`, (4) gives

```text
P_L=5+e,
sum q_L=16-e,
sum |q_L| <= 70-e.
```

The cell identity gives `cut_W(X)<=9`, hence every translation-summed cut
has upper bound `17*9=153`. For an opposite row,

```text
sum q_L=-(hT+Q)=-24,
sum |q_L|<=72,
cut_W(X)<=-12,
```

so its translated upper bound is `-204`.

## 2. Eight distance bins and all translated cuts

For `a=1,...,8`, let `q_L(a)` be the sum of the normalized signed
coefficients on the 17 unordered pairs of cyclic distance `+-a`. Six exact
hard stars force the global forms `M_2=M_4=0`, so every row satisfies

```text
sum a^2 q_L(a)=0 mod 17,
sum a^4 q_L(a)=0 mod 17.                           (6)
```

For a nine-set `X`, put

```text
c_X(a)=|X triangle (X+a)|.
```

The 24,310 nine-sets yield 698 distinct vectors `c_X`; every vector has even
entries in `[0,16]`, coordinate sum 72, and catalog digest

```text
a8ac7349cb601db5163ef1526949587c766914d774fe26858fe93eac1d940708.
```

The exact translation identity is

```text
sum_(t in F_17) cut_W(X+t)=c_X . q_L.              (7)
```

All 698 inequalities from (7) are imposed. This is not another full
24,310-point coefficient model: each row has only eight integer variables.

## 3. Exact row certificates

The primary one-worker CP-SAT model uses only:

```text
integer q in [-l1,l1]^8,
the exact sum and l1 bound,
the two modular moments (6),
all 698 translated upper cuts (7).
```

It deliberately uses no prior energy cap, entry alphabet, aggregate
coordinate bounds, lower cut, row parity, 136-entry coefficient matrix, or
complete-domain Boolean values. An independent model encodes each triple
`(q,|q|,q^2)` by an allowed-assignment table and excludes the next energy.

| row | primary result | sharp energy | primary model SHA-256 | table replay SHA-256 |
|---|---|---:|---|---|
| hard `e=1` | `INFEASIBLE` | -- | `a049d78a70996a1fbb39fc54e82972b865faadeb8a2da3b3f77c6241e4d6c229` | `9b042e511bb05d67022a80784366566d36302e3b480c56940e4b7b99d47e6188` |
| hard `e=2` | `OPTIMAL` | 70 | `4114629c63f660eeb6c5d59451c2fed7aec658a42ad7ec8b0f562b53800cce0a` | `75845543e8cf52119c926afeb464b2cdaf42f26f5ad311f583645f90d7dc8460` |
| hard `e=3` | `OPTIMAL` | 119 | `b1825c0634dda8d5d1be0a3105e668bc6217594df9a2c39e74e061958a725a56` | `50432c39725cb4e131b3be37d4b86e3e2d95156a9ad7c5b6cb5e068ce105caf7` |
| opposite | `OPTIMAL` | 72 | `3f01c71e80cedf0f5d108f1db9c0088ccd6f728537b9c75a61b8e0530c51f92c` | `44515f63b8e2944601cfad513895563ed672a9f2f40b78e2fd2327ab9ea32724` |

Explicit maximizers are

```text
hard e=2:  (1,-2,5,3,-1,2,5,1),
hard e=3:  (6,4,-1,-3,-2,-1,4,6),
opposite:  (-3,-3,-3,-3,-3,-3,-3,-3).
```

For the opposite row, fixed sum `-24` gives the exact identity

```text
8 sum q_a^2 - (sum q_a)^2 = sum_(a<b)(q_a-q_b)^2 >=0.
```

Thus its energy is at least 72, with equality only when every coordinate is
`-3`. The solver upper bound 72 and the displayed row therefore prove
uniqueness, not just existence of one maximizer.

## 4. The common Parseval contradiction

The p17 difference-Radon incidence matrix has 144 displacement columns and
162 rows including zero bins. Its exact Gram matrix is

```text
B^T B=17I+2J-G_parallel.
```

If `C=sum_delta binom(m_delta,2)>=0`, removing the exact-star rows gives the
following nonexact energies for (5):

| hard excess partition | exact hard stars | nonexact Parseval energy |
|---|---:|---:|
| `1+1+1` | 6 | `1287+34C` |
| `2+1` | 7 | `1251+34C` |
| `3` | 8 | `1211+34C` |

The first two partitions contain an infeasible hard `e=1` row. In the last
partition, the aggregate upper bound is

```text
119 + 9*72 = 767,
```

whereas the common graph requires at least 1211. The gap is 444. Hence all
three partitions are impossible and `p=17,k=74` is closed.

## 5. Exact scope and discarded route

No full-solution row count or spectral row cap is part of this certificate.
An exploratory eight-dilate spectral bound omitted an eigenvalue and was
discarded before this proposition was recorded. The broad-domain threshold
models above do not consume it.

The remaining residual-(ii) scope includes critical `p=5,7`, `p=11` at
`k>=50`, `p=13` at `k>=60`, `p=17` at `k>=76`, all `p>=17` layers `t>=4`,
and the generic branch-B `t=3` primes `p=1 mod 4` beginning at `p=29`.

Canonical artifacts:

* `src/e1_gmin_m4_prop15743.py`
* `tests/test_prop15743.py`
* `evidence/e1_gmin_m4_prop15743.json`
