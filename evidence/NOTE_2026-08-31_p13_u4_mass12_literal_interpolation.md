# p=13,t=4,u=4: mass-12 exclusion and literal-root interpolation

**Propositions:** 15.747--15.748
**Status:** proved `P=3` branch exclusion plus exhaustive finite open
reduction of `P=5`

## Mass-12 cut obstruction

For a phase-zero `b=0` mass-12 lift `A=2C`, normalized coefficients satisfy

```text
S=sum W=13Q-51,       cut_W=7Q-27-2C.
```

If `C` is Boolean, it has support 396. For a uniform six-set, one edge is
cut with probability `7/13`, two adjacent edges with `7/26`, and two
disjoint edges with `42/143`. With `E2=sum W_e^2` and squared row-degree sum
`D2`, this gives

```text
E[cut_W^2]=(-7D2+84E2+84S^2)/286.
```

Equating the required two-valued cut moment yields

```text
-7D2+84E2+182Q^2-1428Q+2598=0.
```

The residual is one modulo seven for every integer `Q`; height one is
impossible without the `l1` bound, row parity, third differences, or field
moments.

For height four, dropping field moments restores `S_13` symmetry and permits
one maximum point to be anchored. The exact projected models have 169
variables and 3,526 constraints:

| Q | sum W | l1 | cut interval | model SHA-256 | status |
|---:|---:|---:|---:|---|---|
| 3 | -12 | 58 | `[-14,-6]` | `e8404a5684e033b73750b1f36a338aa13038861d6dbfc614cc99b6f0666423d9` | `INFEASIBLE` |
| 5 | 14 | 56 | `[0,8]` | `8f992368fac869f29c23e6ecd20400228c2c10d5bda4d1001b291242dd6e3941` | `INFEASIBLE` |

Both are deterministic one-worker, no-timeout certificates. Therefore the
omitted-pair `P=3` branch is empty, while every minimum `P=5,Q=3` cell is the
exact `b=12` literal.

## Literal-root interpolation

Let `z` be the number of those minimum literals. They are common roots of
`M2,M4,M6`, and `z>=2` because the seven opposite excesses sum to five.
The hard baseline-pair/all-equal-triple moment alphabet contains exactly 69
triples. When `M2=0`, its nonzero `M4` alphabet is `{7,8,11}`.

- `z>=5` is impossible by the quartic root bound.
- Exact `M4` interpolation gives zero `z=4` and `z=3` cases for either hard
  sign.
- For `z=2`, write `M2=cR2`, `M4=R2Q2`, `M6=R2Q4`. Per sign, 1,554 first-two
  moment candidates and 2,688 `N6` vectors leave exactly 336 distinct
  moment-level survivors.

The raw interpolation payload hash is
`894c087d4acae7ff0722ba236b1fac494984b9b331431e6117b2edbde0afbbec`.
The survivors are not common graphs. They force only

```text
z=2,       positive opposite excess partition=(1,1,1,1,1).
```

The next gate is to couple those 336 survivors per hard sign to the five
excess-one `Q=4` cells and a single common 61-edge graph.

## Artifacts

- `src/e1_gmin_m4_prop15747.py`
- `src/e1_gmin_m4_prop15748.py`
- `scripts/p13_p5_literal_interpolation.py`
- `tests/test_prop15747.py`
- `tests/test_prop15748.py`
- `evidence/e1_gmin_m4_prop15747.json`
- `evidence/e1_gmin_m4_prop15748.json`
