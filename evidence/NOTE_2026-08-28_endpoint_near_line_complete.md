# Complete endpoint-only near-line exclusion

Date: 2026-08-28. This is Proposition 15.673. It excludes both product
signs of every infinity-plus-`(p-2)` boundary for which every directional
odd-fibre count is an endpoint,

```text
b_d in {1,p-2},
```

for every odd prime `p>=17`. Propositions 15.671--15.672 supply the
collinear case. The result does not exclude non-endpoint profiles of the
same boundary size and therefore does not close residual (ii).

## 1. Exact means and endpoint baselines

Put

```text
q=(p-1)/2,    m=(p+1)/2,    P=p+1.
```

For a direction `d`, let `I` be the number of infinity edges, `P_d` the
number of finite selected edges parallel to `d`, and `T` the total
Paley-signed sum of the finite selected edges. Proposition 15.632 gives the
exact directional mean

```text
a_d = I + P P_d - eps_d T - 3p.                         (1)
```

Hence all `m` means of one quadratic direction type have one even residue
`r=2u` modulo `P`. Write

```text
a_d=r+P k_d,    sum_d k_d=m-u.                          (2)
```

At the two endpoint fibre counts, the parity baselines and scaled costs are

| `p mod 4` | phase | `b=1` cost | `b=p-2` cost |
|---:|---:|---:|---:|
| 3 | 0 | `P` | `P` |
| 1 | 1 | `P-2` | `P` |
| 3 | 1 | `P-2` | `P-2` |
| 1 | 0 | `P` | `P-2` |

A nonzero lift above one of these baselines costs at least four scaled
units by Proposition 15.642. Thus a mean only two units above its baseline
is impossible. Combining that fact with (2) leaves only baseline directions
and, in the nonsaturated cases, one `P`-unit mean jump in a type.

## 2. The endpoint geometry needed by the reduction

Let `R` be the number of directions with `b_d=1` and let `s=p-2`. Every
pair of finite boundary points collides in exactly one direction, so

```text
R(s-1) <= sum_d(s-b_d) <= s(s-1),
```

and therefore `R<=s`. Equality deserves care. It makes every fibre contain
at most two boundary points, so the boundary is an `s`-arc, and exactly
three directions are undetermined.

Choose any two of those three points at infinity. Adjoining them to the
`(p-2)`-arc gives a `p`-arc in `PG(2,p)`. Segre's odd-order `p`-arc theorem
puts it on a conic. Repeating this with another pair gives a second conic
sharing the `p-2>=5` finite points, so the conics coincide. The common conic
would then contain all three collinear infinity points, impossible for a
nondegenerate conic. Hence the equality case is impossible.

If at most two directions are determined, the finite set is collinear:
three noncollinear points already determine three distinct directions.
That branch is closed by Propositions 15.671--15.672.

After those two exits,

```text
3 <= R <= p-3 = 2m-4.
```

The same-type residue equation now gives the following exhaustive type
forms. Here `B` means that the baseline directions have `b=1`, `C` means
that they have `b=p-2`, and the listed alternatives allow the unique
exception to have either endpoint kind.

| branch | `B`-type possible `b=1` counts | `C`-type possible `b=1` counts | edge offsets |
|---|---:|---:|---:|
| `p=3 mod 4`, phase 0 | `m` | `0` | `0,0` |
| `p=1 mod 4`, phase 1 | `m-1,m` | `0` | `1,0` |
| `p=3 mod 4`, phase 1 | `m-1,m` | `0,1` | `1,1` |
| `p=1 mod 4`, phase 0 | `m` | `0,1` | `0,1` |

For completeness, in the equal-floor rows a type cannot have baseline
directions of both endpoint kinds. Equal baseline means force equal `P_d`
by (1), whereas the two coefficient congruences in the next section differ
by one modulo `q`. In the unequal-floor rows the table follows directly
from (2): the only possibilities are the saturated endpoint baseline, or
`m-1` low baselines and one direction one period higher.

Pairing two `B` types puts `R>2m-4`; pairing two `C` types puts `R<3`.
Thus this bound forces exactly one `B` type and one `C` type in every branch.
Their total edge offsets are respectively `0,1,2,1`, which is the claimed
two-count normal form.

## 3. Coefficient congruences

On the middle slice write `z_s=2x_s-1`. A baseline with one odd fibre has

```text
eps_d S_H = 4 + sigma z_j,    sigma=(-1)^phase,
```

while a complementary endpoint has

```text
eps_d S_H = 4 + tau z_a z_b,
```

with `tau=-1` in the two rigid-sign rows of the table and `tau=+1` in the
two opposite-sign rows. Comparing cross coefficients gives

```text
b=1:    q | I+P_d-(4+sigma),
b=p-2:  q | I+P_d-4.                                  (3)
```

Let `x` be the common baseline parallel count in the `b=1` type and `y`
the corresponding count in the `b=p-2` type. The one-exception normal forms
and (3) give the complete reduction:

| branch | finite edges `E` | substituted congruences |
|---|---:|---|
| `p=3 mod 4`, phase 0 | `m(x+y)` | `q|y`, `q|(x-1)` |
| `p=1 mod 4`, phase 1 | `m(x+y)+1` | `q|x`, `q|(y-1)` |
| `p=3 mod 4`, phase 1 | `m(x+y)+2` | `q|y`, `q|(x+1)` |
| `p=1 mod 4`, phase 0 | `m(x+y)+1` | `q|x`, `q|(y+1)` |

In all rows

```text
I=4p+1-E >= 1
```

forces `x+y<=7`.

For `p=3 mod 4`, one has `q>=9`. Phase zero therefore has the unique
candidate `(x,y)=(1,0)`. It has

```text
E=m,    I=7q+4,    s+2E=4q+1,
```

contradicting the elementary boundary support inequality `I<=s+2E`.
Phase one has no candidate because `q|(x+1)` but `1<=x+1<=8<q`.

For `p=1 mod 4`, phase one has the unique candidate `(0,1)`. Here

```text
E=m+1,    I=7q+3,    s+2E=4q+3,
```

again impossible. Phase zero has no candidate for `q>8`. Only `q=8`, or
`p=17`, remains.

The boundary support inequality used above is immediate. The `I` infinity
edges toggle `I` finite vertices. The boundary of the `E` finite edges can
toggle at most `2E` vertices, and their symmetric difference is the
`s`-point finite boundary. Thus `I<=s+2E`.

## 4. The exact `p=17` endpoint

At `p=17`, phase zero, the sole arithmetic candidate is

```text
(x,y)=(0,7),    E=64,    I=5.
```

Take one of the baseline complementary directions. If `n_s` counts its
infinity-star endpoints in fibre `s`, then `sum n_s=5`, `P_d=7`, and the
coefficient scalar satisfies `2c=1`. Its signed inter-fibre matrix is

```text
L_st = 1-n_s-n_t+1_{st=ab}.                            (4)
```

Suppose exactly `u` of the 17 counts are positive and put `z=17-u`.
Before the distinguished `ab` correction, separating zero-zero,
zero-positive, and positive-positive pairs gives

```text
C(z,2) + z(5-u) + 5(u-1) - C(u,2).
```

The distinguished correction can lower this by at most one. For
`u=1,...,5`, the resulting lower bounds are

```text
183, 153, 125, 99, 75.
```

The last value is attained by five unit counts including fibres `a,b`, so
the exact minimum is 75. But only

```text
E-P_d=64-7=57
```

transverse selected edges exist, whereas their signed sums must have
entrywise `l1` norm at least 75. This final contradiction closes `p=17`.

## 5. Literature and OEIS context

The geometric equality step uses the classical odd-order `p`-arc
classification, restated in Ball--Lavrauw,
[Planar arcs](https://arxiv.org/abs/1705.10940): Segre classified arcs of
sizes `q` and `q+1`, and an odd-order `q`-arc lies on a conic. The use here
is narrower: two extensions of one endpoint boundary force three collinear
points onto that conic. Targeted searches found no source combining this
arc extension with the Paley type budgets, mean quantization, and additive
inter-fibre coefficient obstruction above.

OEIS searches for the sample infinity-count values and for the final pair
`75,57` returned only unrelated occurrences. The quantities here are direct
linear or finite `l1` formulas, not a proposed integer sequence, and no OEIS
submission is claimed.

## 6. Reproduction and scope

```bash
python src/e1_gmin_m4_prop15673.py
python -m pytest -q tests/test_prop15673.py
```

The generated exact record is
`evidence/e1_gmin_m4_prop15673.json`. The independent test enumerates all
weak compositions of five into 17 fibres and reproduces the exact minimum
75.

The endpoint condition is essential. Directions with
`3<=b_d<=p-4`, large all-finite boundaries, general residual (ii), R1,
global QVAR, Type I, and the limit all remain open.
