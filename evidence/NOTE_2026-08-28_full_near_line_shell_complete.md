# Complete infinity-plus-`(p-2)` shell exclusion

Date: 2026-08-28. This is Proposition 15.674. It removes the endpoint
hypothesis from Proposition 15.673 and proves that an odd-degree boundary
consisting of infinity and `p-2` finite points is impossible for both product
signs and every odd prime `p>=17`, regardless of its directional odd-fibre
profile. Larger infinity-present boundaries remain open.

## 1. All odd-fibre floors fit one narrow band

Put

```text
P=p+1=2m,    q=(p-1)/2,    s=p-2.
```

For every direction, `b_d` is odd and lies between `1` and `s`. Proposition
15.669 gives the complete floor table. The endpoint floors are

| `p mod 4` | phase | `b=1` | `b=s` |
|---:|---:|---:|---:|
| 3 | 0 | `P` | `P` |
| 1 | 1 | `P-2` | `P` |
| 3 | 1 | `P-2` | `P-2` |
| 1 | 0 | `P` | `P-2` |

Every intermediate odd count `3<=b<=p-4` has floor either `2P-8` or
`2P-2`. In particular, every floor is at least `P-2`, every intermediate
floor is strictly above `P`, and every intermediate floor is at most
`2P-2`.

## 2. Exact type sums permit only two residues

For the `m` directions in one quadratic type, the exact mean formula gives

```text
a_d = 2u + P k_d,    0<=u<m,    sum_d k_d=m-u.          (1)
```

The nonzero-lift theorem also forbids a value exactly two above its parity
floor.

If `1<=u<=m-2`, every directional floor exceeds the residue `2u`, so every
`k_d>=1`. That would give `sum k_d>=m>m-u`, impossible.

If `u=0`, every `k_d>=1`. A floor-`P-2` direction cannot take value `P`,
because that is the forbidden two-unit lift, and an intermediate direction
cannot take value `P` because its floor is larger. Since `sum k_d=m`, all
directions must therefore have floor and mean `P`.

If `u=m-1`, then `sum k_d=1`. Exactly `m-1` directions have `k_d=0`, so
they have floor and mean `P-2`. The remaining direction has mean `2P-2`.
It may have any odd-fibre count: an intermediate floor is either met exactly
or raised by six, never by the forbidden two.

Thus an intermediate odd-fibre count can occur only in the unique high
direction of a type. There is at most one such direction per type.

## 3. Geometry forces opposite endpoint baseline types

When both endpoint kinds share one floor, they still cannot both occur among
the baseline directions of one quadratic type. Equal baseline means force
equal parallel-edge counts `P_d`. But the baseline coefficient identities
give

```text
b=1:    q | I+P_d-(4+sigma),
b=s:    q | I+P_d-4,
```

whose right sides differ by one modulo `q`.

Call a type a `B` type if its baseline directions have `b=1`, and a `C`
type if they have `b=s`. The four branches now have the exhaustive forms

| branch | `B` baseline count / exceptions | `C` baseline count / exceptions | total offset |
|---|---:|---:|---:|
| `p=3 mod 4`, phase 0 | `m / 0` | `m / 0` | `0` |
| `p=1 mod 4`, phase 1 | `(m-1) / 1` | `m / 0` | `1` |
| `p=3 mod 4`, phase 1 | `(m-1) / 1` | `(m-1) / 1` | `2` |
| `p=1 mod 4`, phase 0 | `m / 0` | `(m-1) / 1` | `1` |

Two `B` types contain at least `2(m-1)=p-1>s` directions with `b=1`.
Each contributes `s-1` to the pair deficit, contradicting

```text
sum_d(s-b_d) <= s(s-1).
```

Two `C` types have at most two exceptional directions. Since `b=s` means
that the direction contains no pair of boundary points on a common fibre,
the finite set would determine at most two directions. Three noncollinear
points already determine three directions, so the set is collinear. That
case was closed by Propositions 15.671--15.672.

Consequently every branch has one `B` type and one `C` type. The total
finite-edge offsets are exactly `0,1,2,1`, just as in Proposition 15.673.

## 4. The endpoint arithmetic closes the full shell

Let `x` and `y` be the parallel-edge counts of the `B` and `C` baselines.
The same coefficient substitution now applies without any endpoint
hypothesis on the exceptional directions:

| branch | finite edges `E` | congruences |
|---|---:|---|
| `p=3 mod 4`, phase 0 | `m(x+y)` | `q|y`, `q|(x-1)` |
| `p=1 mod 4`, phase 1 | `m(x+y)+1` | `q|x`, `q|(y-1)` |
| `p=3 mod 4`, phase 1 | `m(x+y)+2` | `q|y`, `q|(x+1)` |
| `p=1 mod 4`, phase 0 | `m(x+y)+1` | `q|x`, `q|(y+1)` |

Since `I=4p+1-E>=1`, one has `x+y<=7`. For `p=3 mod 4`, phase zero has
only `(x,y)=(1,0)`, which violates `I<=s+2E`, and phase one has no
candidate. For `p=1 mod 4`, phase one has only `(0,1)`, again violating
the support inequality, while phase zero has no candidate when `q>8`.

At `p=17,q=8`, phase zero again leaves only

```text
(x,y)=(0,7),    E=64,    I=5.
```

The `C` type has at least eight complementary baseline directions. Choosing
one gives exactly the same inter-fibre matrix as in Proposition 15.673,

```text
L_st=1-n_s-n_t+1_{st=ab},    sum n_s=5.
```

Its exact entrywise `l1` minimum is 75, while only `E-P_d=57` transverse
edges are available. This contradiction does not depend on the odd-fibre
count of the exceptional direction.

Therefore

```text
every infinity-plus-(p-2) boundary is impossible for both signs
and every odd prime p>=17.
```

## 5. Independent checks, literature, and OEIS

`tests/test_prop15674.py` independently enumerates all common residues,
floor classes, and relaxed lifts for the first prime in each modulo-four
class. It reproduces the saturated/one-exception dichotomy without calling
the symbolic classifier.

The adjacent finite-geometry literature studies directions determined by
point sets, including Lev, [Point distribution and perfect directions in
`F_p^2`](https://arxiv.org/abs/1903.01518), and the Rédei direction
theorems cited there. Those results are compatible with the elementary
determined-direction exit above, but the searched sources do not combine it
with Paley type sums, parity floors, or the coefficient congruences.

Targeted OEIS searches for the arithmetic-candidate values
`67,81,109,151,165,207`, `59,101,129,143,185,213`, and the final pair
`75,57` returned unrelated prime-filtered, partition, Riordan-array, and
Pythagorean entries. These values are evaluations of linear formulas and a
single finite norm calculation, not a proposed sequence.

## 6. Reproduction and remaining scope

```bash
python src/e1_gmin_m4_prop15674.py
python -m pytest -q tests/test_prop15674.py
```

The generated record is `evidence/e1_gmin_m4_prop15674.json`. The next
infinity-present shell has `p` finite boundary points and different endpoint
floors. The large all-finite range, general residual (ii), R1, global QVAR,
Type I, and the limit remain open.
