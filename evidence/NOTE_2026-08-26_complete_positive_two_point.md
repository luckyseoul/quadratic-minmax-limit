# Complete positive-product infinity-plus-point exclusion

Proposition 15.651 closes `D={infinity,v}`, `c_H=+1` for every odd prime
`p>=5`.  Together with Proposition 15.650, both edge-product signs of this
boundary shape are now closed at every prime in scope.  Other boundary
shapes, residual (ii), R1, global QVAR, and the limit remain open.

## 1. Exact coefficient model

Propositions 15.642--15.643 give, in every affine direction `d`,

```text
epsilon_d S_H(z) = 4 + z_j,
```

where `j` is the fibre containing the finite boundary point.  Translate that
point to zero.  Put `q=(p-1)/2`; let `I` be the infinity-edge count, `P_d`
the parallel finite-edge count, `n_s` the infinity-star count in fibre `s`,
and `K_st` the signed finite-edge count between fibres `s,t`.  Exact
coefficient comparison gives

```text
I   = 5 + q k0,
P_d = q kd,
sum_d kd = 8-k0,
K_st = epsilon_d (k0+kd+delta_sj+delta_tj-n_s-n_t).       (1)
```

Conversely, (1) reconstructs `epsilon_d S_H(z)=4+z_j` on the middle slice,
so no affine score information is lost.  The finite model selects `H`
directly, with exactly `4p+1` edges, boundary `{infinity,0}`, Paley-sign
product `+1`, and every equation in (1).

The exact directional edge budget also imposes

```text
sum_{s<t} |K_st| <= q(8-k0-kd).                           (2)
```

Enumerating the unordered integer fibre-count profiles in (2) is tiny and
is built into the certificate model as an exact allowed-assignment table.

Infinity parity and the aggregate form of (2) leave exactly

```text
p=5:  k0 = 0,1,2,3,4,5,8
p=7:  k0 = 0,2,4,8
p=11: k0 = 0,2,8
p=13: k0 = 0,1,8.
```

## 2. Immediate exact-l1 eliminations

The full fibre-count form of (2) is already decisive in three cases:

- at `p=7,k0=4`, only `kd=0` survives, but `sum kd=4`;
- at `p=11,k0=2`, only `kd=0` survives, but `sum kd=6`;
- at `p=13,k0=1`, only `kd=1` survives in each of 14 directions, but the
  required total is seven.

The `k0=8` endpoint has no finite edges, so its finite boundary is the whole
infinity star rather than the singleton `{0}`.  It is impossible directly.

## 3. The five-star `k0=0` structure

At `k0=0`, `I=5` and there are `4(p-1)` finite edges.  If `kd=0` in one
direction, the triangle lower bound in (2) is exactly the full transverse
edge budget.  Equality forces every inter-fibre term to have one sign.  For
`p>=7`, the special fibre contains a star point, and every finite edge has
quadratic type opposite to the unpopulated direction.

Thus all unpopulated directions have one type and every populated direction
has the other.  The populated multiplicities total eight.  One type has
only six directions at `p=11` and seven at `p=13`, immediately excluding
both primes.

At `p=7` there are two and only two possibilities:

```text
all-one:   kd=1 in all eight directions;
type-split: kd=2 in all four directions of one type,
            kd=0 in all four directions of the other type.
```

## 4. Exact `p=7` exhaustion

For the type split, exact `l1` equality reduces the star search to 238644
five-sets for each populated type.  Exactly 2250 survive all eight fibre
profiles.  Every survivor contains zero.  Square field multiplication and
Frobenius form a 48-element stabilizer and reduce the survivors to 56
orbits per type.  The exact fixed-star edge model certifies all 112 orbit
representatives `INFEASIBLE`, with no feasible or unknown case.

For the all-one profile, square multiplication gives three exhaustive
normalizations:

1. zero is absent from the star;
2. zero is present and some other star point is square, normalized to `1`;
3. zero is present and all four other points are nonsquare, with one
   normalized to `8`.

All three exact models are `INFEASIBLE`.

## 5. Remaining finite cases

The direct additive-coefficient model certifies all seven `p=5` arithmetic
cases infeasible.  It also independently certifies the remaining nonzero
`k0` cases used in the audit.  In total the audit checks 14 direct cases,
112 rigid `p=7` star-orbit cases, and three normalized all-one cases, with
zero feasible and zero unknown result.

Proposition 15.643 supplies every `p>=17`.  Therefore the positive-product
infinity-plus-point branch is closed for every odd `p>=5`.  Proposition
15.650 already closes the negative-product branch, so this entire boundary
shape is now gone.

## 6. Reproduction and archive

Core scripts:

- `scripts/positive_two_point_additive_cpsat.py`;
- `scripts/p7_positive_star_classify.py`;
- `scripts/p7_positive_fixed_star_cpsat.py`;
- `scripts/p7_positive_orbit_certificate.py`;
- `scripts/positive_two_point_certificate_audit.py`.

Permanent raw archive:

```text
/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-26-positive-small/
  positive_two_point_certificate_2026-08-26.tar.gz
SHA256 a507cb917f97a998638e9d4226fb6925a2c358c8c9212fc25d1ea17795cffd26
```

The audit JSON has SHA256

```text
6eb06210b8837b929c52c3b157ec00aba4defe5f7a5399330e47b1a43b8adee8
```

and the 112-row `p=7` orbit certificate has SHA256

```text
d9e01ab4bc37ba77f620111852e833e0ca37fa4808db8a56a1cbe0e05d93614c
```
