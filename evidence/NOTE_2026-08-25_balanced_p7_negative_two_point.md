# Balanced `p=7` negative two-point exclusion

Proposition 15.647 leaves one balanced negative-product profile at `p=7`:

```text
positive/negative baseline parallel counts = (3,3)
exceptional parallel counts                = (4,4)
finite edges E                             = 26
infinity edges I                           = 3
```

There is one exceptional direction of each quadratic type.  For either one,
the exact affine slack is

```text
A(X) = 1 - x_j + 2 B(X),       X in J(7,4),
```

where `B` is nonnegative, integer-valued, and quadratic.  Its exceptional
excess is eight, so `E B=2/7` and therefore `sum_X B(X)=10` over the 35
points of `J(7,4)`.

## Exact lift classification

Write `B(X)=sum_{ij subset X} c_ij`, and set

```text
U_i  = sum_{X contains i} B(X),
T_ij = sum_{X contains i,j} B(X).
```

The pair-incidence matrix of `J(7,4)` has rank 21.  Inverting its Gram
matrix, using `sum_X B(X)=10`, gives

```text
6 c_ij = 2 T_ij - U_i - U_j + 6,
6 B(X) = 2 sum_{ij subset X} T_ij - 3 sum_{i in X} U_i + 36.    (1)
```

Equation (1) is an exact bounded-integer description of all relevant lift
vectors.  Exhaustive CP-SAT enumeration leaves only:

| positive values of `B` | support | labelled vectors |
|---|---:|---:|
| `2,2,2,2,2` | 5 | 56 |
| six `1`s and two `2`s | 8 | 280 |
| eight `1`s and one `2` | 9 | 420 |
| ten `1`s | 10 | 1008 |

Thus every value is at most two, and there are exactly 1764 labelled lift
vectors.

## Infinity-star exhaustion

Normalize the finite boundary point to zero.  The three infinity edges form
a three-subset `S` of `F_49`.  In every baseline direction,

```text
K_st = eps_d (1 - w_s - w_t),
w_s  = |S intersect fibre_s| + 1_{s=0}.
```

Only 23 finite edges are transverse to that direction, so
`sum_{s<t}|K_st|<=23`.  For each of the two type-preserving exceptional-pair
orbits, this removes 210 of `C(49,3)=18424` stars.  The remaining 18214
stars split into 3038 orbits under the six square-semilinear maps fixing the
exception pair.

For every representative, the exact finite model imposes:

- all eight directional parallel-edge counts;
- the six baseline inter-fibre matrices;
- all 70 exceptional pointwise score identities;
- finite-boundary XOR and negative edge-product parity; and
- the exact quadratic-lift classification above.

The main 64-way sweep certified 6049 representatives infeasible and left 27
timeouts.  A second sweep constrained each exceptional lift directly to the
complete 1764-vector table and certified all 27 infeasible.  There are no
feasible or unknown cases.

The raw 77-file certificate archive is:

```text
/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-25-balanced-p7/
  balanced_p7_certificate_2026-08-25.tar.gz
SHA256 3aaf86364daf20ff3727382d49e9a6600a1b1239c3e2b827ab4f249c1b4b3f62
```

Its audited per-file manifest digest is
`2c269d29eba9434a05d5628b0055f0fa101a87b15fcb3bbd2c1deb28080f202d`.

Therefore the balanced profile is empty.  Together with Proposition 15.648,
every `p=7` negative two-point profile is excluded.  This does not address
the `p=5` case, positive-product finite cases, other boundary profiles,
residual (ii) as a whole, or R1.
