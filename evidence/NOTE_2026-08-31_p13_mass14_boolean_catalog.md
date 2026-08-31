# Proposition 15.738: exact `p=13` mass-14 Boolean residual catalog

**Result status:** exhaustive finite certificate.

**Changed premise.**  Correct coefficient comparison in the exceptional
`p=13,t=3` row gives hard parallel count `P=2` or `P=8`, rather than the
previously recorded `P=5`.  The minimum opposite mass-14 cells consequently
occur at `Q=6` and `Q=0`.  This note classifies those cells.  It does not
repeat a finite-prime profile census.

## 1. Exact quadratic evaluation space

Let `Omega=J(13,7)`, of size 1,716.  The 78 pair monomials span every
degree-at-most-two function on the slice because

```text
sum_(j != i) x_i*x_j = 6*x_i,
sum_(i < j) x_i*x_j = 21.
```

Their `1716 x 78` evaluation matrix has rank 78 modulo 101.  For a fixed
4-set and three disjoint swap pairs, the alternating sum on the resulting
eight 7-sets annihilates every pair monomial.  Deterministic modular
elimination examined 62,721 of the 900,900 such cubes and retained 1,638
independent identities.  Their SHA-256 digest is

```text
ee92d6662f0f14523dc4c6620f89b407a66048dd4a6c0962dd9b058800136083
```

The modular rank gives real rank at least 1,638; annihilation of the
78-dimensional evaluation space gives real rank at most `1716-78=1638`.
Thus the common real nullspace is exactly the quadratic evaluation space.

## 2. The missing lift-to-Boolean bridge

For a nonnegative integral quadratic `B` with `4p E[B]=14`, Proposition
15.688 gives, when `H=max B>=2`,

```text
14 >= 2(p+1)-4H = 28-4H,
14 >= 3H.
```

Therefore `H=1` or `H=4`; height one is Boolean.  The remaining height-four
case is excluded separately at `Q=0` and `Q=6` by an exact necessary
residual model.

Write `W_st=epsilon*K_st` for the signed selected-edge coefficient between
fibres.  Since `b=0`, `A=2B`, and on a slice point `X`,

```text
4B(X) = Q-3 + sum_(s<t) W_st - 2 cut_W(X).
```

Averaging fixes

```text
sum W = 13Q-53.
```

The empty odd-fibre boundary makes every coefficient row sum even.  There
are `59-Q` transverse selected edges, so

```text
sum |W_st| <= 59-Q.
```

These are only necessary conditions, which is the safe direction for an
exclusion.  The symmetric group is transitive on `J(13,7)`, so if `H=4` a
maximum set may be moved to the first lexicographic 7-set; the model imposes
`B(X_0)=4` without loss of generality.

Both integer models have 1,811 constraints and were solved with OR-Tools
9.15.6755, presolve, 32 workers, and seed zero:

| `Q` | `sum W` | l1 budget | model SHA-256 | status | wall time |
|---:|---:|---:|---|---|---:|
| 0 | -53 | 59 | `6398ae7282d5bbc95527c1e3f6e80411017c75cf017d012c46533ced933ba2c1` | `INFEASIBLE` | 11.040936858 s |
| 6 | 25 | 53 | `990bee74e7a978df1b8a8f6ed28056849a4609183a28af95c24dfe6831dda2a2` | `INFEASIBLE` | 16.271538156 s |

Timing and parallel-search counters are environment-dependent.  Exact
infeasibility, the model hashes, and the structural counts are the
certificate claims.  Consequently a mass-14 residual cell at either
parallel count is Boolean and has support

```text
1716 * E[B] = 1716 * 7/26 = 462.
```

## 3. Exhaustive support-462 Boolean catalog

There are three explicit families:

```text
78   selected pairs:              x_i*x_j;
156  oriented mixed pairs:        x_i*(1-x_j);
858  signed all-equal triples:     x_i*x_j+x_k-x_i*x_k-x_j*x_k,
                                   equivalently z_i=z_j=-z_k.
```

All 1,092 supports are distinct, each has size 462, and each satisfies all
1,638 third-difference identities.  The full catalog digest is

```text
1609545bd2cddaa5f2389ea0e62b32a6bf62fd750bfb038aa1e3e1ba3ce127f6
```

For exact no-good exclusion, a support point can be moved to the first
slice point.  The catalog is invariant under `S_13`, so after imposing
`f(X_0)=1` it is sufficient to exclude the 294 catalog supports containing
`X_0`: 21 selected pairs, 42 oriented pairs, and 231 signed triples.  Their
digest is

```text
3d723e4171e711c8e8bf4d819edd0cee5a77eed77044466fa58135d0a6e04270
```

The exact model has 1,716 Boolean variables, 1,638 third-difference
equalities, one support equation, the support-point anchor, and 294 filtered
no-goods: 1,934 constraints in total.  Its text-proto SHA-256 is

```text
4b73bd641d500cdf6a0c5edb7f4c8b225903db3a81d13be7e2329df3dbdaed83
```

OR-Tools 9.15.6755 with presolve, 32 workers, and seed zero returned exact
status `INFEASIBLE` in 58.502218189 seconds (1,711,289 branches and 324,869
conflicts).  Therefore the three displayed families exhaust the support-462
Boolean quadratics on `J(13,7)`.

## 4. Coefficient offsets and moments

With signed-target normalization `epsilon*S_H=3+4B`, the families are

```text
x_i*x_j:        4+z_i+z_j+z_i*z_j,                 offset 6;
x_i*(1-x_j):    4+z_i-z_j-z_i*z_j,                 offset 4;
signed triple:  4+z_i*z_j-z_i*z_k-z_j*z_k,         offset 4.
```

Coefficient comparison modulo `q=6` leaves only `x_i*x_j` at both `Q=0`
and `Q=6`.

For that selected-pair target, the slice-kernel scalar is `-1/2` at `Q=0`
and zero at `Q=6`.  The normalized pair patterns differ by the all-one
complete-graph pattern.  Its degree-two and degree-four moments vanish over
`F_13`, so both parallel counts give, for all 78 pairs,

```text
M_2(i,j) = (i-j)^2,
M_4(i,j) = (i-j)^4.
```

Both values are nonzero for `i!=j`.  Proposition 15.738 stops at this exact
cell classification; the cross-direction binary-form argument using these
moments belongs to Proposition 15.739.  No claim that residual (ii) or the
quadratic-minmax limit is closed is made here.

Reproduce with:

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15738.py
PYTHONPATH=src python -m pytest -q tests/test_prop15738.py
```
