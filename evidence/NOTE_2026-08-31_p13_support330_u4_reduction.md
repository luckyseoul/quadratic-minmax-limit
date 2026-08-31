# Proposition 15.746: support-330 equality catalog and the `u=4` reduction

## Status and scope

This is an **exhaustive finite equality classification and proved open
reduction**.  It classifies the sharp mass-ten Boolean lifts on `J(13,7)` and
propagates that classification through the common `p=13,t=4,u=4` edge
ledger.  It does **not** close `u=4`, `u=6`, the full `p=13,k=60` row,
residual (ii), Type I, or the limit.  The exact `k=60` residue set remains
`u in {4,6}`.

## Pointwise bridge into the Boolean problem

Each of the seven hard `u=4` cells has `b=2` and scaled mean 22.  The
all-positive equality quadrature for the `b=2` floor first identifies the
XNOR baseline pointwise.  Only then is its difference a globally
nonnegative integral degree-at-most-two lift `B`, with

```text
4p E[B] = 10.
```

This pointwise step is essential: an averaged baseline value alone would not
satisfy Proposition 15.688's nonnegativity hypothesis.  Once the hypothesis
is live, the sharp mass-ten case of 15.688 forces height one.  Thus `B` is
Boolean and has support

```text
(10/52) * C(13,7) = 330.
```

## Exact support-330 classification

Proposition 15.738 supplies the exact quadratic evaluation space on
`J(13,7)`: rank 78, with nullspace cut out by 1,638 independent
third-difference identities.  The two proposed families are

| family | polynomial | full count | anchored count |
|---|---|---:|---:|
| omitted pair | `(1-x_i)(1-x_j)` | 78 | 15 |
| all-equal triple | `1-x_i-x_j-x_k+x_i*x_j+x_i*x_k+x_j*x_k` | 286 | 55 |

Every one of the 364 supports is distinct, has size 330, and satisfies all
1,638 identities.  The runner checks closure under the twelve adjacent
coordinate transpositions generating `S_13`, and it checks that the fixed
anchor point has full orbit of size 1,716.  The one-point anchor is therefore
an executable symmetry reduction, not an asserted normalization.

After anchoring, the exact no-new-support model has:

```text
1,716 Boolean variables;
1,638 third-difference equalities;
1 support equation;
1 anchor equation;
70 anchored candidate no-goods;
1,710 constraints total.
```

The unsharded exact CP-SAT model returns `INFEASIBLE`; an independent
deterministic one-worker replay returns the same exact status.  Therefore the
364 displayed supports exhaust all support-330 Boolean quadratics.  The two
catalog digests are

```text
full:     4edf1fe1b9c73f05598b667dba121f064807c68421a4df2c8db7090a3e3ff35f
anchored: 84ce6099dcca66f7cc2792dc60bcbb378672f2e9cac2b19e02812f2f20563c7a
```

The optional GPU path checks only the proposed catalog's support sizes and
identity residuals.  It is explicitly not a proof premise; exact CPU integer
checks and exact CP-SAT infeasibility are authoritative.

## Consequence for the seven hard directions

With `z_i=2x_i-1`, a hard cell has target

```text
4 + z_i*z_j + 4B.
```

For an omitted pair,

```text
4B = 1-z_a-z_b+z_a*z_b,
```

so the full target has coefficient offset three.  For an all-equal triple,

```text
4B = 1+z_a*z_b+z_a*z_c+z_b*z_c,
```

so the full target has offset five.  Common `hT`, the isolated-chart
coefficient congruence, and `7P<=61` force

| hard family | common `P` | `hT=14P-61` |
|---|---:|---:|
| omitted pair | 3 | -19 |
| all-equal triple | 5 | 9 |

The hard families cannot mix.  For the seven opposite counts `Q_L`, put

```text
e_L = P + Q_L - 8 >= 0.
```

Since `sum Q_L=61-7P`, one has `sum e_L=5`, and the opposite scaled mean is

```text
a_L = 12 + 14e_L.
```

At least two opposite directions therefore have mean 12.  Their minimum
counts are `Q=5` in the `P=3` branch and `Q=3` in the `P=5` branch.

A phase-zero mean-12 cell is either the exact `b=12` literal
`3+2(1-x_j)=4-z_j`, or a `b=0` cell `A=2C` with `4p E[C]=12`.
Proposition 15.688 leaves exactly height one or height four in the latter
case; height one is Boolean of support 396.  In the `P=3` branch, the literal
would require `6 | Q-3`, which fails at `Q=5`.  Thus each forced mean-12 cell
there is a `b=0` mass-12 lift, with height one/support 396 or height four.
The `P=5,Q=3` branch still permits either the literal or the `b=0` lift.

## Degree-six common constraint in the omitted-pair branch

For each omitted-pair hard cell and `r=1,2,3`, the normalized even moments
have the two-power-sum form

```text
h*M_(2r) = alpha^(2r) + beta^(2r).
```

Newton's identity gives the homogeneous binary sextic

```text
F6 = 2h*M6 + h*M2^3 - 3*M2*M4 = 0.
```

The complete-graph and star gauge sums vanish in degrees 2, 4, and 6.  The
certificate directly checks all `78^2=6,084` possible overlaps of the
baseline pair and omitted lift pair.  Seven hard directions are therefore
seven projective roots of `F6`, so `F6` vanishes identically.  Every forced
`P=3,Q=5,b=0` opposite mass-12 cell must satisfy `F6=0`.  In the `P=5`
all-equal-triple branch, all 22,308 baseline-pair/triple patterns give ranks
`1,2,3` to the weighted feature lists `N2`, `(N4,N2^2)`, and
`(N6,N2*N4,N2^3)`.  Hence no nonzero universal weighted-homogeneous
polynomial identity in these even moments exists through degree six.  This
does not rule out a different invariant.

The next model must account for the opposite sign.  With its natural
normalization `N'_(2r)=(-h)M_(2r)`, the executable constraint is

```text
2N'_6 + (N'_2)^3 + 3N'_2*N'_4 = 0,
```

not the hard-sign Newton formula copied unchanged.

## Next exact gate

Do not launch a broad mass-12 or directional-profile census.  The immediate
finite gate is one `P=3,Q=5,b=0` mass-12 cell under the additional exact
constraint `2N'_6+(N'_2)^3+3N'_2*N'_4=0`, split into:

- height one, hence a support-396 Boolean quadratic; and
- height four.

An exclusion or a classification strong enough to couple those cells would
strictly advance `u=4`.  A standalone support count without `Q=5` and `F6`
would not.

## Reproduction

```bash
PYTHONPATH=src python scripts/p13_support330_boolean_classifier.py --workers 32
PYTHONPATH=src python scripts/p13_support330_boolean_classifier.py \
  --workers 1 --output /tmp/p13_support330_boolean_classifier_one_worker.json
PYTHONPATH=src python src/e1_gmin_m4_prop15746.py
PYTHONPATH=src pytest -q -n 0 \
  tests/test_p13_support330_boolean_classifier.py tests/test_prop15746.py
```

Tracked generated evidence:

- `evidence/p13_support330_boolean_classifier.json`
- `evidence/e1_gmin_m4_prop15746.json`
