# Proposition 15.736: exact `p=11` sharp Boolean-quadratic catalog

**Result status:** exhaustive finite certificate.

**Changed premise.**  The earlier `p=11` reduction used a conditional
restriction--extension input whose source could not be checked locally.  This
certificate instead works directly on all 462 points of `J(11,6)`.  It does
not rerun a finite-prime profile census.

## Exact linear certificate

On `J(11,6)`, the 55 pair monomials span every function of degree at most two:

```text
sum_(j != i) x_i*x_j = 5*x_i,
sum_(i < j) x_i*x_j = 15.
```

Their `462 x 55` evaluation matrix has rank 55 modulo 101.  For a fixed
3-set and three disjoint swap pairs, the alternating sum on the resulting
eight 6-sets kills every pair monomial.  Deterministic modular elimination
examined 8,321 such cube rows and retained 407 independent rows.  Their
SHA-256 digest is

```text
6e17bd62f6ee15bf06065bdadfeeba9e4c4a8f79302c753214b7185ba9b47777
```

The modular rank 407 gives real rank at least 407.  Annihilation of the
55-dimensional quadratic evaluation space gives real rank at most
`462-55=407`.  Hence the real nullspace of the selected identities is exactly
the quadratic evaluation space.  No unproved extension theorem enters.

## Exact Boolean classification

First, Proposition 15.688 bridges the residual equality object to this model.
At `p=11`, equality has scaled mass `4p E[B]=8`, whereas every nonnegative
integral quadratic with `H=max B>=2` has scaled mass at least 12.  Thus a
nonzero equality lift has `H=1` and is Boolean.  Its mean is `2/11`, and the
462-point slice therefore gives support size 84.  The source and tests check
each of these exact values; Booleanity is not an implicit assumption.

The CP-SAT model has:

- 462 Boolean variables `f_X`;
- the 407 exact eight-term identities;
- `sum_X f_X=84`;
- 220 no-goods, one for each known sharp support: 55 omitted-pair forms and
  165 all-equal triple forms.

It has 628 constraints in total.  Its deterministic text-proto SHA-256 is

```text
0070bf67f0891acb502cd55446b7b4c7162188d2f219350a2dc00589fa5a8b04
```

Every listed support has size 84, the 220 supports are distinct, and code
checks that each satisfies all 407 identities.  Their catalog digest is

```text
6f9b55283e78540ec389c2674ebd1b8a93f4b179bca8e42cfd1e6b5f8f1b7535
```

The canonical run used OR-Tools 9.15.6755, presolve, 32 search workers, seed
0, and no orbit anchor.  It returned exact status `INFEASIBLE` in 9.011843513
seconds, with 1,883,810 branches and 107,304 conflicts.  Timing and parallel
search counters are environment-dependent; infeasibility, dimensions,
counts, and digests are the certificate claims.  Therefore every Boolean
quadratic on `J(11,6)` with support 84 is one of the 220 listed forms.

## Residual consequence and live obstruction

The two signed targets and coefficient offsets are

```text
omitted pair:     4-z_i-z_j+z_i*z_j,                 offset 2;
all-equal triple: 4+z_i*z_j+z_i*z_k+z_j*z_k,         offset 4.
```

For `q=5`, the hard-`b=2` branch forces a minimum opposite count `Q=3`.
Neither offset is congruent to 3, so this branch is excluded.  The
hard-`b=p-1` branch forces `Q=4`: omitted pairs are excluded, but all-equal
triples survive.  Its opposite count sum is 26 across six directions, so at
least four minimum directions must carry all-equal triple targets.

This certificate does **not** prove those simultaneous targets compatible or
incompatible.  Thus `p=11`, residual (ii), multi-level Type I, and the final
quadratic-minmax limit remain open.

Reproduce with:

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15736.py
PYTHONPATH=src python -m pytest -q tests/test_prop15736.py
```
