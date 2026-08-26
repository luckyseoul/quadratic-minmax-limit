# `p=7` unsaturated four-finite exclusion modulo seven

Proposition 15.655 closes the unsaturated part of the `p=7` four-finite
boundary branch.  Combined with Proposition 15.654, every four-finite
boundary at `p=7` is excluded for both Paley edge-product signs.  Combined
also with Proposition 15.653, every size-four boundary at `p=7` is closed.

This does **not** close residual (ii): every `p=5` size-four shape and
boundaries of size at least six remain open.  R1, global QVAR, Type I, and
the limit remain open.

## 1. The common exact linear system

Fix a four-finite boundary, a product sign, and the elevated direction in
each unsaturated quadratic type.  Let `x_e` be the indicator of an edge in
the 29-edge graph `H`.  There are `C(50,2)=1225` edge variables.

For each of the eight affine directions and each `X in J(7,4)`, let
`N(d,X)` be the edges whose normalized sign is bad.  If `A_d(X)` is the
exact directional slack, then

```text
epsilon_d S_d(X) = 29 - 2 sum_{e in N(d,X)} x_e
                   = 3 + 2 A_d(X),
```

so every candidate satisfies the integer equation

```text
sum_{e in N(d,X)} x_e = 13 - A_d(X).                 (1)
```

The 280 equations (1), `sum_e x_e=29`, and the fixed distinguished edge
`x_(0,1)=1` form a common `282 x 1225` zero-one coefficient matrix `M`.
Exact elimination over `F_7` gives

```text
rank_F7(M) = 147,
left-null dimension = 282 - 147 = 135.
```

Thus every complete slack-catalog right side must satisfy 135 exact
mod-seven dependencies.  Failure of even one dependency is a rigorous
integer infeasibility certificate; no SAT timeout or numerical tolerance is
involved.

## 2. Complete catalog syndrome joins

The orbit source from Proposition 15.654 contains 518 unsaturated
square-semilinear orbits whose sizes sum to 23,520 boundaries for one
product sign.  Fixing elevated directions yields 2,408 cases.  Every case
has at most two non-singleton catalogs, with the following complete pattern
census:

```text
catalog sizes       cases
1764                 1372
1764 x 36             294
1764 x 1764           112
2233 x 36             294
2233 x 1764           336
                     ----
                     2408
```

For one catalog, the calculation tests whether any syndrome equals the
negative singleton contribution.  For two catalogs, it hashes the first
catalog's 135-coordinate syndrome vectors and performs an exact complementary
lookup for every row of the second.  This counts all tuples without expanding
the largest `2233*1764` Cartesian product.

The five pattern classes contain 1,716,742,440 catalog tuples in total.
The result is

```text
fixed elevation cases             2408
mod-7 infeasible                   2408
surviving cases                       0
missing / duplicate cases          0 / 0
```

## 3. Independent audit

The production certificate and audit deliberately use different routes.

The production program:

- builds equation rows edge by edge;
- obtains a left-null basis by column-pivot Gauss--Jordan elimination;
- relabels complete slack-value catalogs directly; and
- computes each right side as `13-A`.

The audit program:

- rebuilds rows from vectorized normalized sign products;
- discovers dependencies by incremental row-span reduction with witnesses;
- reconstructs catalog scores from independently interpolated
  constant-plus-pair coefficients; and
- independently reconstructs all 518 orbit and 2,408 elevation keys.

It reports rank 147, dependency dimension 135, 2,408 certificate cases,
zero missing or extra cases, and zero recomputed compatible cases.  A third
spot check with SymPy's finite-field `DomainMatrix` gives rank 147 for `M`
and rank 148 for the augmented hard orbit-145/catalog-row-zero system.

## 4. Transfer to the other product sign

Proposition 15.654 already audits the required nonsquare anti-isometry.  It
fixes infinity, finite zero, and the distinguished edge; anti-commutes with
the Paley conference matrix; exchanges the two eigenshells; and preserves
normalized scores.  Since `|H|=29` is odd and a four-finite boundary gives
even infinity degree, it flips the Paley edge product.  The one-sign
mod-seven exclusion therefore transfers bijectively to the other sign.

## 5. Literature and OEIS context check

Targeted searches after the finding located the general literature on
Paley incidence codes, especially Ghinelli--Key, *Codes from incidence
matrices and line graphs of Paley graphs*, Adv. Math. Commun. 5 (2011),
93--108, DOI `10.3934/amc.2011.5.93`.  That work computes `p`-ary ranks of
ordinary Paley graph incidence matrices.  It does not contain this
`282 x 1225` affine-score matrix, its rank 147 over `F_7`, the Johnson slack
catalogs, or the 2,408-case syndrome exclusion.  Johnson-scheme and finite
incidence-code papers remain methodological context.

Exact web searches for the rank description and the tuple total found no
matching theorem.  Direct OEIS searches returned no entry for
`1716742440`, `3939012`, the pattern tuple
`1372,294,112,294,336`, or `282,1225,147,135`.  These checks are context
only; no integer-sequence novelty or priority claim is made.

## 6. Reproduction and permanent archive

Core programs:

- `scripts/p7_unsaturated_modular_catalog_filter.py`;
- `scripts/p7_unsaturated_mod7_batch.py`;
- `scripts/p7_unsaturated_mod7_audit.py`;
- `scripts/p7_unsaturated_slack_catalog.py`.

Permanent archive:

```text
/mnt/storage/e1work/quadratic-minmax-limit-finite/
  2026-08-26-p7-four-point/
  p7_no_infinity_unsaturated_mod7_certificate_2026-08-26.tar.gz
```

Hashes:

```text
archive  ad4efb2c1f3eb53f2ebd089c8b5b6ce959945bf7ce1298c7dfa694eaeda8965d
batch    64a87c509f92144e1dd1956c92145abc3b42deebb7fdb610405963b30d80c50f
audit    3348579efe686302ba5eb9e654d8121f8445e9ffece19047d581413964581350
source   7f7d3cc26077bb40ac096b638c6fc20ddf1a8fe6ddee60641f2fb568bacfd077
```
