# Proposition 15.717: positive p7 infinity-plus-seven z=3 is closed

The exact pair-transversal census leaves 5,488 positive `z=3` boundaries in
ten orbits of the 2,352-element square affine-semilinear group.  Six orbits
have size 392 and four have size 784.

The corrected exact mean ledger has 400 leaves:

```text
360  residue 00
 20  residue 04
 20  residue 40
---
400
```

Because these are unpointed translation orbits, the distinguished-edge row
is omitted.  The resulting 281-by-1,225 integer edge system has rank 146
modulo seven and 135 exact left dependencies.  Complete Johnson-slice
catalog joins reject 398 leaves.  The only survivors are source orbits 92
and 100; each has three complete 1,764-row medium catalogs and exactly four
matching row-index triples modulo seven.  No high-direction relaxation is
used in either survivor.

`scripts/p7_infinity7_positive_z3_multimod_join.py` extracts all eight
triples with their exact catalog row indices.  For each one it reconstructs
the same 281-entry integer right side—edge count 29 followed by eight
35-entry bad-count blocks—and tests complete left-nullspace bases of that
same integer matrix modulo 3, 5, 7, and 11.  The ranks and dependency
dimensions are:

```text
modulus   rank   dependencies
   3       161       120
   5       167       114
   7       146       135
  11       167       114
```

Every extracted tuple passes the required mod-seven replay and every one
fails modulo three.  Thus no integral edge right side can realize either
surviving leaf.  All 5,488 actual positive `z=3` boundaries are excluded.

The positive remainder is now only 56 affine-line boundaries in two `z=7`
orbits.  The projected envelope falls from 212 profiles to the two labelled
line profiles.  The positive endpoint, negative branch, residual (ii), and
every top-level gate remain open.

The mod-seven decision certificate is
`c39fb7f530a6380c09d0bf300d6d249df2304370867066e6ce74de62906f275f`.
The eight extracted row-index tuples have certificate
`ababbe1e75ad4c913d3262c6df58f5ad44b8755889f5d3614978463ecb826a98`,
and their multimodular decisions have certificate
`c7e06a92a4581309b1b3e59a2e1e84f1d48e5988c53337864e73521f29debc82`.

Targeted GitHub-code, MathOverflow, literature, and OEIS searches found no
prior occurrence of this exact exclusion.  Numerical hits for 5,488 and
225,792 were unrelated and have no evidentiary role.
