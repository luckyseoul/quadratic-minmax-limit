# A nonaffine first-defect-shell point at `p=11`

Date: 2026-09-04

Status: **proved exact counterexample to the proposed affine-shell
classification; residual (ii) remains OPEN**.

This note settles one sharply delimited question left open after Proposition
15.763.  A Paley integral `+p` eigenvector with one coordinate equal to `3`
need not come from a union of parallel affine lines.  The counterexample is
an exact adaptation of a published special-directions set, not a small-prime
census.

## Primary construction

Kiss and Somlai construct small sets in `F_p^2` with four special directions,
where a direction is special when the set is not equidistributed over its
parallel lines:

- Gergely Kiss and Gabor Somlai, *Special directions on the finite affine
  plane*, Designs, Codes and Cryptography 92 (2024), 2587--2597,
  <https://doi.org/10.1007/s10623-024-01404-y>;
- arXiv version <https://arxiv.org/abs/2109.13992>, Section 6, Figure 4.

The exact right-hand `11 by 11` Boolean matrix in their Figure 4 is embedded
as `KISS_SOMLAI_ROWS` in `src/e1_gmin_m4_prop15765.py`.  Its rows are indexed
by `y` and columns by `x`.  It defines a 33-point set `E0` with exactly the
four special spatial directions

```text
infinity, 0, 1, -1.
```

The row `y=1` is empty.  Adjoin that full horizontal line, disjointly, and
apply

```text
T(x,y)=(x,2y).
```

The resulting set `E` has size `33+11=44` and exactly the four special
directions

```text
infinity, 0, 2, -2.
```

The primary downloads used for transcription were pinned as follows:

```text
arXiv 2109.13992v3 PDF     e35a184f64f4e7af03744a10a8cd3eed3ecf05308f4d2dda286136810e85df0b
arXiv 2109.13992v3 source  ebbf527ea4e16e42416de3239e8be5fb7635a1ca2458117bcfe2fd7f3a77e525
```

## Paley realization

Write

```text
F_121=F_11[a]/(a^2-2).
```

Since `2` is a nonsquare modulo 11, this is a field.  The quadratic character
of `x+y*a` is the Legendre symbol of its norm `x^2-2y^2`.  The four direction
representatives

```text
(0,1), (1,0), (1,2), (1,-2)
```

all have character `+1`.  Thus all nonconstant Fourier directions of `1_E`
belong to the positive Paley eigenspace.  The tracked verifier also checks
the stronger pointwise integer identity, without floating point or a solver:

```text
Q 1_E = 11 1_E - 4 1,
```

where `Q_uv=chi(v-u)`.  For `D=F_11^2 minus E`, `|D|=77`, and hence

```text
Q 1_D = 11 1_D - 7 1,
3 + 2 sum_{v in D} chi(v-u) = 11(2 1_D(u)-1)   for every u.
```

Equivalently, `D` is a positive intriguing set of the Paley graph
`SRG(121,60,29,30)` with quotient matrix

```text
       D   E
D     40  20
E     35  25
```

and nonprincipal quotient eigenvalue `5=(11-1)/2`.

For the normalized conference matrix

```text
C = [0  1^T]
    [1   Q ],
```

put `y_infinity=3` and `y_u=2 1_D(u)-1`.  The pointwise identity and
`sum_u y_u=33` give

```text
Cy=11y,             ||y||^2=130,
```

with infinity the unique coordinate of magnitude three.  Replacing `3` by
`1` gives a Boolean vector `x`.  Direct exact multiplication yields

```text
Phi=11*122/2=671,   q_C(x)=649,   Phi-q_C(x)=22=2p.
```

This is therefore a genuine nonaffine point on the first Boolean defect
shell.

## Nonaffineness certificate

A nonempty proper union of parallel lines has exactly one special direction:
its own parallel class has a `0/11` profile and every transverse class is
constant.  The exact profiles of `E` are nonconstant in four directions and
constant `4` in the other eight.  In particular, no profile consists of four
entries `11` and seven entries `0`.  Thus `E` is not a union of four parallel
lines, and its complement `D` is not a union of seven parallel lines.

Compact hashes, with coordinates sorted as `(x,y)` and encoded as consecutive
unsigned bytes, are emitted in `evidence/e1_gmin_m4_prop15765.json`.

## Literature boundary

This same object is an `(n-m)`-perfect set in the block graph of
`OA(m,n)=OA(6,11)`, with `n-m=5`.  The closest primary classification source
is Bailey, Cameron, Gavrilyuk, and Goryainov, *Equitable partitions of
Latin-square graphs*, Journal of Combinatorial Designs 27 (2019), 142--160,
<https://arxiv.org/abs/1802.01001>.  Their Theorem 5.4 classifies the
one-Latin-square case `m=3`.  Section 8 explicitly presents classification
of `(n-m)`-perfect sets for several MOLS as a harder prospective problem and
gives coordinate-line and inflation examples; it does not supply the needed
`m=6` theorem.

A tempting later result also does not apply.  Theorem 7.2 of Dikstein, Dinur,
Filmus, and Harsha, *Boolean Function Analysis on High-Dimensional
Expanders*, Combinatorica 44 (2024), 563--620,
<https://doi.org/10.1007/s00493-024-00084-5>, classifies Boolean degree-one
functions only on **proper** simplicial complexes.  In the natural net
complex, vertices are the `6*11=66` affine lines from the selected directions
and the 121 top faces are the six lines through each point.  Its codimension-
one faces are the 6 omissions from each point-face, all distinct, so the top
up-map has domain dimension `6*121=726` and codomain dimension `121`.  It
cannot be injective.  The complex is not proper, exactly where that theorem
requires properness.

Targeted exact-phrase searches and a citation check of the 2019 paper found
no later primary classification of `(n-m)`-perfect sets in general MOLS/OA
block graphs.  More importantly, Proposition 15.765 gives a direct
counterexample to the particular parallel-line classification needed here.

## Scope and replay

This does **not** construct a switching set `H`, align different deletion
representatives, or meet the all-deletions hypotheses of residual (ii).  It
only proves that the full first shell is larger than the affine-alias family,
so Proposition 15.763 cannot be made universal through the proposed
classification.

Replay:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python src/e1_gmin_m4_prop15765.py
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q tests/test_prop15765.py
```

Residual (ii), E1, `L=1/2`, and the original MathOverflow limit remain
**OPEN**.
