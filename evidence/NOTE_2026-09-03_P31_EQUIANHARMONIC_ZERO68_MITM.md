# Exact $p=31$ equianharmonic zero-six/eight exclusion

## Scope and theorem

This note closes one finite fiber and nothing larger.

> **Finite theorem.**  At $p=31$, $r=b=7$, and $k=11$, no choice of
> six positive all-equal atoms and seven compact atoms realizes the constant
> tangent-conic edge-orbit word while making every degree-six and degree-eight
> moment channel zero.

Equivalently, the `p=31,b=7,k=11` constant equianharmonic conic fiber is
`UNSAT` after the seven congruences

\[
 F_{d,j}(s,t)=(s-t)^2(st)^j(s+t)^{d-2-2j},
 \qquad d\in\{6,8\},\quad 0\le j<d/2
\]

are imposed.  This does **not** exclude other equianharmonic primes or
parameters, does **not** construct or exclude the Boolean/common lift in
general, and does **not** close residual (ii) or $L$.

The executable certificate is
`evidence/p31_equianharmonic_zero68_mitm.cpp`; its aggregate is
`evidence/p31_equianharmonic_zero68_mitm_manifest.json`.  It is an independent
integer meet-in-the-middle/exact-cover generator and uses no SAT, CP-SAT, or
MIP library.

## 1. Alignment makes the search finite and small

Orient the 29 nonantipodal conic orbits by the target word
$t_e\in\{\pm1\}$.  For an atom column $c$, let

\[
 s(c)=\sum_e t_e c_e.
\]

A complete enumeration of the 4,495 label triples and 13,485 distinguished
label triples gives the following exact score distributions:

| atom | score counts |
|---|---|
| all-equal | `{-3:9,-2:1,-1:702,0:3071,1:702,2:1,3:9}` |
| compact | `{-2:111,-1:2133,0:8997,1:2133,2:111}` |

Thus an all-equal atom has score at most three and a compact atom has score at
most two.  Any realization would have total score

\[
 t\mathbin\cdot t=29,
\]

whereas the six/seven atom maximum is $6\cdot3+7\cdot2=32$.  Its total
deficit is therefore exactly three.  There are only three disjoint deficit
partitions:

\[
 (3),\qquad(2,1),\qquad(1,1,1).                 \tag{1}
\]

This is the completeness gate: atoms of deficit at least four cannot occur,
and the exhaustive generator treats every AE/compact assignment within (1),
including repeated atoms.

## 2. The maximal half is reconstructed, not searched

The nine score-three AE atoms are the nine disjoint order-three $q$-cycle
triangles, with $q=(1-k)/(1+k)=25$.  They cover 27 target coordinates and
leave the two edges of the broken cycle.  After the nonmaximal atoms and all
compact atoms have been selected, exact target equality forces the
multiplicity of each maximal AE cycle uniquely.

More explicitly, write the three target coordinates of cycle $C_c$ as
$e_{c,0},e_{c,1},e_{c,2}$, and let $x_e$ be the contribution already
selected.  Exact recovery requires

\[
 x_{e_{c,0}}=x_{e_{c,1}}=x_{e_{c,2}},\qquad
 m_c=1-x_{e_{c,0}}\ge0,                              \tag{2}
\]

and contribution one on each broken-cycle edge.  The required number of AE
cycles is also checked exactly.

If $M_c\in\mathbf F_{31}^7$ is the degree-six/eight moment vector of
$C_c$, the moment test is performed additively through

\[
 J=M_{\rm selected}-\sum_c x_{e_{c,0}}M_c,
 \qquad J+\sum_cM_c=0.                                \tag{3}
\]

Equations (2)--(3) are equivalent to direct target and moment replay; they do
not relax either condition.

The 111 score-two compact atoms have a second exact structure.  Five use an
invisible self-antipodal third edge.  Every other one has a single signed
off-target orbit, supported on 51 orbit pairs.  For a prescribed residual
off-vector, every maximal compact completion is therefore generated as:

1. the forced signed multiplicity on each nonzero residual orbit;
2. an arbitrary multiset of cancelling (+/-) pairs on the 51 supported
   orbits; and
3. an arbitrary multiset of the five antipodal-edge atoms.

This parameterization is exhaustive, including multiplicities.

## 3. Exact disjoint partition totals

The run exhausted the following mutually disjoint cases.

| deficit case | prefilter population | exact maximal completions | edge hits | zero-6/8 hits |
|---|---:|---:|---:|---:|
| one deficit-three atom | 5,204 atoms | 13,528,344 | 60 | **0** |
| deficit two plus deficit one | 20,697,666 ordered-by-deficit pairs; 79,918 off-compatible | 87,840,508 | 2,160 | **0** |
| three deficit-one atoms, at least one AE | 2,278,045 multisets; 24,828 off-compatible | 20,465,801 | 392 | **0** |
| three deficit-one compact atoms | 1,619,689,995 multisets; 2,027,542 unsupported-projection 3SUM hits; 1,089,526 completable | 108,480,057 | 14,464 | **0** |
| **total** | — | **230,314,710** | **17,076** | **0** |

The last line is not obtained by iterating 1.62 billion triples.  The 2,133
deficit-one compact atoms are projected to the complement of the 51 off-orbits
available to maximal compact atoms.  The 2,275,911 unordered atom pairs form
1,031,232 exact projection keys.  A third atom survives precisely when its
negative projection is a stored key.  Assigning a triple to its largest atom
index makes the 3SUM shards disjoint.  Only then is the supported residual
completed by the exact generator above and tested by (2)--(3).

The captured half-open shard coverage was:

- deficit `(2,1)`: `[0,3000)`, `[3000,5600)`, `[5600,8200)`, `[8200,9699)`;
- three deficit-one atoms with an AE atom: `[0,750)`, `[750,1450)`,
  `[1450,2133)`;
- three deficit-one compact atoms: `[0,800)`, `[800,1500)`, `[1500,2133)`.

For `(2,1)` the index is the unique deficit-two atom.  In `d1ae` it is the
smaller of the two compact indices in the one-AE/two-compact subtype; the
all-AE and two-AE subtypes occur only in the shard whose lower endpoint is
zero.  In `d1k` it is the largest of the three compact indices.  These
intervals cover each relevant exceptional multiset exactly once.  The
one-deficit-three case is a single complete run.

## 4. Reproduction

The recorded build used

```text
g++ (Ubuntu 15.2.0-16ubuntu1) 15.2.0
```

and the command

```bash
g++ -O3 -std=c++20 -Wall -Wextra -pedantic -pthread \
  evidence/p31_equianharmonic_zero68_mitm.cpp \
  -o /tmp/p31_equianharmonic_zero68_mitm
```

The aggregate header is cheap to replay:

```bash
/tmp/p31_equianharmonic_zero68_mitm manifest
```

A one-process exhaustive replay is:

```bash
/tmp/p31_equianharmonic_zero68_mitm d3
/tmp/p31_equianharmonic_zero68_mitm d2d1 0 9699
/tmp/p31_equianharmonic_zero68_mitm d1ae 0 2133
/tmp/p31_equianharmonic_zero68_mitm d1k 0 2133
```

Each ranged mode may instead be split into disjoint half-open subranges and
its final counters added.  The `d1ae` mode assigns its all-AE and two-AE
subcases only to the range whose lower endpoint is zero.

Pinned SHA-256 values:

```text
14e23138797b8bde9edbbd447c69ee735ec0a85600145b2ecf0988e490c59520  evidence/p31_equianharmonic_zero68_mitm.cpp
a74e57224ece6d4a4b8d609df40d052f8929ca091e238dc75b3d88bdef26a550  /tmp/p31_equianharmonic_zero68_mitm
```

The Python replay in `src/e1_gmin_m4_p31_equi_zero68_mitm.py` independently
reconstructs the full atom score census, the nine-cycle/51-off-orbit
structure, every combinatorial population in the table, and the explicit
OPEN flags.  `tests/test_p31_equi_zero68_mitm.py` pins both the reduction and
the exhaustive C++ source.

## 5. Consequence for the live proof

The edge-only witness in the equianharmonic branch remains a valid
counterexample to odd-zero centrality.  This certificate proves that its
entire `p=31,b=7,k=11` constant-conic fiber cannot also satisfy both even
moment systems.  It is meaningful evidence for how the conic survivor must
be attacked, but it is a finite fiber theorem.  Residual (ii) and (L) remain
**OPEN**.
