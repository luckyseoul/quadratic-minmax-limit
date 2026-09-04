# Sporadic Peisert `G(23^2,2)`: exact OA/PN deduplication record

Date: 2026-09-04

Status: **regularizable; not a residual-(ii) closure**.

The sporadic rank-3 conference graph on 529 vertices is genuinely
nonisomorphic to the Paley and ordinary Peisert graphs, but it is not outside
the balanced projective-direction / linear-OA family.  Therefore it is not a
new target for the order-square conference route.

## Primary construction and exact reconstruction

Lim's Definition 5.2.6 and Table 5.4, line 11, construct
`G(23^2,2)=Cay(F_23^2,S)` with `S` an orbit of
`M_0 = Z(GL(2,23)) o SL(2,3)` of order 264:

- T. K. Lim, *Edge-Transitive Homogeneous Factorisations of Complete
  Graphs*, arXiv:math/0605253, pp. 60--65:
  <https://arxiv.org/abs/math/0605253>.
- The downloadable published adjacency matrix is at
  <https://www.math.mun.ca/distanceregular/graphdata/sporadicpeisert529.am.csv>.

`src/sporadic_peisert529_exact.py` realizes the group over `F_23` from the
three matrices

```text
[ 0 -1]   [ 8 16]   [2 0]
[ 1  0],  [17 14],  [0 2].
```

Their group has order 264, and the orbit of `(1,0)` has size 264.  The exact
translation check gives common-neighbor histogram `{131:264, 132:264}`, hence
`SRG(529,264,131,132)`.  For `S_graph=J-I-2A`, this identity and
`S_graph 1=0` give `S_graph^2=529I-J`; adjoining the normalized all-ones
row and column therefore gives a symmetric conference matrix of order 530
with `C^2=529I`.

A pynauty canonical-certificate comparison verifies that the constructed
matrix and the downloaded matrix are isomorphic.  They are differently
labelled, so their raw CSV hashes are intentionally recorded separately:

```text
published CSV    fcdd847709adc1527374781fb81857bf3ce741c538e35f59126b1c6256e3cda6
constructed CSV  eeba88c7441385065a235fecc62240de25474f22ad15cab3dcc307b2cbf0c3bc
canonical cert   64c08ce8c6cacedbf5201441a5230637b1961b19f43d784985b5a1d57332a3b0
Aut order         139656 = 23^2 * 264
```

## Explicit Boolean eigenshell witnesses

The scalar generator forces the connection orbit to be exactly the union of
the 12 projective directions

```text
infinity, 0, 2, 4, 5, 7, 9, 11, 12, 13, 17, 21.
```

With finite vertices lexicographically labelled `(a,b)` and `x_infinity=1`:

- `Cx=+23x`: put `x=-1` on the 11 parallel lines
  `a-b in {0,...,10}` and `x=+1` elsewhere.  The finite minus count is 253;
  the little-endian `int16` witness hash is
  `19706232e3f181513356b8515cb3897f79b3827d36975be0549ea704a9c260b7`.
- `Cx=-23x`: put `x=-1` on the 12 horizontal lines
  `-b in {0,...,11}` and `x=+1` elsewhere.  The finite minus count is 276;
  the witness hash is
  `1d9873cc4139bb507a16837e4bfeebac36906bd19c2270f90ca01f2076efc62e`.

Both identities are checked coordinate-by-coordinate by the tracked module.
This is the general linear-OA regularization mechanism in concrete form.

Replay:

```text
python src/sporadic_peisert529_exact.py
pytest -q -n 0 tests/test_sporadic_peisert529_exact.py
```

The compact machine-readable output is
`evidence/sporadic_peisert529_exact.json`.  Do not rerun CP-SAT on this class:
the displayed witnesses already settle both Boolean eigenshells exactly.

## Adjacent Mathon boundary

Mathon's pseudocyclic construction has descendant order
`v=(4t+1)(4t-1)^2`.  Since `gcd(4t+1,4t-1)=1`, a square `v=p^2` forces
`4t+1=r^2` and then `p=r(4t-1)=r(r^2-2)`.  For positive `t`, both factors
exceed one.  Thus the square-order members of that family never have prime
square root; its first case is `v=441=21^2`, not a residual-(ii) prime case.
