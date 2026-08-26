# `p=7` saturated four-finite boundary exclusion

Proposition 15.654 closes the `p=7` boundary profiles consisting of four
finite points whose two Proposition 15.632 type costs are both 32.  It
closes both Paley edge-product signs.

At Proposition 15.654 this was a reduction, not completion: exactly 23,520
unsaturated boundaries, or 518 square-semilinear orbits, remained per sign.
They are subsequently closed by 15.655, and `p=5` size four by 15.656.
Larger boundaries, residual (ii), R1, and the limit remain open.

## 1. Boundary and slack classification

For four finite boundary points, each direction has `b=0,2`, or `4` odd
fibres.  At `p=7` the exact scaled floor table is

```text
                  phase 0   phase 1
b=0                   0        14
b=2                   8         6
b=4                   8        14.
```

Each quadratic direction type has exact budget 32.  Complete enumeration of
the `C(49,4)=211,876` finite four-sets leaves 82,320 boundary-only survivors
for either product sign.  Of these, 58,800 use exactly 32 in each type.

Budget saturation forces every directional slack to have minimum mean.  The
`b=0` and `b=2` slacks are pointwise forced.  Exact degree-two evaluation on
`J(7,4)` gives:

```text
b=4, phase 0:  one slack, A(X)=(|X cap B|-2)^2;
b=4, phase 1:  exactly 36 integer degree-two slacks.
```

The phase-one statement is a complete enumeration: the parity vector has
total lift mass eight, all 14 primitive left-kernel equations of the
rank-21 Johnson evaluation map are imposed, and all 36 solutions are
reconstructed by exact integral coefficient vectors.

## 2. Fixed-boundary coefficient exhaustion

Square multiplication in `F_49`, followed optionally by Frobenius, gives a
48-element stabilizer of infinity, finite zero, the distinguished edge,
Paley signs, and direction type.  The 58,800 saturated boundaries form
1,225 orbits.

For each orbit representative the exact finite model selects the 29 edges
of `H`, fixes the distinguished edge, imposes the four requested odd-degree
vertices and `c_H=-1`, and enforces all eight minimum-mean directional score
polynomials.  On `sum z_s=1`, each score polynomial is imposed by 22 sparse
coefficient equations, for 176 equations total.  A phase-one `b=4`
direction selects one of the complete 36-element target catalog by an exact
allowed-assignment constraint.

All 1,225 fixed-boundary models are `INFEASIBLE`:

```text
INFEASIBLE 1,225, UNKNOWN 0, FEASIBLE 0.
```

Their orbit sizes sum to all 58,800 saturated boundaries.

## 3. Exact transfer to the other product sign

Choose a nonsquare `alpha` in `F_49`, let

```text
P(infinity)=infinity,    P(u)=alpha*u,
s_infinity=-1,           s_u=+1 for finite u.
```

The Paley conference matrix satisfies the entrywise identity

```text
s_i C[P(i),P(j)] s_j = -C[i,j].                            (1)
```

The map fixes infinity, finite zero, and hence the distinguished edge.  Its
signed permutation anti-commutes with `C`, so it exchanges the two affine
eigenshells.  For the correspondingly transformed Boolean vector and edge
set, the raw affine score changes sign while the eigensign also changes
sign.  Thus the normalized product `epsilon*S` is preserved.

For any edge set `H`, (1) gives

```text
c(PH)=(-1)^|H| * product_v s_v^deg_H(v) * c(H).
```

Here `|H|=29` is odd.  Since the boundary contains four finite points and
not infinity, `deg_H(infinity)` is even.  Therefore `c(PH)=-c(H)`.  The map
is a bijection between the saturated survivors for the two signs, so the
single-sign 1,225-orbit exhaustion closes both signs.

## 4. Independent audit and archive

The independent audit re-enumerates all 211,876 boundaries for each sign,
reconstructs every square-semilinear orbit, checks the complete phase-zero
and phase-one slack catalogs, validates every result row against its orbit
representative and model scope, and verifies (1), exact affine-shell
exchange, boundary-set transfer, and product-sign transfer.  It reports:

```text
proved true
saturated boundaries per sign 58,800
saturated orbits per sign 1,225
certificate rows 1,225
missing 0, duplicate 0, malformed 0, unknown 0, feasible 0
remaining unsaturated boundaries per sign 23,520
remaining unsaturated orbits per sign 518.
```

Permanent archive:

```text
/mnt/storage/e1work/quadratic-minmax-limit-finite/
  2026-08-26-p7-four-point/
  p7_no_infinity_saturated_certificate_2026-08-26.tar.gz
SHA256 5234f60d246c50dcdb7f5feb51b23185e8eecbffac431f04c279c1cc02153612
```

Independent audit SHA-256:

```text
4d8fbbba46a0f7b19fc2f2241e4a710bf81fc581b5acb893f834ee65c1f547ea
```

Orbit-source SHA-256:

```text
7f7d3cc26077bb40ac096b638c6fc20ddf1a8fe6ddee60641f2fb568bacfd077
```

Core scripts:

- `scripts/p7_size_four_slack_classify.py`;
- `scripts/p7_no_infinity_saturated_cpsat.py`;
- `scripts/p7_no_infinity_saturated_orbit_batch.py`;
- `scripts/p7_no_infinity_saturated_audit.py`.
