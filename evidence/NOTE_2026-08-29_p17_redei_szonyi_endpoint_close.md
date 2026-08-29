# Proposition 15.712: the p17 second all-finite endpoint is closed

Proposition 15.711 leaves fourteen exact profiles. Every one has phase-one
directional profile `{16:9}`. Szőnyi's sharp finite-affine direction theorem
excludes all fourteen at once.

## Nine missing directions

Let `B` be the sixteen-point affine odd boundary. For direction `d`, let
`n_{d,t}` be its seventeen fibre occupancies and let `b_d` count the fibres
of odd occupancy. If `b_d=16`, then

```text
16=sum_t n_{d,t} >= #{t:n_{d,t} odd}=16.
```

Equality forces sixteen singleton fibres and one empty fibre. No two points
of `B` lie in a common `d`-fibre, so `d` is not a direction determined by
`B`. All nine phase-one directions have `b_d=16`; hence `B` determines at
most the nine phase-zero directions.

## Szőnyi's direction bound

Szőnyi proved that a noncollinear `k`-point subset of `AG(2,p)`, for prime
`p` and `k<=p`, determines at least `(k+3)/2` directions. Here

```text
p=17, k=16, ceil((k+3)/2)=10.
```

The upper bound nine contradicts the theorem, so `B` must be collinear.

A sixteen-point subset of an affine line has occupancy profile
`(16,0,...,0)` in the line direction and `(1,...,1,0)` in every other
direction. Thus its phase-labelled profile must be

```text
phase zero: {0:1,16:8}
phase one:  {16:9}.
```

The line direction is necessarily phase zero because every phase-one
direction already has `b=16`. This collinear profile is absent from the
fourteen-row exact ledger. Therefore no row survives.

The theorem used is T. Szőnyi, *On the number of directions determined by a
set of points in an affine Galois plane*, J. Combin. Theory Ser. A 74
(1996), 141--146. The exact `k<=p` statement is also recorded in G. Somlai,
[*A new proof of Rédei's theorem on the number of directions*](https://doi.org/10.1007/s00013-024-01979-x),
Arch. Math. 122 (2024), 575--580.

The exact p17 ledger falls from fourteen profiles to zero. This closes the
`p=17,s=16` second all-finite endpoint. It does not by itself close residual
(ii), R1, Type I, or the limit. No solver or new classification is used.
