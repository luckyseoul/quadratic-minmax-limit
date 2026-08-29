# Proposition 15.700: p=17 second-boundary slack-zero reduction

At `p=17`, the second all-finite even boundary has size `s=16`. This note
records an exact arithmetic and finite-geometric reduction. It does **not**
claim closure of the endpoint.

## Exact arithmetic census

For one quadratic direction type, write the common mean residue as `2u`
modulo `18`. Completion-bounded enumeration uses every even odd-fibre count
`0<=b<=16`, the exact directional floor, and the exclusion of a two-unit
nonzero lift.

The phase-zero relaxed residues are `u=0,2,3,4,5,6,7,8`; phase one has only
`u=0,8`. Proposition 15.688's sharp nonzero integral-quadratic lift floor is
`p-3=14`. For phase-zero `u=2,...,6`, a quotient-zero direction is forced,
its scaled mean is below the least positive-`b` floor `16`, and hence `b=0`.
The resulting nonzero lift has cost `2u<14`, a contradiction. Thus only

```text
phase zero: u = 0,7,8
phase one:  u = 0,8
```

remain. The exact pair budget and four-divisibility of pair slack then give
1,575 phase-labelled profiles. Exactly 247 have pair slack zero, split by
residue pair as

```text
(u0,u1)=(0,8): 234
(u0,u1)=(7,8):   4
(u0,u1)=(8,8):   9
```

## Exhaustive image of the unique 16-arc class

Pair slack zero means the sixteen affine points are a projective 16-arc.
Sticker's exhaustive classification records one PGL class of 16-arcs in
`PG(2,17)`; it is a nondegenerate conic with two points deleted.

Fix the conic `XZ=Y^2`. We enumerate all 307 projective lines as possible
lines at infinity and every deleted pair containing the conic's intersection
with that line. This gives exactly:

```text
external line: 20,808 cases
tangent line:     306 cases
secant line:      153 cases
total:         21,267 cases
```

For every case, all eighteen affine directional odd-fibre counts and their
Paley phases are recomputed. Including the global phase swap induced by a
nonsquare affine determinant gives 53 phase-labelled geometric profiles.
Their intersection with the 247 arithmetic profiles has size two:

```text
phase 0 {0:1, 2:7, 16:1}; phase 1 {2:9}
phase 0 {0:1, 2:8};       phase 1 {2:8, 16:1}
```

Both are tangent-at-infinity conic-minus-two cases and both have residue
pair `(u0,u1)=(0,8)`. Therefore 245 slack-zero profiles are impossible and
the full exact p=17 profile remainder falls from 1,575 to 1,330.

## Honest status and follow-up model

The two surviving conic profiles and every positive-slack profile remain
open. `scripts/p17_slack0_conic_edge_means_cpsat.py` fixes canonical
representatives and imposes the 69-edge count, exact odd-degree boundary,
Paley product sign, exact directional means, and all coefficient identities
forced by floor equality. Five-minute runs on Soulkiller and Nuka both
returned `UNKNOWN`; they are diagnostics only and are not used in
Proposition 15.700.

The literature and OEIS searches after this finding found no prior listing
of these 247-to-2 directional profiles. The external geometric input remains
the already-linked source:

- H. Sticker, *Classification of Arcs in Small Desarguesian Projective
  Planes*, Ghent University PhD thesis (2012), Section 5.3, printed page 119.

Machine-readable evidence is in `evidence/e1_gmin_m4_prop15700.json`; the
reproducible theorem code is `src/e1_gmin_m4_prop15700.py`.
