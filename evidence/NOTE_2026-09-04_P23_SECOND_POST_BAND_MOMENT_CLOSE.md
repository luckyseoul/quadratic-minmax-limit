# The p=23 second post-band endpoint closes with eleven hard roots

Date: 2026-09-04

Scope: the exceptional `p=23,t=10,k=112` layer covered by Proposition
15.770 for larger `p=3 mod 4` primes.  This is a proved all-boundary endpoint
theorem.  It reuses the fixed 33,649-five-set certificate from the preceding
`p=23,t=9` close; it is not a graph census, a later-layer theorem, or a global
closure of residual (ii).

## Exact residue ledger

Put

```text
p=23, q=11, m=12, t=10, k=112, |H|=113.
```

There are `23^2+1-2*113=304` guaranteed isolated vertices, so signed
transport again gives an all-finite chart with `I=0`.  Write

```text
a_L=2u+24 k_L,             sum_L k_L=22-u.
```

The phase-one floors and sharp lift floor 20 leave exactly three residues.

- `u=9`: eleven quotient-one rows of mean 42 and one quotient-two row; the
  low rows are the sharp mass-20 families classified at `t=9`.
- `u=10`: all rows have quotient one and mean 44; relative to either exact
  parity baseline the lift has scaled mass 22.
- `u=11`: the two old exact parity endpoints.

All `u<9` rows are positive lifts below the sharp mass-20 floor.

At `u=11`, the XNOR and complementary-literal ledgers respectively have

```text
(P,hT,hard,opposite)=(4,5,59,54),
(P,hT,hard,opposite)=(3,-19,47,66).
```

They exclude a mass-eight row and force a row of mass `p+9=32`, already
impossible by Proposition 15.752.

## The new mass-22 residue is absent

For the `u=10` lift, the height-at-least-two scaled floor is `p+1=24`, so a
mass-22 lift must be Boolean.  Its density is `11/46`.  The corrected
Johnson influence bound is

```text
L <= 5929/1058 < 6.
```

Thus at most five slice coordinates remain, all patterns extend because
`5<q=11`, and cube influence leaves at most four active coordinates.  The
fixed four-bit density list is

```text
0, 5/23, 6/23, 11/23, 12/23, 17/23, 18/23, 1,
```

which omits `11/46`.  Hence the `u=10` residue is empty without a new finite
slice census.

## The carried sharp branch

Fix one of the four sharp hard families and let its coefficient offset be
`P in {2,3,4,5}`.  The eleven low rows have parallel count `P`.  The common
row identity forces the unique high row to have count `P+1`; therefore

```text
hard edges = 12P+1,
hT = 24P-111,
opposite edges = 112-12P.
```

For an opposite direction of count `Q`,

```text
a(Q)=24Q+hT-69.
```

The row `Q=8-P` has mass 12, below both the phase-zero nonzero-boundary floor
24 and the sharp lift floor 20.  Thus every opposite row has `Q>=9-P`.
The surplus after assigning this value to all twelve directions is exactly
four, so at least eight directions have `Q=9-P` and mass 36.  Every
nonzero-boundary alternative would leave a positive lift of mass 12 below
the sharp floor, hence each such row has boundary zero and is `A=2C` with
`92 E[C]=36`.

The preceding endpoint theorem classifies and globalizes every such `C` as
`F4` or `F5`, with coefficient offsets one and five.  The opposite
congruence `11 | Q-offset` leaves only

```text
P=4, Q=5, C=F5.
```

Consequently every one of the eleven low hard rows is the same
complement-literal plus all-equal-triple family.  Their common sign is `h`.
Because `I=0` and `P=4` equals the hard target offset, the slice-kernel scalar
vanishes and their actual normalized coefficient graphs are
triangle-minus-full-star.  Likewise `Q=5` equals the opposite `F5` offset,
so every forced opposite graph is an actual `K5`.

## Eleven roots still force the two moment forms to vanish

For the global even moments `M_d(L)`, every low hard direction is a root of

```text
G4 = 2h M4-M2^2,
G8 = 24h M8-32M2 M6+5M2^4.
```

The eleven low directions are distinct projective directions.  Since
`11>8`, they force both homogeneous binary forms, of degrees four and eight,
to vanish identically.  The unique high hard row is not used.

On an opposite `F5/K5` row the two identities become

```text
-2S4-S2^2=0,
-24S8-32S2*S6+5S2^4=0.
```

The authoritative exact certificate over all `binom(23,5)=33,649` five-sets
has 1,518 quartic zeros, 2,024 octic zeros, and no simultaneous zero.  Since
at least eight such opposite rows are required, the branch is impossible.

Therefore

```text
p=23,t=10,k=112 is empty for every boundary size.
```

The next `p=23` layer `t=11`, later layers, residual (ii) globally, E1,
`L=1/2`, and the original MathOverflow limit remain open.

## Replay

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_p23_second_post_band_moment_close.py
```

Machine-readable evidence:
`evidence/e1_gmin_m4_p23_second_post_band_moment_close.json`.
