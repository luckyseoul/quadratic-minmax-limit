# Proposition 15.706: both p=17 slack-zero profiles are impossible

The two profiles left by Proposition 15.700 are excluded by a common
solver-free signed-edge identity. This removes the entire slack-zero block;
the `p=17` second all-finite ledger now contains 639 profiles, all of pair
slack at least 20.

## Rigid directions of both quadratic types

Every exact mean allocation in either profile retains a rigid `b=2`
direction of each quadratic type. For one such direction let

```text
I       infinity degree
P_d     number of finite edges parallel to d
eps     quadratic type of d
c       Paley product sign c_H
sigma   eps*c
g_d     integral coefficient gauge
```

Its floor mean is `M_d=17-sigma`, and its unique charged inter-fibre cell
has target sum `sigma`. Summing the 136 rigid cell identities gives

```text
R_d = sigma + 136*g_d - 16*I.
```

The exact directional mean is also

```text
M_d = I + 17*P_d - R_d - 51.
```

Consequently

```text
P_d = 4 + 8*g_d - I.                         (1)
```

## One global finite-edge sign sum

Let `S` be the sum of the Paley signs of all selected finite edges. Every
edge parallel to a direction of type `eps` has Paley sign `eps`. Therefore
the signed cross sum in that direction is

```text
R_d = eps*S - P_d.
```

Combining this with (1) yields

```text
S = c + eps*(4 + 144*g_d - 17*I).            (2)
```

Apply (2) to one rigid direction of each type and equate the two expressions
for the same global `S`:

```text
17*I = 4 + 72*(g_+ + g_-).
```

Since the inverse of 17 modulo 72 is 17, this forces

```text
I = 68 (mod 72).
```

The range `0<=I<=69` leaves only `I=68`. A 69-edge graph with infinity
degree 68 has one finite edge. Its 68 star endpoints have their parity
toggled at zero, one, or two endpoints by that finite edge, so its affine
odd boundary has size 70, 68, or 66. It cannot have the required size 16.

Thus both profiles are impossible for every infinity degree. No numerical
solver result or unproved classification is used beyond Proposition 15.700's
already-audited identification of the two profiles.

Machine-readable evidence is in `evidence/e1_gmin_m4_prop15706.json`; the
proof predicate is `src/e1_gmin_m4_prop15706.py`.
