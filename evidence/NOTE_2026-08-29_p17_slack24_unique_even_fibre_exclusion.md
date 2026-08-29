# Proposition 15.708: all 54 p17 slack-twenty-four profiles are impossible

The pair-slack-twenty-four block has the exact residue split

```text
(u_0,u_1)=(0,8): 45 profiles
(u_0,u_1)=(8,8):  9 profiles.
```

All 54 profiles are excluded analytically. No solver result or new finite-arc
classification is used.

## The 45 `(0,8)` profiles

Every profile retains at least eight rigid phase-one `b=2` directions.
Quotient accounting also retains rigid phase-zero `b=0` directions with
exact lower-bound histogram

```text
retained:  2  3  4   5  6
profiles:  2  6 18  17  2.
```

Comparing one rigid direction of each phase in Proposition 15.706's global
Paley-sign identity gives

```text
17*I = 4 + 72*(g_0+g_1).
```

Thus `I=68` and `g_0+g_1=16`. Parallel-count nonnegativity for the rigid
phase-zero `b=0` direction requires `g_0>=9`; for the rigid phase-one `b=2`
direction it requires `g_1>=8`. This contradicts their sum. Hence all 45
profiles are impossible.

## The nine `(8,8)` profiles

Every profile retains rigid phase-zero `b=16` directions. Their exact lower
bounds are

```text
3, 3, 3, 3, 2, 2, 2, 2, 2.
```

They also retain at least eight rigid phase-one `b=2` directions. The
canonical phase-zero `b=16` floor is `A_d=1-x_j`, where `j` is the unique
even fibre. Its pure-pair target is `-1` exactly on pairs incident with `j`,
and zero otherwise. Comparing the two rigid phases in the global sign
identity gives

```text
17*I = -4 + 72*(g_16+g_2).
```

Consequently

```text
I=4,  g_16=1,  P_16=7,  g_2=0,  P_2=0.
```

Let `z_s` be the number of selected infinity edges ending in fibre `s`, so
`sum_s z_s=I=4`. Let

```text
L_st = sum eps_d*C_e
```

over selected finite edges crossing fibres `s,t`. The exact rigid cell
identity is

```text
L_st = g_16-z_s-z_t-1_{j in {s,t}}.
```

Since `g_16=1`, for every `t!=j`,

```text
L_jt = -z_j-z_t.                              (1)
```

There are 65 finite selected edges. The rigid phase-one identity gives their
normalized Paley-sign sum as `-63`, so exactly 64 edges have phase zero and
one has phase one. For the fixed phase-zero direction `d`, each crossing
phase-zero edge contributes `+1` to `L`, while the unique phase-one edge
contributes `-1`.

Let `N_j>=0` count phase-zero crossing edges incident with fibre `j`, and let
`delta_j` be zero or one according as the unique phase-one edge is not or is
incident with that fibre. Summing (1) over the other sixteen fibres yields

```text
N_j-delta_j = -16*z_j-sum_{t!=j}z_t
            = -15*z_j-I.
```

Therefore

```text
N_j = delta_j-15*z_j-I <= 1-4 = -3,
```

contradicting `N_j>=0`. All nine profiles are impossible.

The p17 ledger falls from 561 to 507 profiles, all of pair slack at least 28.
The p17 endpoint and every top-level gate remain open.

Targeted GitHub code and MathOverflow searches found no prior occurrence of
the distinctive global congruence or unique-even-fibre identity. OEIS count
searches produced only unrelated isolated matches; no sequence result is used
as evidence.
