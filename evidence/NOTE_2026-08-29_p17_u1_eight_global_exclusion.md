# Proposition 15.709: every remaining p17 profile with `u_1=8` is impossible

After Proposition 15.708, the exact p17 ledger contains 507 profiles. They
split by residue pair as

```text
(u_0,u_1)=(0,8):  66
(u_0,u_1)=(8,8): 214
(u_0,u_1)=(0,0): 181
(u_0,u_1)=(7,0):   9
(u_0,u_1)=(8,0):  37.
```

The two solver-free contradictions in Proposition 15.708 exclude the first
two blocks uniformly, independent of pair slack.

## The 66 `(0,8)` rows

Every row retains at least eight rigid phase-one `b=2` directions. It also
retains a rigid phase-zero `b=0` direction; the exact lower-bound histogram
is

```text
retained b=0:  3   4   5
profiles:     10  27  29.
```

The global finite-edge Paley-sign comparison gives

```text
17I = 4 + 72(g_0+g_1).
```

Thus `I=68` and `g_0+g_1=16`. Nonnegative parallel counts require
`g_0>=9` for the `b=0` direction and `g_1>=8` for the `b=2` direction, an
impossible sum. All 66 rows are excluded.

## The 214 `(8,8)` rows

Again every row retains at least eight rigid phase-one `b=2` directions. It
also retains rigid phase-zero `b=16` directions, with lower-bound histogram

```text
retained b=16:  2   3   4   5   6   7   8
profiles:       4  30  36  36  36  36  36.
```

Comparing those two floors forces `I=4`, `g_16=1`, and `g_2=0`. Proposition
15.708's unique-even-fibre identity applies without reference to the pair
slack:

```text
N_j-delta_j = -15z_j-I,
```

where `N_j>=0`, `delta_j<=1`, and `z_j>=0`. Hence `N_j<=-3`, impossible.
All 214 rows are excluded.

The 280 excluded profiles occupy every slack block from 28 through 128. The
ledger falls from 507 to 227 profiles. Every survivor has `u_1=0`; their
pair slack is at least 96. The p17 endpoint and all top-level gates remain
open. No solver or new classification is used.

Targeted GitHub code and MathOverflow searches found no prior version of this
full-ledger rigid-anchor exclusion. OEIS count searches produced only
unrelated isolated matches and are not used as evidence.
