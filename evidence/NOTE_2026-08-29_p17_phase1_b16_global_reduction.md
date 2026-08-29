# Proposition 15.710: complementary global identities leave nineteen p17 profiles

Proposition 15.709 leaves 227 profiles, all with phase-one residue zero. In
every one, all nine phase-one directions have a rigid `b=16` floor. Two
global Paley-sign comparisons exclude 208 rows without a solver.

## 176 rows with a rigid phase-zero `b=0` direction

The exact lower-bound histogram for retained phase-zero `b=0` directions is

```text
retained:  1  2  3  4  5  6
profiles:  8 26 32 56 38 16.
```

For a rigid phase-zero `b=0` direction, the normalized global finite-sign
sum is

```text
-3-144*g_0+17*I.
```

For a rigid phase-one `b=16` direction it is

```text
21+144*g_16-17*I.
```

Equating them gives

```text
17*I = 12 + 72*(g_0+g_16).
```

Thus `I=60` and `g_0+g_16=14`. Nonnegative parallel counts require

```text
P_0  = 3+8*g_0-I  >= 0  => g_0>=8,
P_16 = 5+8*g_16-I >= 0  => g_16>=7.
```

Their minimum sum is 15, a contradiction. All 176 rows are impossible.

## Thirty-two rows with rigid `b=16` directions in both phases

These rows have phase-zero rigid-`b=16` lower-bound histogram
`{1:4,2:4,3:4,4:4,5:4,6:4,7:4,8:4}`. Comparing their phase-zero and phase-one
global identities gives

```text
17*I = 4 + 72*(g_16^-+g_16^+).
```

Thus `I=68` and the gauge sum is 16. Parallel nonnegativity requires the
phase-zero gauge to be at least 9 and the phase-one gauge at least 8, whose
sum is at least 17. These thirty-two rows are also impossible.

## Exact remainder

Nineteen profiles survive, all with phase-one profile `{16:9}`. Their residue
split is `(0,0):5`, `(7,0):9`, `(8,0):5`, and their slack histogram is

```text
96:3, 100:4, 104:4, 108:3, 112:3, 116:1, 128:1.
```

The five `(0,0)` rows are the profiles at slack 96, 100, 104, 108, and 112
with phase-zero signatures `{0:9}`, `{0:7,2:2}`, `{0:5,2:4}`, `{0:3,2:6}`,
and `{0:1,2:8}`. The p17 ledger falls from 227 to nineteen profiles. The
endpoint and every
top-level gate remain open. No solver or new classification is used.

Targeted GitHub code and MathOverflow searches found no prior version of the
complementary global-sign argument. OEIS count/slack searches returned only
unrelated arithmetic sequences and are not used as evidence.
