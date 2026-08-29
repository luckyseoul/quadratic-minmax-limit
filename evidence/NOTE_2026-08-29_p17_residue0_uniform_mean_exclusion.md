# Proposition 15.711: the five p17 residue-zero profiles are impossible

Proposition 15.710 leaves nineteen profiles, including five with
`(u_0,u_1)=(0,0)`. Their phase-one profile is always `{16:9}`; their
phase-zero profiles are

```text
{0:9}, {0:7,2:2}, {0:5,2:4}, {0:3,2:6}, {0:1,2:8}.
```

All five are excluded analytically.

## Avoiding the previous anchor forces uniform mean 18

If any phase-zero `b=0` direction remains at quotient zero, its floor is
rigid and Proposition 15.710 already excludes the witness. In the five rows,
the phase-zero quotient sum is nine. The minimum quotient is zero for each
`b=0` direction and one for each `b=2` direction. Raising every `b=0`
direction therefore spends every free increment exactly once. Thus every
phase-zero direction has quotient one and scaled mean 18. The nine
phase-one `b=16` directions also have quotient one and mean 18.

## Four global candidates

Let `I` be the infinity degree, `P_+` and `P_-` the common parallel-edge
counts in the two phases, and `Sbar=c_H*S` the normalized finite Paley-sign
sum. The directional mean equations are

```text
18*P_+ = 69-I+Sbar,
18*P_- = 69-I-Sbar.
```

Hence `69-I=9(P_++P_-)`. Infinity is not in the sixteen-point odd boundary,
so `I` is even. Therefore

```text
I = 6+18k,  k=0,1,2,3,
P_++P_- = 7-2k.
```

For a rigid phase-one `b=16` direction,

```text
P_+ = 5+8g-I.
```

Modulo eight, `P_+` equals the upper bound `7-2k`. Thus

```text
P_+=7-2k,  P_-=0,  g=1+2k.
```

Every finite edge consequently has phase one.

## The fibre-capacity contradiction

Fix a phase-one direction and let `j` be the unique even fibre. Write `z_s`
for the number of infinity edges ending in fibre `s`. The rigid cell identity
is

```text
L_st = g-z_s-z_t+1_{j in {s,t}}.
```

Because every finite edge has phase one, every `L_st` is a nonnegative
crossing-edge count. If `m=max_{t!=j}z_t`, the ordinary-fibre cells and the
cells incident with `j` imply

```text
I <= g+1+15*min(m,g-m)
  <= g+1+15*floor(g/2).
```

The four candidates fail as follows:

```text
I   g   P_+   upper bound on I
6   1    7             2
24  3    5            19
42  5    3            36
60  7    1            53.
```

Thus all five residue-zero profiles are impossible. Fourteen p17 profiles
remain, with residue split `(7,0):9,(8,0):5` and slack histogram
`{96:2,100:3,104:3,108:2,112:2,116:1,128:1}`. The endpoint and every
top-level gate remain open. No solver or new classification is used.

## Context search

Targeted GitHub-code and MathOverflow searches found no prior occurrence of
the exact identities `69-I=9(P_++P_-)`,
`L_st=g-z_s-z_t+1_{j in {s,t}}`, or
`I<=g+1+15 floor(g/2)`, nor of the resulting Paley fibre-capacity argument.
Related literature concerns Paley spectra, character sums, and equitable
partitions, but not this exclusion. OEIS matches for the candidate lists
`6,24,42,60` and `2,19,36,53` are unrelated digit-sum and fractional-part
sequences. These searches provide context only, not mathematical evidence or
a formal priority claim.
