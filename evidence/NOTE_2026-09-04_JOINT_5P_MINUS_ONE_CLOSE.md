# Proposition 15.773: the joint generic layer k=5p-1 is empty

Date: 2026-09-04.

Classification: **proved infinite-family theorem**. For every prime
`p>=29`, residual (ii) is empty at `t=q=(p-1)/2`, equivalently
`k=5p-1`, for every boundary size. This advances both generic congruence
classes by one layer. It does not close any new small-prime endpoint,
residual (ii) globally, E1, or the original limit problem.

The changed premise is the completed local and common-row theory through
15.772. The genuinely new branch has all hard means `2p`, but its common
row identity closes it without classifying those local equality cells.
No new equality catalog, graph search, prime census, or mesh search is
part of the proof.

## 1. Isolated chart and exhaustive quotient split

Put `q=(p-1)/2`, `m=q+1`. A witness at this layer has `|H|=5p`, so
at least `p^2-10p+1>0` vertices are isolated. Signed PSL transport puts
one at infinity with `I=0`. Every directional boundary size is even.
Choose the hard sign `h` to have phase one. All signs below refer to the
transported chart: recompute `c_H` there, rather than retaining its value
before transport. For even boundary and `I=0`, the exact phase rule is
`(-1)^eta=epsilon*(-1)^((5p-3)/2)*c_H`. Hence the hard sign is
`h=c_H` for `p=1 mod4`, and `h=-c_H` for `p=3 mod4`.

Here is the row identity without any equality or offset normalization.
Define `T=sum_(e in H) C_e` in this chart. Each edge parallel to a hard
direction has `h*C_e=1`. Consequently the constant term in `h*S_H` is
the actual parallel count `P_L`, and the sum of its off-fibre quadratic
coefficients is `hT-P_L`. On the middle slice, with `z_i=2x_i-1`,
`E[z_i*z_j]=-1/p` for distinct fibres. Averaging `h*S_H=3+2A_L`
therefore gives `2p E[A_L]=(p+1)P_L-hT-3p`. Replacing `h` by `-h`
gives the opposite identity. This derivation works for both signs and
every row, not merely a classified low row.

The exact quotient equations, with all `k_L` nonnegative integers, are

```text
a_L=2u+(p+1)k_L,       0<=u<=q,
sum_hard k_L=p-u=2q+1-u,
hT=(p+1)P_L-3p-a_L.                                  (1)
```

Every phase-one floor is at least `p-1`. If `u<q`, quotient zero is
therefore impossible, and at least

```text
2m-(p-u)=u+1
```

hard quotients equal one. This is one fewer low row than at `t=q-1`,
but is always positive. At `u=q`, the quotient sum equals `m`.
Either at least one quotient is zero, or every quotient is exactly one.
These are exhaustive alternatives; the latter is the new flat branch.

For `p=1 mod4`, the low-row classifications already proved in 15.772 give:

| Residue | Low quotient | Low mean | Possible low families / offsets |
|---|---:|---|---|
| `0` | 1 | `p+1` | complement literal / 5 |
| `1,...,q-4` | 1 | `p+1+2u` | none |
| `q-3` | 1 | `2p-6` | complement triple / 2 |
| `q-2` | 1 | `2p-4` | XNOR plus sharp lift / 3,5 |
| `q-1` | 1 | `2p-2` | gap-four or literal-plus-sharp families / 4,6 |
| `q`, with a zero | 0 | `p-1` | exact XNOR / 4 |

In particular the punctured complement-triple gap-two theorem, its
gap-four equality classification, and the p1 mass-`p-1` exclusion retain
their precise roles. Do not assert that the complement-triple difference
is nonnegative on its omitted layer.

For `p=3 mod4`, the low mean is below `2p` whenever `u<q`, so the
only possible boundaries are `b=2,p-1`, both with pointwise baseline
mean `p-1`. The genuine half-difference has mass `2u+2`. Thus:

| Residue | Low quotient | Low mean | Possible low families / offsets |
|---|---:|---|---|
| `0,...,q-3` | 1 | `p+1+2u` | none: positive lift mass below `p-3` |
| `q-2` | 1 | `2p-4` | either baseline plus sharp lift / 2,3,4,5 |
| `q-1` | 1 | `2p-2` | none: mass-`p-1` exclusion in 15.770 |
| `q`, with a zero | 0 | `p-1` | complement literal / 3, or XNOR / 4 |

The four sharp offsets are respectively obtained from baseline offsets
3,4 and lift increments -1,+1. Their classification is the existing
15.769 theorem, with arbitrary support overlaps allowed.

## 2. Common-row normalization before selecting an offset

Let a classified low row have quotient `ell`, mean
`a=2u+(p+1)ell`, and parallel count `P`. Equation (1), applied before
any offset substitution, gives

```text
P_L=P+k_L-ell,
hard edges=m(P-ell)+(p-u),
hT=(p+1)P-3p-a.                                       (2)
```

For `ell=1,u<q`, nonnegative opposite count implies

```text
P<=floor((5p-(p-u)+m)/m)<=9<q.
```

For `ell=0,u=q` it gives `P<=floor((5p-m)/m)<=8<q`.
The isolated coefficient congruence now forces `P` to equal its
classified offset. Thus different offsets cannot mix among low rows,
while families of the same offset may mix. No high-row equality
classification is used.

For the p3 sharp residue, `p-u=m+2`; hence at least `m-2` rows are low,
and the total carry above quotient one is two. It may occur as one row
of quotient three or two rows of quotient two. In either case (2) gives
exactly `hard edges=mP+2`. A unique-high-row assumption would be wrong.

## 3. Complete carried-row ledger

Write `Q` for an opposite parallel count. Its mean and the opposite
total are

```text
a_opp(Q)=(p+1)Q+hT-3p,
sum Q=5p-hard edges.                                  (3)
```

Equations (2)--(3) give the following full ledger. An entry with several
values of `P` denotes one row for each value.

| Class and low family | P | Hard edges | hT | Forbidden Q / mass | Next Q / mass | Next rows at least |
|---|---:|---|---|---|---|---:|
| p1 complement literal | 5 | `3p+2` | `p+4` | `2 / 6` | `3 / (p+7)` | 4 |
| p1 complement triple | 2 | `p+4` | `8-3p` | `6 / 14` | `7 / (p+15)` | 8 |
| p1 sharp XNOR | 3,5 | `mP+2` | `(p+1)P-5p+4` | `(8-P) / 12` | `(9-P) / (p+13)` | 7 |
| p1 gap-four / literal-plus-sharp | 4,6 | `mP+1` | `(p+1)P-5p+2` | `(8-P) / 10` | `(9-P) / (p+11)` | 6 |
| p3 sharp baseline-plus-lift | 2,3,4,5 | `mP+2` | `(p+1)P-5p+4` | `(8-P) / 12` | `(9-P) / (p+13)` | 7 |
| either class, zero-quotient XNOR | 4 | `5m` | 5 | `3 / 8` | `4 / (p+9)` | 5 |
| p3 zero-quotient literal | 3 | `4m` | `4-p` | `4 / 8` | `5 / (p+9)` | 5 |

In each row, smaller `Q` gives negative mean. The positive forbidden
mass is at most 14, below both the nonzero phase-zero boundary floor
and the sharp zero-boundary lift floor `p-3`. Hence every opposite
count reaches the displayed next value. Subtracting that baseline from
the total in (3) leaves surplus `m-d`, where `d` is the last column.
At least `d` directions therefore attain the next value exactly.

These are precisely the prior row contradictions with opposite surplus
increased by one. All displayed lower counts remain positive, including
at the first primes 29 and 31.

## 4. The flat mean-2p branch needs no local classification

At `u=q`, suppose every quotient is one. All hard means are `2p`.
Equation (1) makes their actual parallel count `P` common immediately:

```text
hard edges=mP,           0<=P<=floor(5p/m)=9,
hT=(p+1)P-5p,
a_opp(Q)=(p+1)(P+Q)-8p.                               (4)
```

No coefficient offset or equality classification is needed to obtain
this integer `P`. If `P+Q<=7`, the mean is negative. If `P+Q=8`,
it is the forbidden positive mass eight. Therefore `P+Q>=9`, or
`Q>=9-P`. This argument also covers `P=9`: the formal mass-eight
index is then `Q=-1`, not an actual opposite row, and the inequality
already follows from `Q>=0`.

The surplus over the common minimum is independent of `P`:

```text
sum(Q-(9-P))=5p-mP-m(9-P)=(p-9)/2=m-5.               (5)
```

At least five opposite rows have `Q=9-P`, and every such row has
mean `p+9`. This is exactly the old local forbidden mass, not a new
mean-`2p` equality problem.

## 5. Exclude every forced opposite row and conclude

At any forced mass `p+c`, where `c` is one of `7,9,11,13,15`, the
phase-zero floors leave only `b=0,2,p-1` in its indicated prime class.
At `b=2`, subtract the pointwise XOR parity minimum. At `b=p-1`,
subtract the pointwise omitted-coordinate bit of the required parity
(`1-x_j` in the p1 phase-zero case, `x_j` in the p3 case).
After division by two the remaining quadratic is nonnegative and
integral, with positive scaled mass below `p-3`, impossible by 15.688.

At `b=0`, write `A=2F`; then `4p E[F]=p+c`. The required local
exclusions are already proved:

- `p+7`: Proposition 15.751;
- `p+9`: Proposition 15.752, valid for both classes here;
- `p+11` in the p1 class: the local theorem used by 15.772;
- `p+13`: 15.770 in the p1 class and 15.769 in the p3 class;
- `p+15` in the p1 class: Proposition 15.768.

Thus every carried branch and the new flat branch is impossible.
All residues and boundary sizes are exhausted, proving the claim for
every prime `p>=29`.

The common generic frontier is now `t>=q+1`, equivalently `k>=5p+1`.
The next layer is not attacked here. The p23 frontier remains `t>=12`;
the critical p5/p7 cases, p11 and p13/p17/p19 later layers, the positive
`p=7,z=7` branch, and all global acceptance gates remain open.

Dependencies: `src/e1_gmin_m4_prop15769.py`,
`src/e1_gmin_m4_prop15770.py`, `src/e1_gmin_m4_prop15772.py`, and
`evidence/NOTE_2026-09-04_COMPLEMENT_TRIPLE_PUNCTURED_GAP.md`.
The executable theorem is `src/e1_gmin_m4_prop15773.py`.
Independent mesh arithmetic and regression checks corroborate the proof;
they do not replace it or enlarge its range.
