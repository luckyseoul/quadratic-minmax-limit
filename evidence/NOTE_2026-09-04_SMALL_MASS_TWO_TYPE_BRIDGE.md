# Proposition 15.774: small mass and two-type capacity

Date: 2026-09-04.

Classification: **proved infinite-band theorem**, using the local
small-mass spectrum stated below. It advances residual (ii) and the
separate minimal-four-gap implication bridge, but does not prove global
residual closure, E1, or the limit. The protected 15.766/15.767 records
are neither edited nor treated as new inputs requiring a rerun.

## 1. Local input and its phase-independent consequence

The core local lemma is proved in
[the sharp small-mass note](NOTE_2026-09-04_SHARP_SMALL_MASS_SPECTRUM.md):
for every prime `p>=29`, a nonzero nonnegative integer-valued quadratic `C`
on `J(p,(p+1)/2)` with

```text
0 < M_C=4p E[C] < L_p:=2p-10
```

is Boolean, and its mass belongs to `{p-3,p+1}`. The strict upper
endpoint matters: mass `2p-10` is not excluded or classified here.

Consequently, if `A` is a nonnegative integral quadratic whose parity
is an affine parity with an even active-coordinate set on this slice, then

```text
a=2p E[A] < L_p  =>  a in {0,p-3,p-1,p+1}.              (1)
```

This is a union over the possible parity phases and boundaries, not an
assertion that every listed mass occurs at each fixed boundary. Constant
odd parity costs at least `2p`. Constant even parity gives `A=2C` and
the core lemma. For nonconstant parity, the established phase floors
below `L_p` leave only the pointwise two-bit or omitted-bit baselines,
of mass `p-1` or `p+1`. A positive half-difference from either baseline
has mass below `p-3`, contradicting 15.688. The nonminimal
complement-triple floor `2p-6` is above `L_p`, so no punctured
subtraction is used in this deduction.

## 2. Universal signed-shell identities

Let `r in {3,4,5}` and suppose a finite edge set `H` satisfies

```text
|H|=rp+2t,       t>=0,
T_H^epsilon(y)>=r  for every y in E_epsilon, epsilon=+1,-1.
                                                               (2)
```

Here `T_H^epsilon=epsilon*sum_(uv in H) C_uv y_u y_v` for the Paley
conference matrix. The signed frame mean excludes `|H|<rp`, so (2)
covers every size needing consideration at baseline `r`.

All ranges below have an isolated vertex. Transport it to infinity and
recompute the finite signed total `T=sum_(uv in H) C_uv`; then `I=0`.
Put `m=(p+1)/2` and `M=p+1=2m`. For a direction of quadratic type
`epsilon`, let `P_d` be its parallel H-edge count and define

```text
A_d=(T_H^epsilon-r)/2,       a_d=2p E[A_d].
```

The functions `A_d` are nonnegative integral quadratics. Their parity is
affine: the product of the edge signs equals a constant times the
product of the vertex bits on the odd-degree boundary, hence a parity
of the affine fibre bits. Shifting the baseline changes only its phase.
The active-coordinate set is even because the transported odd-degree
boundary is even. In particular the argument below uses the union over
phases; at even baseline four it does not assume the two types have
opposite phases or that either type has phase one.

A parallel edge has signed expectation one; a nonparallel edge has
signed expectation `-epsilon*C_uv/p`. Therefore, without an equality
classification,

```text
a_d=M P_d-epsilon T-rp,
sum_(d:type epsilon) P_d=(|H|+epsilon T)/2,
sum_(d:type epsilon) a_d=M t.                           (3)
```

Every `a_d` is a nonnegative even integer. Each of the `m` directions
of one type has the same residue, so write uniquely

```text
a_d=2u_epsilon+M k_d,
0<=u_epsilon<m,       k_d in Z_(>=0),
sum_(d:type epsilon) k_d=t-u_epsilon.                  (4)
```

Adding the two congruences from (3) gives the coupling

```text
u_+ + u_- = r (mod m).                                 (5)
```

## 3. The uniform averaged-mass band

Suppose `p>=29` and `0<=t<=p-6`. Each type has average `2t<=2p-12`,
strictly below `L_p`. At least one row in each type therefore satisfies
(1), forcing its common residue into

```text
u_epsilon in {0,m-2,m-1}.                              (6)
```

Their pair sums modulo `m` lie in `{0,m-4,m-3,m-2,m-1}`. Since
`m>=15`, none equals `r=3,4,5`. This contradicts (5).
Thus (2) is impossible through

```text
r=3: |H|<=5p-12;
r=4: |H|<=6p-12;
r=5: |H|<=7p-12.                                      (7)
```

Isolation is valid throughout: at the largest endpoint its margin is
`p^2+1-2(7p-12)=p^2-14p+25>0` for `p>=29`.

## 4. Sharper capacity, rather than only an average

The strict spectrum (1) supplies a lower bound on every quotient in
(4), depending only on its residue:

| Residue u | Necessary k_min(u) |
|---|---:|
| `0,m-2,m-1` | 0 |
| `1,...,m-7` | 2 |
| `m-6,...,m-3` | 1 |

For `1<=u<=m-7`, both `2u` and `2u+M` are positive, below `L_p`,
and absent from (1). For `m-6<=u<=m-3`, the zero-quotient mass is
forbidden but `2u+M>=L_p`; the endpoint `u=m-6` must be retained.
Consequently the whole type obeys

```text
t >= u + m*k_min(u).                                  (8)
```

Now let `p>=37` and fix `s_r=ceil(r/2)-1`, so `s_r=1,1,2` for
`r=3,4,5`. For every `t<=2m+s_r=p+ceil(r/2)`, (8) confines each
residue to

```text
{0,...,s_r} union {m-6,...,m-1}.                       (9)
```

We use the fixed upper bound here; an earlier actual `t` only removes
possibilities. In particular `u=0` is not lost when `t<2m`.
The three pair types all contradict (5):

- Low plus low is at most `2s_r<r`.
- Low plus high has residue at most `s_r-1`, or at least `m-6>r`.
- High plus high has residue at least `m-12>r` and at most `m-2`.

The last strict inequality uses `m>=19`. Thus the sharper theorem is

```text
r=3: |H|<=5p+4;
r=4: |H|<=6p+4;
r=5: |H|<=7p+6,           for every prime p>=37.        (10)
```

At the largest endpoint, isolation has margin
`p^2+1-2(7p+6)=p^2-14p-11>0` for `p>=37`.
These are necessary integer-capacity contradictions, not classifications
of the high-mass cells or assertions that every allowed residue is realizable.

The scalar threshold is sharp for this relaxation. At the first uncovered
layer `t=p+1+ceil(r/2)`, take `u_+=floor(r/2)` and
`u_-=ceil(r/2)`. Both types have all quotients two, except for one
quotient-three row in the lower-residue type when `r` is odd. Set each
actual parallel count to `P_d=r+k_d` and the signed totals to
`epsilon*T=r-2u_epsilon`. These satisfy (3)--(5), the exact total edge
budget, and the strict-spectrum capacity bounds; low rows have `P=r+2`.
All their masses are at least `L_p`. They are exact scalar survivors,
not local quadratic cells or graph realizations. Thus another use of only
these same mass residues and quotas cannot pass this threshold.

## 5. Complete the residual endpoint at p=29,31

The two missing primes in the `r=3` bound (10) are supplied by the
already audited continuation of 15.773. The following argument is valid
for every prime `p>=29`, but only p29/p31 are needed to complete (10).
Use the original residual layer notation

```text
t_original=q+s,       |H|=5p+2s,       s=1,2,
q=(p-1)/2,            sum_hard k_L=p+s-u.
```

Choose the phase-one hard sign using the transported edge product:
`h/c_H=(-1)^(q+s)`, with `c_H` recomputed after transport. The
universal identity remains
`hT=M P_L-3p-a_L`, with `a_L=2u+M k_L`.

For `u<q`, all quotients are positive. If a quotient-one row occurs,
its local classification is unchanged from 15.773. Apply the common-row
identity before an offset: `P_L=P+k_L-1`. The low count has `P<=9<q`,
so its coefficient congruence fixes its offset. The opposite total rises
by `s` and the old forced counts decrease by `s`:

```text
p1 branches: [4,8,7,7,6,6,5] - s;
p3 branches: [7,7,7,7,5,5] - s.                       (11)
```

The last entries record the quotient-zero endpoints discussed below.
Every count stays positive. Existing local exclusions at `p+11,p+13,p+15`
remain available for the classified carries; no high-row equality
classification or unique-high-row assumption is made.
The bound `P<=9` is on the chosen low count only; elevated `P_L` can
exceed nine and are not excluded by that bound. Nor is a quotient-one
row inferred from a possibly nonpositive counting lower bound: this
branch explicitly assumes one occurs.

If `u<q` and no quotient is one, all are at least two. Necessarily
`u<s`. Put `e=s-u-1`; then at least `m-e` rows have quotient two.
Their common actual parallel count `P` requires no equality classification,
and direct row accounting gives

```text
hard edges=mP+e,       0<=P<=9,
hT=M P-5p-2u-2,
a_opp(Q)=M(P+Q)-8p-2u-2.                              (12)
```

Writing `R=P+Q`, the value `R=8` has forbidden mass `6-2u`, while
`R<=7` is negative. Hence `R>=9`. The next mass is `p+7-2u`, attained
in at least `4-s-u` opposite rows. These are exactly

| (s,u) | Forced mass | Number of rows at least |
|---|---|---:|
| `(1,0)` | `p+7` | 3 |
| `(2,0)` | `p+7` | 2 |
| `(2,1)` | `p+5` | 1 |

At `u=q`, a zero quotient has the old exact baseline, with `P=4` in
either class or `P=3` in the p3 class. It forces `5-s` mass-`p+9`
opposite rows. With no zero quotient, at least `m-s` rows have quotient
one and mean `2p`. Their actual common low `P<=9` gives

```text
hard edges=mP+s,       hT=M P-5p,
a_opp(Q)=M(P+Q)-8p.
```

The same `R>=9` argument forces `5-s` mass-`p+9` rows. In both uses
of actual `P`, the case `P=9` does not invent a `Q=-1` row: formulate
the lower bound as `P+Q>=9` throughout.

All newly required masses `p+5,p+7,p+9` are below `L_p` for p29/p31
and absent from the local C-spectrum. For a forced phase-zero cell,
the only nonzero-boundary possibilities have pointwise XOR or omitted-bit
parity minima; their half-differences have positive mass below `p-3`.
The zero-boundary case is `A=2C` and is excluded by the local spectrum.
The same pointwise reduction applies to the classified carries.
Thus `s=1,2` are impossible at p29/p31. Together with (10) and the
previously closed layers, residual (ii) is empty for every prime `p>=29`
through separator size `|H|=5p+4`, equivalently even deletion size
`k<=5p+3` in its official range.

## 6. Exact minimal-four-gap and limit consequences

For an inclusion-minimal four-gap H, all deletions have gap two.
Proposition 15.764 gives the exact parity alternatives. An even H that
hits signed level two enters the already closed Type-I unit; otherwise
both signed shells have floor four. An odd H with no signed level-three
row has floor five in both shells; a level-three row supplies all
official residual-(ii) hypotheses, including the frozen distinguished
edge and the phase normalization.

The present theorems therefore imply:

- Odd minimal four-gap `H` is impossible through `5p+4` for `p>=29`.
  In this range the no-bridge alternative is excluded by (7)/(10),
  while the critical-row alternative is excluded by the residual result.
- Even minimal four-gap `H` is impossible through `6p-12` for all
  `p>=29`, and through `6p+4` for `p>=37`.
- The odd no-bridge branch is impossible through `7p-12` for all
  `p>=29`, and through `7p+6` for `p>=37`. This last statement alone
  does not exclude an odd witness that does have a critical row beyond
  the proved residual range.

The global acceptance functions must remain false. A minimal witness
can in principle have more edges than these bands, and no theorem here
bounds every minimal witness by a covered size. That missing all-size
quantifier is the remaining obstacle to deriving E1 from this work.

For the limit, finite small-prime bases are not intrinsically necessary:
eventual E1 for every sufficiently large prime would already give
`alpha_(p^2+1)->1/2`; the proved prime-order denseness lemma would then
give `L=1/2`. But the present bounded-support exclusions are not eventual
E1 over all signings. They give neither a global witness-localization
theorem nor an `o(p^3)` bound on the possible improvement over Paley.
No global residual, E1, or limit closure is claimed.
