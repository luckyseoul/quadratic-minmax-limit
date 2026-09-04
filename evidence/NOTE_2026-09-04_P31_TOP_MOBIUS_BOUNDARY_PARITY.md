# Adaptive kernel-selector obstruction at the all-active top endpoint

Date: 2026-09-04

**Status:** let `p=4r+3` be a prime for which the balanced all-active branch-C
top ledger applies (in the present proof chain, `p>=31`), and put
`m=(p+1)/2=2r+2`.  An exact adaptive-selector argument excludes the `j=0`
top localized-Mobius completion with one half in each hard direction and an
opposite fixed edge.  The argument is independent of all nonzero half
centres, all choices of representatives in the kernel selectors, all
relative auxiliary scales, and whether the sole cancellation unit is
realized as `1:1` or `2:1`.  This closes that precise localized-Mobius
endpoint branch.  It does **not** prove that every residual-(ii) source must
be a sum of localized halves, and therefore does not close residual (ii)
globally.

## 1. The localized discrepancy

Work over `V=F_p^2`.  Let `L,M` be independent functionals and let
`e_1,e_2` be their dual vector basis.  For a nonzero centre `c`, the `p-1`
raw Mobius edges are indexed by `t in F_p minus {-1}` and have endpoints

\[
 u_t=c\left(e_1+{t\over t+1}e_2\right),\qquad
 v_t=ct(e_1+e_2).                                      \tag{1}
\]

Put

\[
 \tau_t=\eta\bigl(Q(u_t-v_t)\bigr)\in\{+1,-1\}.
\]

The normalized forced half selects the edge

\[
                         \{-\tau_tu_t,-\tau_tv_t\}.     \tag{2}
\]

Its degree discrepancy from the hard affine star is

\[
 d(L,M,c)=\partial E(L,M,c)+\mathbf 1_{\{L=c\}}
          \quad\hbox{in }\mathbb F_2^V.                \tag{3}
\]

The map `t -> t/(t+1)` bijects `F_p minus {-1}` with
`F_p minus {1}`.  Hence the points `u_t`, together with `v_1`, are exactly
the affine line `L=c`.  Also

\[
 \tau_{-t}=\tau_t                                      \tag{4}
\]

whenever both parameters occur.  The selected `v_t` endpoints for `t` and
`-t` are therefore antipodal.  For each `u_t`, (3) either cancels `u_t`
when `tau_t=-1`, or contains both `u_t` and `-u_t` when `tau_t=+1`.
The unpaired `v_1` is handled in exactly the same way by the missing affine
line point.  Thus

\[
                         d(L,M,c)(x)=d(L,M,c)(-x).       \tag{5}
\]

The origin is in this discrepancy.  Moreover,

\[
                         d(L,M,c)=c\,d(L,M,1),           \tag{6}
\]

because every endpoint scales by `c` and the Paley edge sign is unchanged
by nonzero scalar multiplication.

## 2. Kernel selectors

For every projective functional `N`, choose

\[
 S_N=\{0\}\ \cup\ \{\hbox{one representative of each pair }
                         \{x,-x\}\subset\ker N\}.       \tag{7}
\]

It has `m` points.  Define

\[
 q_N(L,M)=\langle\mathbf 1_{S_N},d(L,M,c)\rangle.       \tag{8}
\]

By (5), the pairing counts the origin plus whole antipodal pairs in
`ker N`; consequently it is independent of the representatives in (7).
By (6), dilation merely permutes those pairs, so (8) is also independent of
the nonzero centre `c`.

There is an explicit formula.  For `t != -1`, put

\[
                         N(t)=-tL+(t+1)M.                \tag{9}
\]

The only selected-half endpoints on `ker N(t)` are the origin and the pair
containing `u_t`.  The affine-line contribution gives

\[
 q_{N(t)}(L,M)=\mathbf 1_{\{\tau_t=-1\}}.               \tag{10}
\]

The two projective directions missed by (9) are `L` and `L-M`.  Directly,

\[
 q_L(L,M)=1,
 \qquad
 q_{L-M}(L,M)=\mathbf 1_{\{\tau_1=-1\}}.                \tag{11}
\]

For the second identity, `ker(L-M)` contains all `v_t`; the `2r` ordinary
parameter pairs `t,-t`, the origin, and the unpaired `t=1` term
give

\[
 1+2r+\mathbf 1_{\{\tau_1=+1\}}
       =\mathbf 1_{\{\tau_1=-1\}}\pmod 2.               \tag{12}
\]

Finally, the exact Paley identities are

\[
 \tau_0=\epsilon_M,\qquad \tau_1=\epsilon_L,
 \qquad
 \prod_{t\ne-1}\tau_t=\epsilon_M\epsilon_L.            \tag{13}
\]

Indeed, after `t=0` and the unpaired `t=1` are removed, (4) pairs every
remaining factor with an equal factor.  Multiplying the signs encoded by
(10)--(11) over all `p+1` projective directions yields the key identity

\[
 (-1)^{\sum_Nq_N(L,M)}
   =-\tau_1\prod_{t\ne-1}\tau_t=-\epsilon_M,
\]

or equivalently

\[
 \boxed{\ \sum_Nq_N(L,M)
       =\mathbf 1_{\{M\text{ is hard}\}}\pmod2.\ }      \tag{14}
\]

This is the selector theorem.  Notice that it is an all-direction identity;
the pairing for one preselected kernel need not have a fixed parity.

## 3. What the all-prime top ledger forces

For one half, its parallel count in direction `N` is odd exactly for
`N=L` and `N=M`; every other count is zero or two.  Therefore its raw
parallel profile satisfies

\[
 R_N\equiv \#\{i:L_i=N\}+\#\{i:M_i=N\}\pmod2.          \tag{15}
\]

At the top value

\[
 t=t_{\max}=4r^2-2r-5,
 \qquad s=(t+1)\bmod(p+1)=m+2,                           \tag{16}
\]

the balanced final quotas have multisets

\[
 \begin{array}{c|c}
 \text{hard}&(2r)^{m-2},(2r+1)^2\\
 \text{opposite}&(2r+1)^3,(2r+2)^{m-3}.
 \end{array}                                             \tag{17}
\]

Put `v=P mod 2 + 1_hard`.  Its support has `m-2` hard and three opposite
directions, hence weight `m+1`.  At `j=0` the exact Hamming/support ledger
has one fixed antipodal edge and no unused doubled orbit.  If its direction
is `F`, parallel parity gives the auxiliary incidence parity

\[
                              a=v+e_F.                   \tag{18}
\]

The parity support of `m` auxiliary occurrences has weight at most `m`.
Thus `F` lies in `supp(v)`, `wt(a)=m`, and all `m` auxiliary directions are
distinct.  In the opposite-fixed branch, they consist of exactly `m-2`
hard and two opposite directions.  In particular,

\[
                  \sum_N\sum_{i=1}^{m}q_N(L_i,M_i)
                     =m-2=2r=0\pmod2.                   \tag{19}
\]

Distinct auxiliary directions also make the `t=0` origin-edge orbits of
the `m` halves distinct: equality of two such orbits would force equality of
their projective auxiliary kernels.  The sole `j=0` cancellation unit is
therefore nonorigin.  These SDR and nonorigin facts are consequences of the
top ledger, not additional search assumptions.

This is precisely the first `j=0` residue left open by the earlier
parallel-parity endpoint theorem: that theorem excludes `4<=s<=m+1`, while
the top value has `s=m+2`.  The selector identity adds the missing
target-boundary information at this adjacent residue.

## 4. The adaptive contradiction

Let the sole cancellation unit occur on a nonfixed inversion orbit
`{e,-e}` whose two endpoints are nonzero.  Whether its multiplicities are
`1:1` or `2:1`, the difference between the raw half boundary and the final
boundary is

\[
                         C_e=\partial e+\partial(-e).    \tag{20}
\]

Its support consists of two nonzero antipodal point pairs.  Each pair lies
on exactly one projective kernel, hence

\[
                         \sum_N\langle S_N,C_e\rangle=0
                         \pmod2.                         \tag{21}
\]

For the antipodal fixed edge `f={z,-z}`, exactly one projective kernel
contains its endpoints.  Therefore

\[
                         \sum_N\langle S_N,\partial f\rangle=1
                         \pmod2.                         \tag{22}
\]

For each direction define the necessary boundary pairing

\[
 B_N=\sum_iq_N(L_i,M_i)+\langle S_N,C_e\rangle
                       +\langle S_N,\partial f\rangle.  \tag{23}
\]

If a common graph in this construction family existed, its compact target
atoms would have zero vertex boundary, so every `B_N` would be zero.
But (19), (21), and (22) give

\[
                         \sum_NB_N=1\pmod2.              \tag{24}
\]

At least one projective kernel selector is therefore a left-null
contradiction witness.  The witness direction can depend on the design;
requiring the fixed direction itself to work is false.

More generally, arbitrary simplification of overlapping half coefficients
changes the mod-two boundary by whole inversion-orbit pairs.  Every
nonorigin pair contributes zero after summing over all selectors, while an
origin-containing pair contributes one.  If `h` is the number of hard
auxiliaries, `f_0` the number of fixed edges, and `z_0` the parity of the
origin-containing orbit-pair adjustments (including unused origin doubles),
the aggregate necessary identity is

\[
                              h+f_0+z_0=0\pmod2.          \tag{25}
\]

Equation (25) must not be used to exclude arbitrary `j>0` completions:
origin doubles can change `z_0`.  At the `j=0` top endpoint above,
`h=m-2` is even, `f_0=1`, and auxiliary distinctness proves `z_0=0`, giving
the contradiction.

## 5. Exact p31 finite regression

The independent catalogue replay used

```text
/tmp/resii_auxiliary_component_10000_v1.json
sha256 196588c21a37c7788565b64c5b2a7dbfcafaedbd864dadf7e51b8b278895ae5b
```

It contained 2,969 designs admitting a clean collision seed and 148 distinct
scaled halves.  For all 148 halves, all thirty nonzero centres gave the same
pairing on every one of the 32 projective kernel selectors, and every tested
discrepancy was centrally symmetric.  The preselected fixed-direction sum
split almost exactly in half:

```text
fixed-direction half-sum even: 1486 designs
fixed-direction half-sum odd:  1483 designs
```

Thus a fixed-selector claim would be demonstrably wrong.  Allowing the
selector direction to adapt caught every design.  The number of obstructing
directions per design was:

```text
 7:   2
 9:  14
11: 139
13: 555
15: 749
17: 866
19: 442
21: 152
23:  49
25:   1
```

There were zero uncaught designs.  Independently, the full centre/fixed-edge
GF(2) systems for the same 2,969 records were all inconsistent, with rank
pair `225/226` in every case.  The catalogue itself was not component-
exhausted, so those counts are regression evidence only; the proof of the
whole stated endpoint branch is (14)--(24).

## 6. Scope boundary

This theorem uses all of the following structural hypotheses:

1. a prime `p=4r+3` in the balanced all-active branch-C range and the top
   value `t=t_max` at cancellation offset `j=0`;
2. `m=(p+1)/2` direction-localized Mobius halves, one for each hard target
   direction, with all hard centres nonzero;
3. the exact balanced top quotas (17);
4. an opposite fixed edge.

The prior support and parallel-parity theorems then force one fixed edge, no
unused doubled orbit, the `m-2` hard plus two opposite auxiliary SDR, and one
nonorigin cancellation unit.  The present theorem excludes all auxiliary
choices, relative scales, centres, and overlap realizations within that
branch, including both possible one-unit multiplicities.  It does not
establish that an arbitrary residual-(ii) source admits this localized-
Mobius representation.  The correct repository status is:

```text
j=0 top localized-Mobius opposite-fixed branch: CLOSED by adaptive GF(2)
arbitrary antisymmetric lifts and residual (ii): OPEN
```
