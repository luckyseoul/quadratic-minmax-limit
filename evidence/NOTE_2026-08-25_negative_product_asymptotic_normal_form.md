# Asymptotic normal form of the negative infinity-point boundary

Date: 2026-08-25. This is Proposition 15.644. It reduces, but does not
exclude, the residual branch `D={infinity,v}`, `c_H=-1` for all sufficiently
large odd primes. Small primes, other boundaries, residual (ii), Type I, R1,
global QVAR, and the limit remain open.

## 1. Exactly one exception of each type

Proposition 15.642 writes every directional slack as

\[
A_d=1-x_{s_d(v)}+2B_d,                                  \tag{1}
\]

where `B_d` is a nonnegative integer-valued polynomial of degree at most two
on `J(p,(p+1)/2)`. Each quadratic direction type has exact surplus `p+1`:

\[
\sum_{d:\epsilon_d=\tau}4p\,\mathbb E B_d=p+1.          \tag{2}
\]

The degree-two case of the near-optimal polynomial distance theorem of
Amireddy--Behera--Srinivasan--Sudan gives, with `q=(p-1)/2`,

\[
\Pr(B_d\ne0)\ge(q/p)^2(1-q^{-\varepsilon})              \tag{3}
\]

for some absolute `epsilon>0`, once `p` is sufficiently large. Since a
nonzero nonnegative integer-valued `B_d` is at least one on its support,
two nonzero lifts would cost at least

\[
8p(q/p)^2(1-q^{-\varepsilon})>p+1                       \tag{4}
\]

for all sufficiently large `p`, contradicting (2). The surplus is positive,
so each type has exactly one exceptional direction. In that direction (2)
forces

\[
4p\,\mathbb E B_d=p+1,qquad a_d=(p-1)+(p+1)=2p.         \tag{5}
\]

This use of the asymptotic theorem is qualitative: it proves existence of a
threshold but does not claim an explicit numerical cutoff.

## 2. Baseline parallel counts

Let `I` be the number of infinity edges, `E=4p+1-I` the number of finite
edges, and `P_d` the number of finite edges parallel to direction `d`. In a
baseline direction, the argument of Proposition 15.643 with target `4-z_j`
gives

\[
{2(I+P_d-3)\over p-1}\in\mathbb Z                     \tag{6}
\]

and the signed sum of the additive inter-fibre matrix gives

\[
|pP_d-E+2|\le E-P_d.                                    \tag{7}
\]

The lower half of (7) is automatic; the upper half is

\[
(p+1)P_d\le2E-2.                                        \tag{8}
\]

Because infinity lies in the odd boundary, `I>=1`, so `E<=4p` and (8)
gives `P_d<=7`. For sufficiently large `p`, `q>7`; equation (6) therefore
forces the same value `rho in {0,...,7}` in all `p-1` baseline directions.
Write

\[
I=3-\rho+qk_0.                                          \tag{9}
\]

## 3. Exceptional parallel counts

In an exceptional direction, the exact mean formula and (5) give a signed
transverse contribution

\[
5p-I-pP_d.
\]

Its absolute value is at most the `E-P_d` transverse edges. The two sides of
that inequality yield

\[
1\le P_d\le {9p+1-2I\over p+1}<9.                       \tag{10}
\]

Let `U,V` be the positive- and negative-type exceptional parallel counts.
Summing all direction counts and using (9) gives

\[
U+V=(8-k_0-2\rho)q+2+\rho.                              \tag{11}
\]

By (10), `2<=U+V<=16`. Hence, once `q>14`, the coefficient of `q` in
(11) must vanish:

\[
k_0=8-2\rho,qquad U+V=2+\rho.                           \tag{12}
\]

## 4. Unique surviving arithmetic profile

Substitution into (9) gives

\[
I=3-\rho+q(8-2\rho).                                    \tag{13}
\]

The infinity degree is odd, so `rho` is even. Nonnegativity leaves
`rho=0,2,4`. The value `rho=4` makes `I=-1`. If `rho=0`, then
`I=4p-1` and only two finite edges remain. The finite boundary of the
infinity star differs from the desired singleton `{v}` in at least `I-1`
vertices, whereas two finite edges can toggle at most four vertices. Thus
`rho=0` is impossible.

Only `rho=2` remains, and (12)--(13) become

\[
\boxed{I=2p-1,\quad E=2p+2,\quad P_d=2
\text{ in every baseline direction},\quad U+V=4.}       \tag{14}
\]

There are `q` baseline negative-type directions, contributing `2q` negative
Paley edges, an even number. Since `c_H=-1`, the negative-type exceptional
count `V` is odd. With `U,V>=1` and `U+V=4`,

\[
\boxed{(U,V)=(3,1)\text{ or }(1,3).}                    \tag{15}
\]

## 5. Interpretation and remaining obstruction

The number `I=2p-1`, two exceptional directions, and two edges in every
other direction strongly suggest the union of two nonparallel affine lines
through `v`: such a union has exactly `2p-1` points and meets every other
parallel class in one fibre of size one and `p-1` fibres of size two. This
is a geometric lead, not part of the proved classification. The next task is
to derive that two-line structure from the complete inter-fibre `l1`
constraints, then test whether its finite edge graph can satisfy all signed
baseline identities.

## 6. Literature check

The external input is [A Near-Optimal Polynomial Distance Lemma over Boolean
Slices](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2025.11).
The exact rank-resilience problem for the critical pair-inclusion matrix was
also compared with Grosu--Person--Szabo and Plaza--Xiang; their general
inclusion-matrix resilience theorems do not directly supply the diagonal
`n=2r+1` minimum-support statement. No duplicate of the arithmetic normal
form (14)--(15) was found.

## Reproduction

- arithmetic verifier: `src/e1_gmin_m4_prop15644.py`;
- generated evidence: `evidence/e1_gmin_m4_prop15644.json`;
- tests: `tests/test_prop15644.py`.
