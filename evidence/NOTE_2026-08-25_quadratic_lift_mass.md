# Minimum mass of a nonzero quadratic lift

Date: 2026-08-25. This is Proposition 15.642. It sharpens the affine
boundary reduction of Proposition 15.632. It does **not** close the
non-Walsh residual, Type I, R1, global QVAR, or the limit.

## 1. Exact stabilizer lemma

Let `p` be odd, `m=(p+1)/2`, and let `B` be a nonzero,
nonnegative, integer-valued polynomial of degree at most two on the middle
slice `J(p,m)`. A self-contained stabilizer argument gives

\[
 \mathbb E B\ge
 \begin{cases}
  1/p,&p\equiv3\pmod4,\\[2mm]
  r/((r+1)p),&p=4r+1.
 \end{cases}                                             \tag{1}
\]

Choose `X_0` with `B(X_0)>=1`, and average `B` under the stabilizer
of `X_0`. The result is a quadratic `q(t)`, where
`t=|X intersect X_0|` lies in `{1,...,m}`, with `q(t)>=0` and `q(m)>=1`.
For two uniform middle sets,

\[
 \mathbb Et={m^2\over p},\qquad
 \mathbb E[t(t-1)]={m^2(m-1)^2\over p(p-1)}.             \tag{2}
\]

Put `a=floor((m-1)/2)`. Matching the three moments of every
quadratic gives the following exact dual certificates.

If `p=4r+3`, then

\[
 \mathbb E q(t)= {2(2r+1)\over4r+3}q(r+1)
                 +{1\over4r+3}q(m).                     \tag{3}
\]

If `p=4r+1`, then

\[
 \mathbb E q(t)=
 {r(2r+1)\over(r+1)(4r+1)}q(r)
 +{2r+1\over4r+1}q(r+1)
 +{r\over(r+1)(4r+1)}q(m).                              \tag{4}
\]

All displayed weights are nonnegative, so (1) follows immediately. It is
sharp for the real-valued quadratic

\[
 q_*(t)={(t-a)(t-a-1)\over(m-a)(m-a-1)},                \tag{5}
\]

which is nonnegative at every integer `1<=t<=m`, vanishes at the
two consecutive central integers, and equals one at `m`. Sharpness here
is for the stabilizer-averaged real quadratic relaxation; no assertion is
made that every `q_*` lifts to a global integer-valued `B`.

## 2. Exact polynomial-distance reinforcement

Lemma 2 of Amireddy--Behera--Srinivasan--Sudan applies to every nonzero
degree-two polynomial on `J(p,m)`, for every `p>=5`, and gives

\[
 \Pr(B\ne0)\ge {\binom{p-4}{m-2}\over\binom pm}
 ={p^2-1\over16p(p-2)}.                                  \tag{6}
\]

Here `B` is nonnegative and integer-valued, so `B>=1` on its support and
therefore `E B >= Pr(B!=0)`. This exact all-parameter support lemma, not the
paper's stronger asymptotic main theorem, is what is used below. Combining
it with (1), the certified mass floor is the larger of the stabilizer floor
and (6).

## 3. Even scaled cost

In the parity decomposition of Proposition 15.632, write

\[
 A=P+2B,
\]

where `P` is the prescribed `{0,1}`-valued pointwise parity baseline and
`B` is a nonnegative integer-valued quadratic. A nonzero `B` adds

\[
 \Delta a=4p\,\mathbb E B.                              \tag{7}
\]

Both `a` and the baseline cost are even integers. Combining this with the
stronger of (1) and (6), a nonzero lift costs at least the even ceiling of
that lower bound. Sample floors are

\[
\begin{array}{c|rrrrrrrr}
p&5&7&11&13&17&19&23&31\\ \hline
\Delta a&2&4&4&4&6&6&8&10.
\end{array}
\]

In particular, for every `p>=5`, four nonzero lifts cost strictly more than
`p+1`: already (6) gives

\[
4\Delta a\ge16p\,{p^2-1\over16p(p-2)}
={p^2-1\over p-2}>p+1.                                  \tag{8}
\]

## 4. Boundary `D={infinity,v}`

Now take the residual size `|H|=4p+1` and suppose the odd-degree boundary
is exactly `D={infinity,v}`. In every affine direction, the finite point
`v` occupies one distinguished fibre `s_d(v)`. Proposition 15.632's
product formula reduces to

\[
 A_d(x)\equiv
 \begin{cases}
  x_{s_d(v)},&c_H=+1,\\
  1-x_{s_d(v)},&c_H=-1
 \end{cases}\pmod2.                                     \tag{9}
\]

Each quadratic direction type contains `m` directions and has exact
budget `m(p+1)`.

For `c_H=+1`, the baseline `x_{s_d(v)}` costs `p+1` in every
direction and therefore consumes the entire type budget. Thus there is no
lift at all:

\[
 \boxed{A_d(x)=x_{s_d(v)}\quad\hbox{pointwise for every }d.} \tag{10}
\]

For `c_H=-1`, the baseline `1-x_{s_d(v)}` costs `p-1` per direction,
leaving exactly `p+1` scaled units in each type. Consequently (8) proves
that **each type has at most three directions with a nonzero quadratic
lift**, uniformly for every `p>=5`. At `p=7`, the exact even floor improves
this to two.

Equation (10) and the uniform three-exception bound are an all-prime
rigidity/sparsity reduction, not an exclusion of the boundary. The next
attack must combine the pointwise
quadratic identities across projective directions, or control the sparse
exceptional set in the negative-product branch.

## 5. Literature and OEIS check

The exact support bound (6) is Lemma 2 in
Amireddy--Behera--Srinivasan--Sudan,
[A Near-Optimal Polynomial Distance Lemma over Boolean Slices](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2025.11).
Their main theorem improves the balanced-slice support asymptotically toward
`1/4`, but uses an unspecified absolute exponent and is not used as an
effective finite-`p` statement here. The independent stabilizer certificate
(3)--(4) is stronger at the smallest primes. The resulting formulas are
closed rational functions rather than an unexplained integer sequence, so
no OEIS sequence claim is made.

## Reproduction

- exact formulas and machine-readable certificates:
  `src/e1_gmin_m4_prop15642.py`;
- tests of all odd `p<=103`: `tests/test_prop15642.py`;
- generated evidence: `evidence/e1_gmin_m4_prop15642.json`.
