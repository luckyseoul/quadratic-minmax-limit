# All-prime negative-product infinity-point exclusion for `p>=17`

Date: 2026-08-25. This is Proposition 15.647. It replaces the asymptotic
input of Propositions 15.644--15.646 with exact all-prime arithmetic and
excludes `D={infinity,v}`, `c_H=-1` for every odd prime `p>=17`. Together
with Proposition 15.643, both edge-product signs of this two-point boundary
are now excluded for `p>=17`. The four primes `5,7,11,13`, other boundaries,
residual (ii), Type I, R1, global QVAR, and the limit remain open.

## 1. A global signed-mean identity

Let `I` and `E` be the numbers of infinity and finite edges, `P_d` the
number of finite edges parallel to direction `d`, and

\[
 W=\sum_{e\text{ finite}}C_e.
\]

Every finite edge parallel to a direction of type `epsilon_d` has sign
`epsilon_d`. Hence its signed transverse sum in direction `d` is
`W-epsilon_d P_d`. Proposition 15.632's exact directional mean becomes

\[
 a_d=I+(p+1)P_d-\epsilon_dW-3p.                         \tag{1}
\]

For two directions `d,e` of the same type, (1) gives

\[
 a_d-a_e=(p+1)(P_d-P_e).                                \tag{2}
\]

## 2. Exactly one exception of each type for every `p>=7`

In the negative-product branch,

\[
 A_d=1-x_{s_d(v)}+2B_d,
 \qquad r_d:=a_d-(p-1)\ge0.                             \tag{3}
\]

The exact type budget says

\[
 \sum_{d:\epsilon_d=\tau}r_d=p+1.                     \tag{4}
\]

Proposition 15.642 gives at most two exceptional directions per type at
`p=7` and at most three for every larger prime. Each type contains
`m=(p+1)/2` directions, so for every `p>=7` there is at least one baseline
direction `e` with `r_e=0`. Equation (2) then shows

\[
 r_d=(p+1)(P_d-P_e).                                    \tag{5}
\]

Every nonzero `r_d` is therefore a positive multiple of `p+1`. Combined
with (4), this proves that each type has exactly one exception and

\[
 r_d=p+1,\qquad a_d=2p.                                 \tag{6}
\]

This is exact for all odd primes `p>=7`; no near-optimal asymptotic slice
distance theorem is needed.

## 3. Baseline count parametrization

Let `x,y` be the common positive- and negative-type baseline parallel
counts. Equation (5) says the two exceptional counts are `x+1,y+1`.
Counting directions gives

\[
 E=m(x+y)+2,\qquad I=4p-1-m(x+y).                       \tag{7}
\]

The baseline additive-matrix comparison from Proposition 15.643 gives

\[
 {2(I+P_d-3)\over p-1}\in\mathbb Z.                    \tag{8}
\]

Substitute (7) into (8) first with `P_d=x`, then with `P_d=y`. With
`q=(p-1)/2=m-1`, the two conditions reduce to

\[
 q\mid y,\qquad q\mid x.                                \tag{9}
\]

## 4. Exclusion for `p>=17`

Since `I>=1`, equation (7) gives

\[
 x+y\le {4p-2\over m}<8.                                \tag{10}
\]

For `p>=17`, `q>=8`. Equations (9)--(10) force `x=y=0`.
Then (7) gives

\[
 E=2,\qquad I=4p-1.                                     \tag{11}
\]

The infinity star contributes odd degree at `I` finite vertices. To leave
only `v` in the finite odd boundary, the finite graph must toggle at least
`I-1=4p-2` vertices. Two finite edges can toggle at most four. This is
impossible, proving the branch empty for every odd prime `p>=17`.

The same arithmetic leaves only the following negative-product count pairs:

- `p=7`: `(x,y)=(0,3),(0,6),(3,0),(3,3),(6,0)`;
- `p=11`: `(x,y)=(0,5),(5,0)`;
- `p=13`: `(x,y)=(0,6),(6,0)`.

The prime `p=5` needs separate treatment because the uniform three-exception
bound does not guarantee a baseline direction in each three-direction type.

## Reproduction

- arithmetic verifier: `src/e1_gmin_m4_prop15647.py`;
- generated evidence: `evidence/e1_gmin_m4_prop15647.json`;
- tests: `tests/test_prop15647.py`.
