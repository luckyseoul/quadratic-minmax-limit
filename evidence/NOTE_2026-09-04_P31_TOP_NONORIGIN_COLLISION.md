# Exact pair loci and a p31 top triple-overlap countermechanism

Date: 2026-09-04

**Status:** proved a narrowly scoped collision theorem for the localized-
Mobius top endpoint.  With distinct auxiliary directions, an isolated
opposite pair collision lies on one of two exact algebraic loci and the pair
shares only one inversion orbit.  However, the sole top cancellation unit
need not come from an isolated pair: three hard-target halves with distinct
auxiliaries can meet on one nonorigin orbit with orientation multiplicities
`2:1`.  This invalidates pair-only exact models.  It does not construct all
sixteen halves, meet the transverse target, close branch C, or close residual
(ii).

## 1. Pair normalization

For two distinct target directions, put

\[
 X=L_1/j_1,\qquad Y=L_2/j_2.
\]

Normalize the two projective auxiliary directions as

\[
 M_1\sim \alpha X+Y,\qquad M_2\sim X+\beta Y.             \tag{1}
\]

They are distinct precisely when

\[
                         \alpha\beta\ne1.                 \tag{2}
\]

In the existing two-half normal form this means

\[
 A=q+\alpha,\qquad B=r+\beta.                             \tag{3}
\]

An opposite-orientation common edge has either the direct or the swapped
endpoint matching.

### Direct matching

The remaining endpoint equations are

\[
 t+rs=0,\qquad qt+s=0.
\]

The zero solution would force `alpha=beta=1`, contrary to (2).  Therefore
the direct locus is exactly

\[
 qr=1,\qquad \alpha+\beta=2,                              \tag{4}
\]

with

\[
 r=q^{-1},\quad
 t={\alpha-1\over q+1},\quad
 A=q+\alpha,\quad B=r+\beta.                              \tag{5}
\]

Its common nonorigin edge and parallel direction are

\[
 \left\{(1,-1),(t,qt)\right\},\qquad
 C\sim(1+q\alpha)X+(q+\beta)Y.                            \tag{6}
\]

If `C~cX+Y`, (6) recovers

\[
                              q={1-c\beta\over c-\alpha}, \tag{7}
\]

so fixed target, auxiliary, and cancellation directions leave at most one
direct candidate.  The exceptional projective cases are already covered by
(6); (7) is only the affine chart `C~cX+Y`.

### Swapped matching

Put

\[
 u=\alpha+1,\qquad v=\beta+1,\qquad z=qr.
\]

Eliminating `q,r,A,B` gives the exact quadratic

\[
 uvz^2-(u+v+1)z+1=0,                                     \tag{8}
\]

or equivalently

\[
 (\alpha+1)(\beta+1)z^2-(\alpha+\beta+3)z+1=0.           \tag{9}
\]

Every admissible root uniquely recovers

\[
 q=1-uz,\quad r=1-vz,\quad A=q+\alpha,\quad B=r+\beta.   \tag{10}
\]

The common edge and its parallel direction are

\[
 \left\{(1,1/r),(-1/q,-1)\right\},\qquad
 C\sim q(r+1)X-r(q+1)Y.                                  \tag{11}
\]

Thus fixed target and auxiliary data leave at most two swapped candidates,
before (11) imposes the prescribed cancellation direction.  The only point
with two opposite-orientation common orbits is the already-known rigid point

\[
 q=r={1\over2},\qquad A=B={3\over2}.
\]

There `alpha=beta=1`, so `M_1=M_2`; top auxiliary distinctness excludes it.
Consequently every distinct-auxiliary pair shares at most one
opposite-orientation orbit.

## 2. Why a pair-only top model is incomplete

One cancellation unit at an orbit `O` is

\[
 \kappa_O={n_O-|c_O|\over2}=1,
 \qquad |c_O|\le1.
\]

There are exactly two possibilities:

1. `n_O=2`, orientation multiplicities `1:1`, and final coefficient zero;
2. `n_O=3`, orientation multiplicities `2:1`, and final coefficient `+/-1`.

The second alternative is not merely formal.

Fix independent source endpoints `x,y`.  For every
`theta not in {0,-1}`, define functionals by

\[
 \begin{array}{c|cc}
       &x&y\\ \hline
 L_\theta&1&\theta\\
 M_\theta&\theta/(\theta+1)&\theta.
 \end{array}                                              \tag{12}
\]

Equivalently, in the dual basis,

\[
 L_\theta=(1,\theta),\qquad
 M_\theta=\left({\theta\over\theta+1},\theta\right)
            \sim(1,\theta+1).                            \tag{13}
\]

The parameter-`k` edge of this half is

\[
 E_\theta(k)=
 \left\{
 \left({\theta+1\over k+1},
       {k-\theta\over\theta(k+1)}\right),
 \left(0,{k\over\theta}\right)
 \right\},\qquad k\ne-1.                                \tag{14}
\]

At `k=theta`, this is the common edge `{x,y}`.

For distinct `theta,phi`, a same-orientation direct comparison first gives

\[
                         (\phi-\theta)(z-1)=0,
\]

so `z=1` and the common edge is exactly `{x,y}`.  In the opposite
orientation, the remaining equation is

\[
                         (\theta+\phi)(1-z^2)=0.          \tag{15}
\]

The cases `z=1`, `z=-1`, and `phi=-theta` force respectively
`theta=-1`, `phi=-1`, or both Mobius parameters to the excluded value `-1`.
A swapped matching is impossible because the affine endpoint in (14) has
nonzero first coordinate.  Hence two distinct members of (12) share exactly
the one inversion orbit of `{x,y}` and no other orbit.

## 3. Explicit p31 ternary triple

Take

\[
 x=(1,0),\qquad y=(0,1),\qquad
 (\theta_1,\theta_2,\theta_3)=(1,2,3).                   \tag{16}
\]

The targets

\[
 (1,1),(1,2),(1,3)
\]

are all hard.  The auxiliary directions are

\[
 (1,2),(1,3),(1,4),                                      \tag{17}
\]

which are distinct and have signs hard, hard, opposite.  Use centers
`(+1,+1,-1)`.  On the common orbit the three localized trades have
coefficients

\[
                              (-1,-1,+1),                 \tag{18}
\]

whose sum is `-1`.  Every other orbit occurs in exactly one trade.  The full
three-trade sum is therefore ternary, and

\[
 3(p-1)=90\quad\longrightarrow\quad 88
\]

used inversion orbits: exactly one cancellation unit.  The common edge is
nonorigin and its parallel direction is

\[
                              C=(1,1),                    \tag{19}
\]

which is hard at `p=31`.  The local auxiliary sign counts `2 hard + 1
opposite` fit inside either allowed top auxiliary sign profile; this is only
a type-count observation, not a sixteen-half extension theorem.

## 4. Parallel pullback and exact scope

Whether the unique cancellation is a `1:1` pair or a `2:1` triple, its raw
parallel ledger is

\[
 \sum_iP_D(E_i)=P_D^{\rm target}-{\bf1}_{D=F}
                              +2{\bf1}_{D=C}.             \tag{20}
\]

After Paley weighting,

\[
 \sum_iS(L_i,M_i)=4-p-\epsilon_F+2\epsilon_C.            \tag{21}
\]

At `p=31` the four `(F,C)` type cases are `-26,-30,-24,-28` for
`(hard,hard)`, `(hard,opposite)`, `(opposite,hard)`, and
`(opposite,opposite)` respectively.  Equation (21) does not distinguish the
pair and triple mechanisms.

The exact next model must therefore allow both alternatives and impose the
actual transverse cells of all sixteen halves.  The theorem here proves no
full top common graph and no all-prime obstruction.

Replay with

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_p31_top_nonorigin_collision.py

PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python \
  src/e1_gmin_m4_p31_top_nonorigin_collision.py
```

Branch C, residual (ii), E1, and `L=1/2` remain **OPEN**.
