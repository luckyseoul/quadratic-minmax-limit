# Degree-four bivector preordering no-go

**Status:** proved all-orders method barrier for the natural squared-row
encoding.  This note does not prove or disprove the outgoing-half theorem and
does not close the MathOverflow limit.  It strengthens the single-row moment
obstruction in Proposition 6.5d: even after adding arbitrary affine SOS
localizers and every degree-four product of two distinct bivector rows, the
uniform Boolean pseudoexpectation remains feasible for every `n>=45`.

Let `A` be an order-`n` signing, let `M=Phi(A)`, and retain the notation

\[
 b_s=b_{xy}={x\wedge y\over2}\in\{0,\pm1\}^{N},\qquad
 d_s=M-{|Q_A(x)+Q_A(y)|\over2},\qquad N={n\choose2},
\]

where `s=(x,y)`.  At the critical doubling threshold put

\[
 t=(\sqrt2-1)M,\qquad w_s=t+d_s,
 \qquad g_s(r)=w_s^2-\langle b_s,r\rangle^2.                 \tag{1}
\]

Because `d_s>=0`, the exact orientation system is

\[
 r_e^2=1\quad(e\in[N]),\qquad g_s(r)\ge0\quad(s=(x,y)).     \tag{2}
\]

Indeed, `g_s>=0` is exactly
`|<b_s,r>|<=t+d_s`.  Thus (2), not a relaxation of the row inequalities,
is the exact squared formulation of the nonlinear cover gate.

## 1. The degree-four Pluecker-aware preordering

Consider the full degree-four preordering in this squared-row encoding.  A
refutation in this specified system would be an identity

\[
 -1=\sigma_0+
       \sum_s\sigma_s g_s+
       \sum_{s\le u}\lambda_{su}g_sg_u+
       \sum_e h_e(r)(r_e^2-1),                             \tag{3}
\]

where

- `sigma_0` is a sum of squares of degree at most four;
- each `sigma_s` is a sum of squares of affine polynomials, so
  `deg(sigma_s g_s)<=4`;
- `lambda_su>=0` and the row-product terms have degree four; and
- the Boolean-ideal terms have degree at most four.

The coefficients `b_s` in (3) are the actual decomposable Boolean
bivectors.  Consequently every Pluecker identity among their coordinates,
and every cross-row identity derived from those relations, is already
available when forming (3).  The terms `g_sg_u` are the first level at which
two different row inequalities interact in this squared encoding.

## 2. Exact one-row localizing spectrum

Let `E_0` denote expectation for a uniformly random
`r in {+1,-1}^N`.  Fix a row `b`, write

\[
 S=\langle b,r\rangle,qquad m=\|b\|_2^2=|\operatorname{supp}b|,
 \qquad \delta=w^2-m.
\]

For every affine polynomial `q=a_0+sum_e a_e r_e`, the localizing form
`E_0[(w^2-S^2)q^2]` has no constant--linear cross term.  Its constant block
is `delta`.  On inactive coordinates its linear block is `delta I`, while
on `supp(b)` it is

\[
 (\delta+2)I-2bb^T.                                       \tag{4}
\]

The exact eigenvalues are therefore

\[
 \delta,\qquad \delta+2=w^2-m+2,
 \qquad w^2-3m+2,                                         \tag{5}
\]

with the evident multiplicities.  In particular

\[
 E_0[gq^2]\ge0\quad\hbox{for every affine }q
 \quad\Longleftrightarrow\quad w^2\ge3m-2                \tag{6}
\]

when `m>=1` (and the zero row is immediate).  This is the exact degree-four
localizing threshold, not a Khintchine estimate.

## 3. Exact cross-row product

For two rows `b,c`, put

\[
 m_b=\|b\|_2^2,\qquad m_c=\|c\|_2^2,\qquad
 \gamma=\langle b,c\rangle,\qquad
 h=|\operatorname{supp}b\cap\operatorname{supp}c|,
\]

and `delta_b=w_b^2-m_b`, `delta_c=w_c^2-m_c`.  Direct fourth-moment
expansion on the Boolean cube gives

\[
 \boxed{
 E_0[g_b g_c]
   =\delta_b\delta_c+2\gamma^2-2h.}                        \tag{7}
\]

For the decomposable rows, `gamma` is the determinant pairing

\[
 \gamma={
 (x\mathbin\cdot u)(y\mathbin\cdot v)
 -(x\mathbin\cdot v)(y\mathbin\cdot u)\over4},           \tag{8}
\]

and `h` is the intersection size of the two cut supports.  Thus (7) retains
the exact cross-row bivector geometry rather than replacing the rows by
arbitrary vectors.  The only negative contribution is the Boolean
same-coordinate correction `-2h`.

## 4. Critical widths dominate the whole degree-four system

A Boolean bivector supported on the cut where `x_i y_i` changes has

\[
 m_s=d_H(x,y)(n-d_H(x,y))\le m_*:=\left\lfloor{n^2\over4}\right\rfloor.
                                                                    \tag{9}
\]

Proposition 5.2 gives `M>=n sqrt(n-1)/pi`, and hence

\[
 t^2\ge {3-2\sqrt2\over\pi^2}n^2(n-1).                  \tag{10}
\]

For every `n>=45`,

\[
 {3-2\sqrt2\over\pi^2}(n-1)>{3\over4}.                 \tag{11}
\]

For a completely rational check of the endpoint, use
`sqrt(2)<5657/4000`, `pi^2<10`, and therefore

\[
 44{3-2\sqrt2\over\pi^2}
 >44{343\over20000}={15092\over20000}>{3\over4}.
\]

Equations (9)--(11), together with `w_s>=t`, show

\[
 w_s^2>3m_*\ge3m_s.                                      \tag{12}
\]

Thus every one-row affine localizer in (3) is positive semidefinite by
(6).  Moreover

\[
 \delta_s=w_s^2-m_s>2m_*,
\]

so (7) and `h<=m_*` give

\[
 E_0[g_sg_u]
 >4m_*^2-2m_*>0.                                         \tag{13}
\]

The uniform functional is an actual expectation, so it is nonnegative on
`sigma_0` and annihilates every Boolean-ideal term.  It is nonnegative on
all remaining terms of (3) by (6) and (13).  Applying it to (3) would give
`-1>=0`, a contradiction.

Hence:

\[
 \boxed{\text{For every }n\ge45\text{, system (2) has no degree-four
 squared-row preordering refutation of the form (3).}}     \tag{14}
\]

This conclusion holds whether (2) is feasible or infeasible and for every
signing `A`; it is an all-orders tail obstruction to this proof system.

## 5. Verdict

The changed premise over Proposition 6.5d was to use genuine degree-four
information: affine SOS multipliers for each exact row constraint, all
pairwise products of distinct constraints, and the exact intersection
geometry of decomposable Boolean bivectors.  It still does not reach the
critical scale.  At that scale `w_s^2` has order `n^3`, whereas every row
has only order `n^2` support, leaving the uniform degree-four pseudomodel
strictly interior from order 45 onward.

This does **not** rule out higher degree, a different lifted encoding, or a
constructive rounding theorem.  It does show that a continuation through
the natural squared-row SOS/preordering must use degree at least six, and a
fixed-degree continuation would need information not already dominated by
the critical `n^3` widths.  The outgoing-half theorem and the original MO
limit remain open.
