# The mod-p Paley completion code and the mismatch-one tangent

**Date:** 2026-09-02
**Status:** proved all-prime structure and route barrier; residual (ii)
remains open

Let (q=p^2), (n=q+1), and let (C) be the symmetric Paley conference
matrix on (mathbf P^1(mathbf F_q)).  This note determines the reduction of
(C) modulo (p), its first Bockstein, and what those data can say about an
all-bad three-spike completion.  The result is an exact obstruction
coordinate, but not an exclusion: the bounded global alphabet, not the
linear code, remains the missing input.

## Rank and Smith form

Use homogeneous representatives (r_P=(X_P,Y_P)), and put

\[
 N={q-1\over2},\qquad m={p+1\over2}.
\]

Over (mathbf F_q), and hence after scalar extension from (mathbf F_p),

\[
 \overline C_{P,Q}
 =\det(r_P,r_Q)^N
 =\sum_{k=0}^N(-1)^k{N\choose k}
 X_P^{N-k}Y_P^kY_Q^{N-k}X_Q^k.                 \tag{1}
\]

The two base-(p) digits of (N) are both ((p-1)/2).  Lucas' theorem says
that the coefficient in (1) is nonzero exactly when

\[
 k=k_0+pk_1,\qquad 0\le k_0,k_1\le{p-1\over2}.
\]

There are (m^2) such exponents.  Their homogeneous evaluation vectors are
independent: after setting (Y=1), a linear relation would be a polynomial
of degree at most (N<q) vanishing at every element of (mathbf F_q).
Thus

\[
 \boxed{\operatorname {rank}_{\mathbf F_p}\overline C=m^2.} \tag{2}
\]

Since (C^2=p^2I), every Smith invariant of (C) divides (p^2), and
(|\det C|=p^n).  Equation (2) counts the unit invariants; the determinant
then counts the other two types.  Therefore

\[
 \boxed{\operatorname {SNF}(C)=
  1^{m^2}\oplus p^{(p-1)^2/2}\oplus(p^2)^{m^2}.}             \tag{3}
\]

This is compatible with, but different from, Proposition 15.629's Smith
form for the Gram matrix of the integral (+p) eigenlattice.

## The first Bockstein

Put

\[
 K=\ker\overline C,\qquad I=\operatorname {im}\overline C,
 \qquad \mathcal H=K/I.
\]

The identity (overline C^2=0) gives (I\subset K).  Symmetry gives
(K^\perp=I), so the ordinary dot product descends to a nondegenerate form
on (mathcal H).  From (2),

\[
 \dim\mathcal H=n-2m^2={(p-1)^2\over2}.                    \tag{4}
\]

For (z\in K), take an integral lift, write (Cz=pw), and define

\[
 \tau[z]=[w]\quad\hbox{in }\mathcal H.                    \tag{5}
\]

Changing the lift changes (w) by an element of (I), so (5) is
well-defined.  The identities (C^2=p^2I) and (C^T=C) show that
(	au^2=1) and that (	au) is self-adjoint.

Let (L_\pm=\ker_{\mathbf Z}(C\mp pI)).  Reduction embeds
(L_\pm/pL_\pm) into (K), with dimension (n/2).  Their intersection is
exactly (I).  Indeed, (pe_i\pm Ce_i\in L_\pm), so (I) lies in the
intersection.  Conversely, if (a\in L_+), (b\in L_-), and
(a-b=pz), then (Cz=a+b\equiv2a\pmod p), so their common residue lies in
(I).  Consequently

\[
 \mathcal H=\mathcal H_+\mathbin\perp\mathcal H_-,\qquad
 \dim\mathcal H_+=\dim\mathcal H_-={(p-1)^2\over4},        \tag{6}
\]

and these are exactly the (+1) and (-1) Bockstein eigenspaces.

## Why neither restriction nor Schur powers exclude mismatch one

Fix a square augmented affine line (Gamma).  Its incidence word, and the
incidence word of every augmented line in a square direction, lies in
(L_+).  The circle sublattice already has a surjective integral restriction
map

\[
 L_+\longrightarrow\mathbf Z^\Gamma.                       \tag{7}
\]

To see this, a line parallel to (Gamma) restricts to (e_\infty).  For
each finite (u\inGamma), a line through (u) in a second square direction
restricts to (e_\infty+e_u); subtracting the parallel line gives (e_u).
Thus every prescribed integral trace on (Gamma), including a trace with
three entries (3), one entry (-1), and all remaining entries (1), has
an exact integral (+p)-eigenvector extension.  Its off-circle entries need
not belong to the required small alphabet, which is precisely the point:
no congruence or local restriction theorem can supply that missing bound.

There is an equally sharp Schur-power barrier.  If (U_+) is the reduction
of (L_+), then

\[
 \boxed{\operatorname {span}(U_+\ast U_+)=\mathbf F_p^n.}   \tag{8}
\]

Indeed, the coordinatewise product of two augmented square lines in
different directions through (u) is (e_\infty+e_u).  Summing these over
the (p) finite points of one line and subtracting its incidence word
isolates (e_\infty) modulo (p), and then every (e_u).  Hence a
Lucas/Reed--Muller or AG-code argument that only places a Schur product in a
linear evaluation code has no nonzero parity check left; higher Schur powers
are full as well.

## Exact mismatch-one coordinate

Now work in the signed gauge of the near-pencil reduction.  Let

\[
 v=\mathbf1+2\mathbf1_E,\qquad Kv=pv,
\]

where (|E|=3), and orient the unique square-circle word (d) so that
(Kd=pd), (operatorname {supp}d=Gamma(E)), and (d=1) on (E).
If (mu) is the number of negative entries of (d) on
(Gamma(E)\setminus E), then

\[
 \|v\|^2=p^2+25,\qquad \|d\|^2=p+1,
 \qquad \langle v,d\rangle=p+7-2\mu.                       \tag{9}
\]

It follows exactly that

\[
 \begin{aligned}
 \|v-5d\|^2&=p(p+15)+20(\mu-1),\\
 \langle v-5d,d\rangle&=-4p-2(\mu-1).                     \tag{10}
 \end{aligned}
\]

For the live primes (p\ge53), (10) gives the intrinsic code-theoretic
description

\[
 \boxed{\mu=1\iff [v-5d]\in[d]^\perp
 \text{ is isotropic in }\mathcal H_+.}                    \tag{11}
\]

This does not contradict (6).  The vector ([d]) is anisotropic because
(|d|^2\equiv1pmod p); its orthogonal complement in the nondegenerate
space (mathcal H_+) has dimension ((p-1)^2/4-1\ge3), and therefore
contains nonzero isotropic vectors over (mathbf F_p).  Thus the first
Bockstein sees mismatch one as a tangent point of its finite orthogonal
quadric, but that tangent locus is genuinely populated.

**Verdict.**  The all-prime rank, Smith form, Bockstein splitting, integral
restriction map, and Schur closure are now explicit.  They rule out a pure
mod-(p), higher-lift-without-bounds, or linear Schur-code exclusion of
(mu=1).  A close must use the simultaneous global alphabet
(v_i\in\{1,3\}) in the signed gauge (equivalently
(y_i\in\{\pm1,\pm3\}) before switching), together with the common graph
or another nonlinear condition; the exact remaining modular locus is the
isotropic tangent in (11).
