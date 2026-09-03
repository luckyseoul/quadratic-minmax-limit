# Removed-pencil-line circle: exact \(\mu=1\) cut audit

**Date:** 2026-09-02
**Status:** proved local reduction and explicit obstruction; not a residual-(ii)
close

This note treats only the uniform affine exception
\(\Gamma=U_i\cup\{\infty\}\) left by
`NOTE_2026-09-02_CIRCLE_BOUNDARY_GEOMETRY.md`.  It uses no prime census.

Work in the signed chart of
`NOTE_2026-09-02_TWO_HALF_NEAR_PENCIL_REDUCTION.md`, so

\[
 v={\bf1}+2{\bf1}_E,\qquad Kv=pv,
\]

and orient the sparse \(+p\) circle vector \(d\), supported on \(\Gamma\),
to be positive on the three-set \(E\).  Put

\[
 N=\{j\in\Gamma\setminus E:d_j=-1\},\qquad
 P=\Gamma\setminus(E\cup N).
\]

Here \(|N|=1\), \(|P|=p-3\), and the circle signing gives
\(K_{uv}=d_ud_v\) for distinct \(u,v\in\Gamma\).  Thus every edge inside
\(E\cup P\) is positive and every circle edge incident with \(N\) is
negative.

## Complete circle-supported cut formula

The eigenvector equation gives the exact \(K\)-row sums

\[
 (K{\bf1})_u=
 \begin{cases}
 3p-4,&u\in E,\\
 p-6,&u\in P,\\
 p+6,&u\in N.
 \end{cases}                                             \tag{1}
\]

For \(T\subseteq\Gamma\), let

\[
 e=|T\cap E|,\qquad
 a=|T\cap(E\cup P)|,\qquad
 \nu={\bf1}_{N\in T}.
\]

Subtracting twice the signed internal edge sum from the row-sum total
proves

\[
 \boxed{W_K(T)=a(p-5-a)+2e(p+1)+\nu(p+6+2a).}             \tag{2}
\]

For \(p\ge53\), (2) is negative for exactly the following
circle-supported cuts (up to complementation):

\[
\begin{array}{c|c}
T&W_K(T)\\ \hline
P&-2p+6,\\
P\setminus\{u\},\ u\in P&-p+4,\\
P\cup\{e\},\ e\in E&8-p.
\end{array}                                               \tag{3}
\]

Indeed, when \(\nu=0\), the concave quadratic in (2) is negative only at
the last two admissible values of \(a\) for \(e=0\), and at the last value
for \(e=1\); its endpoint values are positive for \(e=2,3\).  When
\(\nu=1\), the added term makes both endpoint values nonnegative.

The lower half of the full Max box therefore gives, in particular,

\[
 Z_H(P)\le-(p-3),\quad
 Z_H(P\setminus\{u\})\le-{p-3\over2},\quad
 Z_H(P\cup\{e\})\le-{p-7\over2}.                         \tag{4}
\]

These are simultaneous restrictions, not three independent hereditary
families.

Let \(q=|P\setminus D|\).  The strict odd-boundary surplus lemma and the
first inequality in (4) imply an additional signed statement.  If
\(c_+(P)\) is the number of positive \(K\)-edges of \(H\) crossing \(P\),
then

\[
 \boxed{c_+(P)\le\left\lfloor{2s-q\over2}\right\rfloor,}
 \qquad s=7\text{ or }9.                                 \tag{5}
\]

To see this, write \(c=c_++c_-\) and \(Z=c_+-c_-\).  The boundary lemma
gives \(c\le |D\cap P|+2s=p-3-q+2s\), while (4) gives
\(Z\le-(p-3)\).  Hence \(2c_+=c+Z\le2s-q\).

## The phase gate does not remove the line

In the original Paley chart the canonical sparse vector on the square
line \(U_i\cup\{\infty\}\) is constant.  After a global normalization,
the \(\mu=1\) condition consequently says that the Boolean signing on
\(\Gamma\) is negative exactly at \(N\).

If \(j\ne i\) and
\(q_i=U_i\cap(b+U_j)\), the two-outlier boundary satisfies

\[
 D\cap\Gamma=\Gamma\setminus\{\infty,a,q_i\}.             \tag{6}
\]

Thus the necessary phase \(\prod_Dx=1\) factors exactly as

\[
 \boxed{\prod_{v\in D\setminus\Gamma}x_v
       =(-1)^{{\bf1}_{N\in D}}.}                          \tag{7}
\]

Equation (7) prescribes the outside-circle phase; it is not a contradiction.
In particular, neither the spike positions nor the value \(\mu=1\)
determines that outside bit.

## Exact directional-row check

Let \(L=\operatorname{dir}(U_i)\).  The removed line is the old pencil
fibre, while \(j=b+U_i\) is the replacement exceptional fibre.  The two
circle/complement Boolean eigenvectors give the quadratic-character profiles
\(x^+,x^-\in J(p,(p+1)/2)\); exactly one contains \(j\).  Proposition
15.758 writes the hard directional cell as

\[
 A_L(x)=x_j+2B_L(x),\qquad B_L(x)\ge0,
 \qquad S_H(y)=3+2A_L(x).
\]

Consequently

\[
 \boxed{S_H(d+w)+S_H(d-w)
 =8+4\bigl(B_L(x^+)+B_L(x^-)\bigr)\ge8.}             \tag{8}
\]

This is the same direction as the phase-pair lower bound, not an upper
contradiction.  It is genuinely locally feasible: on three nonexceptional
fibres with patterns \(001\) and \(110\), respectively, the compact atom

\[
 B_0=x_ax_b+x_c-x_ax_c-x_bx_c
\]

has value one at both profiles and has the exact Proposition 15.758
\(p+1\)-atom mass, offset, coefficient sum, and Eulerian Radon parity.
Copies of this atom accommodate every hard excess \(e_L\).  Thus neither
the one-row Radon data nor its parity can reverse (8); only common
cross-direction simple-graph compatibility can do so.

## Explicit local obstruction

Equations (3)--(5) are sharp enough to show why the requested ingredients
do not yet close this exception.  Match each required odd vertex of
\(P\cap D\) to a distinct outside-circle vertex along a negative \(K\)-edge.
If \(q\) is even, add \(q\) further negative crossing edges in pairs at
vertices of \(P\cap D\); if \(q\) is odd, add \(q+1\).  This makes every
vertex of \(P\) have the prescribed parity and gives

\[
 Z_H(P)=-p+3\quad(q\ {\rm even}),\qquad
 Z_H(P)=-p+2\quad(q\ {\rm odd}).                          \tag{9}
\]

The extra incidences can be kept bounded (at most three at a chosen
vertex), and the remaining vertices of \((E\cup N)\cap D\) can each be
matched negatively outside.  Such distinct negative outside neighbours
exist: (1), together with the known circle signing, leaves
\(\gg p\) negative outside neighbours at every circle vertex, whereas at
most \(p+1\) are requested.  For every circle-supported cut in (3), this
local graph satisfies (4); for all other circle-supported cuts its signed
cut is nonpositive while the upper bound in (2) is nonnegative.  The
opposite half of the Max interval is of order \(p^3\) and is automatic for
this \(O(p)\)-edge local piece.

This is a local feasibility obstruction, not a construction of the full
common graph: its outside endpoints still have to be joined to the fixed
outside boundary, exact direction quotas, total edge count, and total
signed score.  It proves that exact circle signing, the boundary parity and
surplus, the phase bit, and every Max cut supported on this one circle are
mutually compatible.  Excluding the removed-line \(\mu=1\) case therefore
requires a genuinely simultaneous cut using vertices outside \(\Gamma\)
(or an independent classification of the Boolean completion fibre).
