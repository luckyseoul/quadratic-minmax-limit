# Centered-circle `mu=1` saturation ledger

**Date:** 2026-09-02
**Status:** proved exact necessary condition; residual (ii) remains open

This note treats only the uniform centered finite-circle exception left by
`NOTE_2026-09-02_CIRCLE_BOUNDARY_GEOMETRY.md`.  It is not a finite census.

Use the gauge of equations (18)--(29) in
`NOTE_2026-09-02_TWO_HALF_NEAR_PENCIL_REDUCTION.md`.  Thus

\[
 K_{ij}=1\quad(i,j\in E),\qquad v={\bf1}+2{\bf1}_E,
 \qquad Kv=pv .
\]

Suppose that the unique square circle through \(E\) is the centered norm
circle

\[
 \Gamma=\{z:N(z)=R\},\qquad \eta(R)=-1,
\]

and that its oriented sparse word has exactly one mismatch \(n\).  Put

\[
 P=\Gamma\setminus(E\cup\{n\}),\qquad a=|P|=p-3.
\]

The sparse-circle signing makes every edge inside \(P\) positive in \(K\).
Equation (24) gives the exact ambient cut weight

\[
 W_K(P)=-2a.                                             \tag{1}
\]

Let \(D=\partial H\), let \(s=|H|-|D|/2\) (so \(s=7\) in the ordinary
two-outlier geometry and \(s=9\) in the triple geometry), and set

\[
 q=|P\setminus D|.
\]

On the \(H\)-cut of \(P\), write \(U\) and \(N\) for the numbers of
positive and negative \(K\)-edges, and put \(I=|H[P]|\).  The lower half of
the complete Max cut box and (1) say

\[
 Z_H(P)=U-N\le-a.
\]

There is therefore a unique integer \(g\ge0\) such that

\[
 Z_H(P)=-a-g,\qquad N=a+g+U.                            \tag{2}
\]

For every vertex define its nonnegative even degree excess

\[
 e(v)=\deg_H(v)-{\bf1}_D(v).
\]

Counting degrees on \(P\) in two ways gives the exact identity

\[
\begin{aligned}
 a-q+\sum_{v\in P}e(v)
   &=\sum_{v\in P}\deg_H(v)\\
   &=(U+N)+2I
     =a+g+2U+2I,
\end{aligned}
\]

and hence

\[
 \boxed{\sum_{v\in P}e(v)=q+g+2(U+I).}                 \tag{3}
\]

Since \(\sum_v e(v)=2|H|-|D|=2s\), every centered-circle \(\mu=1\)
realization must satisfy

\[
 \boxed{q+g+2(U+I)\le2s.}                              \tag{4}
\]

In particular,

\[
 g\equiv q\pmod2,qquad
 U+I\le s-{q+g\over2},qquad q+g\le2s.                \tag{5}
\]

Thus almost the entire circle cut is forced: it has at least \(p-3+g\)
negative crossing edges, while positive crossing edges and edges internal
to \(P\) share a budget of at most \(s-(q+g)/2\).  This is strictly stronger
than the earlier inequality \(q\le2s\), which discarded the signs and the
internal edges.

For completeness, the geometric input \(q\) is explicit.  If
\(L_i=b+u_i\mathbf F_p\) are the two replacement lines, then

\[
 O:=\Gamma\setminus D=(\Gamma\cap L_1)\mathbin\triangle
                         (\Gamma\cap L_2),
 \qquad q=|O\setminus(E\cup\{n\})|.                    \tag{6}
\]

Writing \(\beta_i=b/u_i\), the intersection \(\Gamma\cap L_i\) is governed by

\[
 t^2+\operatorname{Tr}(β_i)t+N(β_i)-{R\over N(u_i)}=0,          \tag{7}
\]

so its size is \(1+\eta(\Delta_i)\), with

\[
 \Delta_i=\operatorname{Tr}(β_i)^2-4N(β_i)
             +{4R\over N(u_i)}.                         \tag{8}
\]

Equations (2)--(8) are the exact remaining algebraic condition supplied by
the centered-circle geometry and the full Max lower cut.  They do not by
themselves give a contradiction: even at \(q=g=0\), the allowed budgets
\(U+I\le7\) and \(U+I\le9\) are nonempty.  A uniform close must now force
\(q+g+2(U+I)>2s\), or obtain an independent lower bound on \(U+I\)
from another Max cut or the boundary phase.  No such lower bound is asserted
here.
