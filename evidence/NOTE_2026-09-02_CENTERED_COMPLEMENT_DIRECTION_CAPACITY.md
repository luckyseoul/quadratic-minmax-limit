# Centered complement pair: exact directional sign capacity

**Date:** 2026-09-02
**Status:** proved local directional barrier; no common-graph lift claimed

Let \(p=4r+1\), let
\(\Gamma=\{z:N(z)=R\}\), with \(\eta(R)=-1\), and use the canonical
centered-circle words

\[
 w(z)=\eta(N(z)-R)\quad(z\notin\Gamma),\qquad
 c|_\Gamma\in\{\pm1\},\qquad y_\pm=w\pm c .
\]

Here \(w\) is zero on \(\Gamma\), \(c\) is zero off \(\Gamma\), and
\(y_\pm\) are the two Boolean \(+p\) eigenvectors.  The calculation below
is invariant under orienting \(c\), which only interchanges \(y_+\) and
\(y_-\).

Fix an affine direction \(U\).  A direct quadratic-character sum on every
line parallel to \(U\) gives the following multiset of numbers of \(+1\)
entries of either \(y_+\) or \(y_-\):

\[
\begin{array}{c|c}
\chi(U)=+1&
 ((p-3)/2)^{(p-1)/4},\
 ((p-1)/2)^{(p-3)/2},\
 ((p+1)/2)^{(p-1)/4},\
 (p-1)^1,\ p^1\\[2mm]
\chi(U)=-1&
 ((p+1)/2)^p .
\end{array}                                                   \tag{1}
\]

For reference, (1) follows by substituting \(z=z_0+tu\) in
\(N(z)-R\).  The discriminant separates external, tangent, and secant
lines.  Resolving each zero with the sparse sign \(c\) gives the five
square-direction values in (1); in a nonsquare direction the two
resolutions balance and every line sum is one.

Consequently the number of edges parallel to \(U\) on which the switched
conference feature

\[
 \chi(U)y_\pm(z)y_\pm(z')
\]

is negative is, for either sign and either direction type,

\[
 A={p(p-1)^2\over4}.                                      \tag{2}
\]

Indeed, in a square direction the negative edges join opposite signs, and
\(\sum_{\ell\parallel U}(\sum_\ell y_\pm)^2=2p^2-p\).
In a nonsquare direction negative features join equal signs, and every
line sum is one.  Both substitutions give (2).

There is also an exact simultaneous table.  In the order
\((\operatorname{sgn}_{y_+},\operatorname{sgn}_{y_-})\), the counts among
all finite edges parallel to one direction are

\[
\begin{array}{c|rrrr}
&(--)&(-+)&(+-)&(++)\\ \hline
\chi(U)=+1&
 {p(p-1)(p-3)\over4}&{p(p-1)\over2}&{p(p-1)\over2}&
 {p(p-1)^2\over4}\\[1mm]
\chi(U)=-1&
 {p(p-1)(p-3)\over4}+1&{p(p-1)\over2}-1&
 {p(p-1)\over2}-1&{p(p-1)^2\over4}+1 .
\end{array}                                                   \tag{3}
\]

To see the mixed columns, \(y_+\) and \(y_-\) differ on exactly one
endpoint precisely for edges crossing \(\Gamma\).  A square direction has
two tangents and \((p-1)/2\) secants, hence \(p(p-1)\) such edges.  A
nonsquare direction has no tangents and \((p+1)/2\) secants, hence
\(p(p-1)-2\).  Equality of the two marginal counts (2) splits these totals
equally.  The diagonal columns then follow from the total
\(p^2(p-1)/2\).

At the lower branch-B endpoint, Proposition 15.758 gives

\[
 P_U\in\{r+2,r+3\}\quad(\chi(U)=+1),\qquad
 P_U=r\quad(\chi(U)=-1).                                  \tag{4}
\]

For every live \(r\ge13\), the \((--)\) entry of (3) is larger than the
corresponding quota in (4).  Hence the sharp lower bound obtainable from
each directional occupancy separately is

\[
 S_{H_U}(y_+)+S_{H_U}(y_-)\ge-2P_U,                       \tag{5}
\]

and (5) is locally attainable by choosing all \(P_U\) edges from the
\((--)\) class.  Summing (5) gives only \(-2|H|\), whereas the paired
full-Max condition requires \(S_H(y_+)+S_H(y_-)\ge8\).

Thus the exact branch-B parallel quotas provide no useful lower bound on
the centered complement-word scores, even simultaneously.  The
off-diagonal edge-Radon equations may still obstruct making all these
local choices in one nonnegative simple graph; (3)--(5) do **not** claim
such a lift.  Any successful directional score bound must therefore use
that common cross-direction lift, not the per-direction quotas or the
one-row character distribution.
