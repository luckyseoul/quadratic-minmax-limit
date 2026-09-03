# Circle/complement Boolean pair and boundary phase

**Date:** 2026-09-02
**Status:** proved paired full-Max constraint; residual (ii) remains open

This note records a second exact consequence of either uniform circle
exception in `NOTE_2026-09-02_CIRCLE_BOUNDARY_GEOMETRY.md`.  It does not
classify one-spike integral eigenvectors and does not use the conditional
affine aliases of Proposition 15.755.

Work in the \(x,\epsilon\) gauge of
`NOTE_2026-09-02_TWO_HALF_NEAR_PENCIL_REDUCTION.md`, so

\[
 K=\epsilon\operatorname{diag}(x)C\operatorname{diag}(x),
 \qquad Kd=pd,
\]

where \(d\) is the signed sparse word supported on the square circle
\(\Gamma\).  Transporting Proposition 15.634's corresponding signed
complement word into the same \(x,\epsilon\) gauge supplies \(w\) with

\[
 Kw=pw,\qquad w|_\Gamma=0,\qquad w_j\in\{\pm1\}
 \quad(j\notin\Gamma).
\]

The two disjoint-support sums

\[
 \boxed{y_+=d+w,\qquad y_-=d-w}                    \tag{1}
\]

are therefore Boolean \(+p\) eigenvectors of \(K\).  This pair is
unconditional for every square circle; it is not an assertion that an
arbitrary one-spike vector belongs to Proposition 15.755's affine family.

Let

\[
 S_H(y)=\sum_{uv\in H}K_{uv}y_uy_v,
 \qquad \sigma_D(y)=\prod_{v\in D}y_v,
 \qquad D=\partial H.
\]

At the lower endpoint \(p=4r+1\), the base \(H\)-score is \(8-3p\),
\(|H|=4r^2+6r+5\), and the product of its \(K\)-edge signs is
\((-1)^r\).  Hence

\[
 \sigma_D(y)=+1\Longrightarrow S_H(y)\equiv1\pmod4,
 \qquad
 \sigma_D(y)=-1\Longrightarrow S_H(y)\equiv3\pmod4.       \tag{2}
\]

Every Boolean \(+p\) eigenvector has \(S_H(y)\ge3\) by the global
\(B\)-maximum bound.  Combining this with (2) sharpens the score floor to

\[
 \boxed{\sigma_D(y)=+1\Longrightarrow S_H(y)\ge5,
 \qquad \sigma_D(y)=-1\Longrightarrow S_H(y)\ge3.}         \tag{3}
\]

If \(T_\pm=\{j:(y_\pm)_j=-1\}\), measured from the all-one base vector,
then

\[
 Z_H(T_\pm)={8-3p-S_H(y_\pm)\over2}.                       \tag{4}
\]

Thus a phase-\(+1\) member of the pair obeys
\(Z_H(T)\le-3(p-1)/2\), while a phase-\(-1\) member obeys
\(Z_H(T)\le-(3p-5)/2\).

Finally, (1) gives \(y_+y_-=1\) on \(\Gamma\) and \(-1\) off \(\Gamma\).
Since every graph boundary has even size,

\[
 \boxed{\sigma_D(y_+)\sigma_D(y_-)
       =(-1)^{|D\cap\Gamma|}.}                            \tag{5}
\]

For a removed pencil line, \(|D\cap\Gamma|=p-2\) is odd.  The two
phases are therefore opposite, so (3)--(4) give the paired constraints

\[
 S_H(y_+)+S_H(y_-)\ge8,
 \qquad Z_H(T_+)+Z_H(T_-)\le-(3p-4).                       \tag{6}
\]

For a centered norm circle, (5) reduces the phase relation to the parity
of its intersections with the two replacement lines.  These are genuine
additional full-Max cuts, but their switching shores contain roughly half
the affine plane because of \(w\).  The odd-boundary surplus estimate is
therefore of order \(p^2\), while (4) asks only for order-\(p\) signed cuts.
Equations (1)--(6) do not alone contradict either exception.  The missing
input is a bound on the two complement-word scores \(S_H(y_\pm)\), or an
equivalent simultaneous cut estimate which uses the directional structure
of \(H\).
