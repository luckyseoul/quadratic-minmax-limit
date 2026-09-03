# Directed outgoing half-cuts: exact random criterion and first-moment barrier

**Status:** proved exact first-moment criterion and proved obstruction to that
specific certificate. The outgoing-half paving problem, multiplier two, and
the original limit remain **OPEN**. In particular, this note does not prove
that random orientations fail; it proves only that their literal
bad-event first moment is exponentially too large.

## 1. Relation to the current doubling reductions

Proposition 6.5c already proves the sparse opposite-diagonal identity

\[
 4D_{\to}(A,S)=
 \max_{x,y}\bigl(|Q_A(x)-Q_A(y)|+|x^TRy|\bigr),
 \qquad R=A\circ S,
\]

and shows that filling the missing corresponding-copy matching costs at most
\(n\). Thus

\[
 D_{\to}(A,S)\le {M\over\sqrt2}+o_{\rm Dini}(n^{3/2}),
 \qquad M=\Phi(A),
\]

would close the multiplier-two ray. The result below does not repeat that
reduction. It gives its exact independent-random-orientation criterion and
shows why the direct first-moment proof cannot reach the sharp constant.

## 2. Exact cut coordinates

For \(T\subseteq[n]\) and \(s\in\{\pm1\}^n\), put

\[
 C_{T,s}=\sum_{u\in T,\ v\notin T}A_{uv}s_us_v,\qquad
 K_{T,s}=\sum_{u\in T,\ v\notin T}R_{uv}s_us_v.
\]

The signed outgoing and incoming halves of this cut are

\[
 F_\to(T,s)={C_{T,s}+K_{T,s}\over2},\qquad
 F_\leftarrow(T,s)={C_{T,s}-K_{T,s}\over2}.
\]

Complementing \(T\) interchanges them, so

\[
 D_\to(A,S)=
 \max_{T,s}{|C_{T,s}|+|K_{T,s}|\over2}.                       \tag{1}
\]

If \(D_T\) changes the signs on \(T\), then, for \(x=s\) and \(y=D_Ts\),

\[
 Q_A(x)-Q_A(y)=2C_{T,s},\qquad x^TRy=2K_{T,s}.                 \tag{2}
\]

In particular \(|C_{T,s}|\le M\). The zero-error outgoing-half target is
therefore exactly

\[
 |K_{T,s}|\le b_{T,s}:=\sqrt2M-|C_{T,s}|
 \quad\hbox{for every }(T,s).                                 \tag{3}
\]

This is a simultaneous one-sided paving of every signed cut by one
tournament. A pointwise balance for a fixed cut or a fixed state does not
imply (3).

## 3. Exact independent-orientation first moment

Choose the entries \(R_{ij}\), \(i<j\), independently and uniformly from
\(\{\pm1\}\), then extend by skew symmetry. For fixed \((T,s)\),

\[
 K_{T,s}\ \buildrel d\over=\ S_m:=\epsilon_1+\cdots+\epsilon_m,
 \qquad m=|T|(n-|T|).                                         \tag{4}
\]

For \(m\ge1\), define the exact lattice tail

\[
 {\rm Tail}_m(t)
 =2^{-m}\sum_{\substack{0\le j\le m\\|2j-m|>t}}{m\choose j}.
\]

Global spin reversal and \(T\leftrightarrow T^c\) generate free
four-element orbits on nontrivial constraints. Hence

\[
 {\cal F}(A):=
 \sum_{[(T,s)]}{\rm Tail}_{|T|(n-|T|)}(b_{T,s})<1              \tag{5}
\]

is sufficient for an orientation satisfying (3): \({\cal F}(A)\) is exactly
the expected number of violated constraint orbits.

Equivalently, the map \(x=s,\ y=D_Ts\) is a bijection before quotienting,
and independent negations of \(x,y\) give the four representatives. Thus
(5) is

\[
 {1\over4}\sum_{\substack{x,y\\0<d<n}}
 {\rm Tail}_{d(n-d)}
 \left(\sqrt2M-\tfrac12|Q_A(x)-Q_A(y)|\right)<1,
 \qquad d=d_H(x,y).                                            \tag{6}
\]

The swap \((x,y)\leftrightarrow(y,x)\) is one additional redundancy that
is not removed in (5)--(6).  Quotienting it would divide \({\cal F}(A)\) by
two, affecting only a constant prefactor and none of the exponential
conclusions below.

Hoeffding gives the more elementary sufficient condition

\[
 \sum_{[(T,s)]}2\exp\left(
 -{(\sqrt2M-|C_{T,s}|)^2\over2|T|(n-|T|)}\right)<1.            \tag{7}
\]

Failure of (7) alone would be inconclusive, because it sums upper bounds.
The next section instead lower-bounds the exact quantity (5).

## 4. Exact first-moment obstruction

Let \(d=\lfloor n/2\rfloor\) and \(h=d(n-d)=\lfloor n^2/4\rfloor\).
There are \(2^n{n\choose d}\) ordered pairs at this Hamming distance.
Since every threshold in (6) is at most \(\sqrt2M\), monotonicity of the
tail gives

\[
 {\cal F}(A)\ge
 {2^n\over4}{n\choose d}\,{\rm Tail}_h(\sqrt2M).                \tag{8}
\]

The needed lower tail is uniform. For fixed \(c,C>0\), if \(m\ge cn^2\)
and \(0\le t\le Cn^{3/2}\), then

\[
 \Pr(|S_m|>t)\ge
 \exp\left(-{t^2\over2m}-O_{c,C}(\log n)\right).                \tag{9}
\]

Indeed, choose the least \(k>t\) with the parity of \(m\). A single binomial
atom and uniform Stirling bounds give

\[
 2^{-m}{m\choose(m+k)/2}
 \ge c_0m^{-1/2}\exp[-mI(k/m)],
\]

where

\[
 I(a)={(1+a)\log(1+a)+(1-a)\log(1-a)\over2}
 ={a^2\over2}+O(a^4).
\]

Here \(k/m=O(n^{-1/2})\), so
\(mI(k/m)=t^2/(2m)+O(1)\), proving (9).

Write \(M=\alpha_n n^{3/2}\), with \(\alpha_n=O(1)\). Stirling for the
central binomial coefficient and (8)--(9) give the exact-first-moment lower
bound

\[
 \boxed{
 {\cal F}(A)\ge
 \exp\bigl((\log4-4\alpha_n^2-o(1))n\bigr).
 }                                                            \tag{10}
\]

For an optimal signing, the known
\(\limsup\alpha_n\le1/2\) makes this at least

\[
 \exp\bigl((\log4-1-o(1))n\bigr),
 \qquad \log4-1=0.386294\ldots>0.                              \tag{11}
\]

Thus even the fictitious most favorable assignment
\(C_{T,s}=0\) on every central constraint would require

\[
 \alpha_n>\sqrt{\log2/2}=0.588705\ldots
\]

for the first moment to decay, already above the optimal-signing upper scale
\(1/2+o(1)\). At the universal lower-floor scale
\(\alpha_n\to1/\pi\), the exponent in (10) is
\(\log4-4/\pi^2=0.981009\ldots\).

The conclusion is deliberately narrow. Equation (11) prevents a proof that
finds an orientation by making the expected number of bad constraints less
than one. It does not show that the probability of at least one bad
constraint tends to one, because these events have extreme overlap.
Dependency-sensitive process bounds, correlated sampling, and
\(A\)-dependent orientations are not excluded.

## 5. Duplication boundary and verification

The two-half note
NOTE_2026-09-02_TWO_HALF_RANDOM_SKEW_CRITERION.md treats Proposition 6.5's
equal-endpoint **sum-energy** diamond. The present result applies the exact
tail calculation to Proposition 6.5c's outgoing-half **difference-energy**
formulation. The entropy mechanism is shared and is cited rather than
presented as a new construction; the new content is the exact outgoing
criterion (1)--(7) with its fourfold quotient and the explicit transfer of
the exact first-moment obstruction to that alternative live lift.

The focused test tests/test_prop65h_directed_halfcut.py verifies (1)--(4),
the fourfold quotient, the exact binomial multiplicities, and the constants
in (10)--(11). It is an identity/arithmetic replay, not a signing census.
