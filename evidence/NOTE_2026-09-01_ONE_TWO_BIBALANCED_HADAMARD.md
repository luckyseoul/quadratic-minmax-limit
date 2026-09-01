# Bi-balanced Hadamard shield for the 1:2 composition

**Status:** proved infinite-family reduction (Proposition 6.8). This does not
close the multiplier-three ray. It replaces the three-state tetrahedral gate
by a two-state rectangular gate and proves that gate outside a fixed
pair-central/high-joint-energy residue.

## 1. Exact rectangular diamond

Let `A` be an order-`n` signing, let `B` be an order-`2n` signing, and let
`C` be an `n` by `2n` sign matrix. Put

\[
 {\cal J}=\begin{pmatrix}A&C\\ C^T&B\end{pmatrix}.
\]

For Boolean `x,y`, flipping all signs of `x` fixes both quadratic terms and
reverses the rectangular term. Therefore

\[
 \Phi({\cal J})=
 \max_{x,y}\left(|Q_A(x)+Q_B(y)|+|x^TCy|\right).                 \tag{1}
\]

Take `A,B` optimal and write

\[
 M=m_n,\qquad N=m_{2n},\qquad
 T=(M^{2/3}+N^{2/3})^{3/2}.                                    \tag{2}
\]

Thus `Phi(J)<=T+e_n` gives

\[
 H(3n)\le H(n)+H(2n)+O(e_n/n^{1/2}).                            \tag{3}
\]

In particular, `e_n=n^(3/2) omega(n)` gives an `O(n omega(n))` error in
`H`. This is the exact 1:2 target.

## 2. Uniform headroom

For positive `a,b`, set

\[
 f(a,b)=(a^{2/3}+b^{2/3})^{3/2}-a-b.
\]

Both partial derivatives are positive. Proposition 5.2 gives

\[
 M\ge {n\sqrt{n-1}\over\pi},\qquad
 N\ge {2n\sqrt{2n-1}\over\pi}.
\]

Consequently

\[
 {T-M-N\over n^{3/2}}\ge d_0-o(1),\qquad
 d_0={3\sqrt3-1-2\sqrt2\over\pi}
     =0.4353604839\ldots .                                      \tag{4}
\]

This is uniform over the actual values of `M,N`; no doubling estimate is
used in (4).

## 3. Four-anchor pairings

Choose Boolean states

\[
 Q_A(z_+)=P(A),\quad Q_A(z_-)=-N(A),\qquad
 Q_B(w_+)=P(B),\quad Q_B(w_-)=-N(B).
\]

Gauge by `z_+` and `w_+`, and put

\[
 r_A=z_+\circ z_-,\qquad r_B=w_+\circ w_-.
\]

Partition the `n` row coordinates into the two level sets of `r_A`, and
pair coordinates arbitrarily inside each level set. Similarly partition
the `2n` column coordinates into the two level sets of `r_B` and pair inside
each level set. There are

\[
 r\le n/2\quad\hbox{row pairs},\qquad
 s\le n\quad\hbox{column pairs},                                \tag{5}
\]

and at most two unpaired coordinates on either side.

For a gauged Boolean row state `xi=z_+ circ x`, define

\[
 u_a={\xi_i-\xi_{i'}\over2}\in\{0,\mathord\pm1\}
\]

on row pair `a={i,i'}`, and define `v_b` analogously from
`eta=w_+ circ y` on each column pair. Put

\[
 k_A(x)=|\operatorname{supp}u|,\qquad
 k_B(y)=|\operatorname{supp}v|.                                 \tag{6}
\]

Every one of `+-z_+,+-z_-` has `k_A=0`, and every one of
`+-w_+,+-w_-` has `k_B=0`. Distance `h` from the corresponding projective
anchor implies `k_A<=h` or `k_B<=h`.

## 4. Hadamard core and exact product shield

The effective prime-number theorem in the progression `3 mod 4` supplies,
for every sufficiently large `n`, a prime `p` with

\[
 n\le p\le n+n\exp(-c\sqrt{\log n}),\qquad p=3\pmod4.
\]

Let `q=p+1`. A Paley Hadamard matrix of order `q` exists. Take any
`r` by `s` submatrix `E`; then

\[
 \|E\|_{\rm op}\le\sqrt q,
 \qquad {q\over n}=1+O(\exp(-c\sqrt{\log n})+n^{-1}).            \tag{7}
\]

On the paired core, put the `2` by `2` tile

\[
 E_{ab}\begin{pmatrix}1&-1\\-1&1\end{pmatrix}                  \tag{8}
\]

on row pair `a` and column pair `b`. Fill every entry incident to an
unpaired row or column arbitrarily by signs, and undo the two gauges. This
produces an `n` by `2n` sign matrix `C`.

The core contribution has the exact identity

\[
 x^TCy=4u^TEv+\operatorname{border}.                             \tag{9}
\]

There are at most `2(2n)+2n=6n` border entries. Hence, for every Boolean
pair,

\[
 \boxed{|x^TCy|\le4\sqrt{q\,k_A(x)k_B(y)}+6n.}                 \tag{10}
\]

This is deterministic. In particular, the cross term is `O(n)` whenever
either state is one of the selected positive or negative extrema.

## 5. Closed region and exact complement

Combining (1), `|Q_A+Q_B|<=M+N`, (4), and (10) proves, for every
sufficiently large `n`,

\[
 k_A(x)k_B(y)\le {n^2\over100}
 \quad\Longrightarrow\quad
 |Q_A(x)+Q_B(y)|+|x^TCy|\le T,                                 \tag{11}
\]

because `4 sqrt(1/100)=0.4<d_0` and the errors in (7), together with
`6n/n^(3/2)`, tend to zero with a dyadic-Dini envelope.

The sharper statewise shield is

\[
 4\sqrt{q\,k_A(x)k_B(y)}+6n
 \le T-|Q_A(x)+Q_B(y)|.                                        \tag{12}
\]

Thus this construction leaves only pairs satisfying both

\[
 k_A(x)k_B(y)>{n^2\over100},                                   \tag{13}
\]

and the strict reverse of (12). In particular it closes fixed linear
Hamming strips around all positive and negative extremizers on either side:
using `k_B<=n`, a row-side distance at most `n/100` suffices; using
`k_A<=n/2`, a column-side distance at most `n/50` suffices.

## 6. Why this is not Proposition 6.7 in disguise

If `B` is itself the equal-endpoint doubling frame

\[
 B=\begin{pmatrix}A&P+D\\-P+D&A\end{pmatrix}
\]

and the small block is `-A`, then grouping the three layers as `1+2` turns
(1) into Proposition 6.7's `K_3`, up to the three diagonal matching terms.
The two branches in `K_3` are the two possible internal layer flips inside
`B`. Thus nested doubling does not simplify tripling.

The construction above escapes that equivalence because `B` is an
independently optimal order-`2n` signing and `C` is built from its actual
positive and negative extremizers. The remaining gate (13) plus the reverse
of (12) is genuinely a two-state gate.

## 7. Multi-anchor extension

For any fixed finite lists of row and column anchors, refine the coordinate
partition by their complete relative-sign signatures, pair inside every
signature class, and repeat (8). Every listed anchor then has zero paired
difference. With `R` row anchors and `S` column anchors there are at most
`L_A<=2^(R-1)` and `L_B<=2^(S-1)` leftovers, respectively. The spectral
core in (10) is unchanged, while its literal `6n` border is replaced by at
most

\[
 (2L_A+L_B)n\le(2^R+2^{S-1})n.
\]

For fixed `R,S` this is still `O(n)`, so the `1/100` shield survives after
increasing the finite threshold. This is the exact halving/refinement
mechanism available for the residual states; it is not a finite-state
census.
