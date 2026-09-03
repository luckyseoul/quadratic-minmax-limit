# Weighted multi-anchor rounding by a Gaussian `L^1` body

**Status:** proved all-orders integral orientation lemma and a conditional
advance of the directed-half-cut gate.  This does not close multiplier two:
it controls every constraint incident with a prescribed anchor family, but
does not prove that the remaining diffuse constraint graph has a small vertex
cover.  No finite-order census is used.

Let `A` be an order-`n` signing, put `M=Phi(A)`, and let
`x^(1),...,x^(k)` be Boolean anchor states.  The signature-cell construction
in Proposition 6.5f controls an arbitrary `k`-anchor family through the number
of realized coordinate signatures, which in the worst case limits the
critical application to `k=O(log n)`.  The following vector-balancing
construction instead controls an arbitrary weighted family of linear size.

## 1. Weighted integral orientation lemma

Choose positive error probabilities `q_1,...,q_k` with

\[
 \sum_{a=1}^k q_a\le {1\over2},
 \qquad
 \rho_a=n\sqrt{2/\pi}+\sqrt{2n\log(1/q_a)}.                 \tag{1}
\]

Let `B_1,...,B_k>0`.  If

\[
 \boxed{
  50\sum_{a=1}^k {\rho_a^2\over B_a^2}\le1,
 }                                                           \tag{2}
\]

then there is a skew signing `R` such that

\[
 \boxed{
  \|R x^{(a)}\|_1\le B_a\qquad(1\le a\le k).
 }                                                           \tag{3}
\]

In particular, this is an exact integral orientation: every upper-triangular
entry of `R` is `+1` or `-1`.  There is no fractional remainder or parity
loss.

### Proof

For every edge `e={i,j}`, `i<j`, introduce a vector with `k` blocks in
`R^(nk)`.  Given positive weights `lambda_a`, its block `a` is

\[
 v_{ij}^{(a)}=\lambda_a
       \bigl(x_j^{(a)}e_i-x_i^{(a)}e_j\bigr).                \tag{4}
\]

Thus

\[
 \|v_{ij}\|_2^2=2\sum_a\lambda_a^2=:L^2.                   \tag{5}
\]

For a standard Gaussian `g` in `R^n`, the map `g -> ||g||_1` is
`sqrt(n)`-Lipschitz and has mean `n sqrt(2/pi)`.  Gaussian concentration and
(1) give

\[
 \Pr\{\|g\|_1>\rho_a\}\le q_a.                             \tag{6}
\]

Consequently the symmetric convex body

\[
 {\cal K}=\{(z^{(1)},\ldots,z^{(k)}):
                 \|z^{(a)}\|_1\le\rho_a\text{ for every }a\}
                                                                    \tag{7}
\]

has standard Gaussian measure at least `1/2`, by the union bound.
Banaszczyk's `5K` vector-balancing theorem, applied to the Euclidean-unit
vectors `v_ij/L`, supplies signs `r_ij in {+1,-1}` for which

\[
 \sum_{i<j}r_{ij}{v_{ij}\over L}\in5{\cal K}.              \tag{8}
\]

Define `R_ij=r_ij` for `i<j`, `R_ji=-r_ij`, and `R_ii=0`.  The `i`th
coordinate of block `a` in the unnormalized sum in (8) is exactly
`lambda_a (R x^(a))_i`.  Hence

\[
 \|R x^{(a)}\|_1\le {5L\rho_a\over\lambda_a}.              \tag{9}
\]

Take `lambda_a=rho_a/B_a`.  Then

\[
 L=\left(2\sum_a{\rho_a^2\over B_a^2}\right)^{1/2},
\]

and (2) makes the multiplier `5L` in (9) at most one.  This proves
(3).

The only deep external input is the following standard normalization of
Banaszczyk's theorem: Euclidean-unit vectors and a symmetric convex body of
Gaussian measure at least `1/2` admit a signed sum in `5K`.  Equivalently,
vectors of norm at most `1/5` admit a signed sum in `K`.

## 2. Exact insertion into the outgoing-half target

For any anchor `x`, define its uniform incident-pair budget

\[
 B_A(x)=(2\sqrt2-1)M-|Q_A(x)|.                              \tag{10}
\]

This is positive because `|Q_A(x)|<=M`.  If (3) holds with
`B_a=B_A(x^(a))`, then every pair incident with an anchor obeys the exact
opposite-diagonal target:

\[
\begin{aligned}
 |Q_A(x^{(a)})-Q_A(y)|+|(x^{(a)})^TRy|
 &\le |Q_A(x^{(a)})|+M+\|Rx^{(a)}\|_1\\
 &\le2\sqrt2 M                                             \tag{11}
\end{aligned}
\]

for every Boolean `y`.  By Proposition 6.5c, (11) is precisely the required
bound on all signed outgoing-half constraints incident with those anchors.
The smallest value of (10), attained at an extremizer, is

\[
 B_{\min}=2(\sqrt2-1)M.                                    \tag{12}
\]

Thus a simple uniform version follows by taking `q_a=1/(2k)`:

\[
 \boxed{
 \|Rx^{(a)}\|_1\le
 10n\sqrt{k/\pi}+10\sqrt{kn\log(2k)}
 \quad\hbox{for all }a.
 }                                                           \tag{13}
\]

Indeed the right side is `5 sqrt(2k) rho` with
`rho=n sqrt(2/pi)+sqrt(2n log(2k))`.  Therefore every anchor-incident pair
is shielded whenever

\[
 10n\sqrt{k/\pi}+10\sqrt{kn\log(2k)}
 \le2(\sqrt2-1)M.                                          \tag{14}
\]

Using the all-orders arcsine floor

\[
 M\ge {n(n-1)\over\pi}\arcsin{1\over\sqrt{n-1}},          \tag{15}
\]

(14) holds, for all sufficiently large `n`, for every

\[
 k\le(c_*-\delta)n,
 \qquad
 c_*={3-2\sqrt2\over25\pi}
     =0.0021845336957706\ldots,                             \tag{16}
\]

where `0<delta<c_*` is arbitrary.  More generally, if `k=o(n)`, the loss in
(13) is `o(n^(3/2))`; if `k=O(n^(1-delta))`, that normalized loss has a
dyadically summable polynomial decay.

## 3. What changed, and what remains open

This is not the already recorded fractional minimax: (8) lands exactly at a
vertex of the edge-sign cube.  It is also not independent random rounding;
the `5K` theorem selects one globally correlated sign vector.  Compared with
Proposition 6.5f, the worst-case anchor capacity improves from logarithmic in
`n` to a positive linear fraction, and (2) additionally spends less capacity
on nonextremal anchors through their larger budgets (10).

The result still does not orient every constraint.  A completion would need
to prove that all potentially bad pairs have a vertex cover satisfying the
weighted capacity (2), or incorporate the non-anchor pairs into a different
large-Gaussian-measure convex body without losing the critical constant.
Precisely, for the orientation furnished above put

\[
 \Gamma(R)=\max_{x,y}|x^TRy|,
 \qquad
 E_R=\bigl\{\{x,y\}:|Q_A(x)-Q_A(y)|
                 >2\sqrt2M-\Gamma(R)\bigr\}.               \tag{17}
\]

Define the tournament `S=A circ R`.  If the chosen anchors form a vertex
cover of `E_R`, then (11) handles its
incident pairs, while every nonincident pair is safe by the definition of
`E_R`; hence Proposition 6.5c gives `D_to(A,S)<=M/sqrt(2)`.  Proving that a
capacity-(2) cover can be chosen for the *same* resulting `R` is the still-open
cover condition.  The theorem above does not bound `Gamma(R)`, so this
condition cannot be silently inferred from the linear-size anchor capacity.
Accordingly multiplier two and the original MathOverflow limit remain open.

## Reference

W. Banaszczyk, *Balancing vectors and Gaussian measures of n-dimensional
convex bodies*, Random Structures & Algorithms **12** (1998), 351--360,
doi:10.1002/(SICI)1098-2418(199807)12:4<351::AID-RSA3>3.0.CO;2-S.

For the exact unit-vector/`5K` normalization used above, see Theorem 1.1 of
N. Bansal, D. Dadush, S. Garg, and S. Lovett, *The Gram--Schmidt Walk: A
Cure for the Banaszczyk Blues*, Theory of Computing **15** (2019), Article
21, 1--27, doi:10.4086/toc.2019.v015a021.  That paper states explicitly that
the absolute constant in this normalization may be taken to be `5`.
