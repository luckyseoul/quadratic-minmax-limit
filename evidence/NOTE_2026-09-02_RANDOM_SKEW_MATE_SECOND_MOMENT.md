# Random skew mates, square defects, and the spectral upper-bridge no-go

**Status:** proved an all-orders approximate-mate construction and proved
that its natural spectral conversion is too weak.  Near the universal
`1/pi` floor, the square/commutator rigidity required by Proposition 6.5e is
automatically satisfiable.  This removes a possible contradiction route but
does not control the outgoing-half cube norm and does not settle the MO
problem.

Let `A` be an order-`n` signing, `n>=3`, and let `R` be a uniformly random
skew signing.  Equivalently, choose a uniform tournament `S` and put
`R=A circ S`.  Define

\[
 P=A^2-R^2,\qquad C=AR-RA,\qquad
 \Sigma(A,R)=\sum_{i\ne j}(P_{ij}^2+C_{ij}^2).              \tag{1}
\]

## 1. Exact expectation

For every ordered pair `i!=j`, expand over `k notin {i,j}`.  The tournament
edge variables are independent and centered, so

\[
 \mathbb E(R^2)_{ij}=0,\qquad
 \mathbb E(R^2)_{ij}^2=n-2.                                \tag{2}
\]

Consequently

\[
 \mathbb E P_{ij}^2=(A^2_{ij})^2+n-2.                      \tag{3}
\]

Likewise

\[
 C_{ij}=\sum_{k\ne i,j}A_{ik}A_{kj}(S_{kj}-S_{ik}),
\]

and the summands are centered and mutually orthogonal in `L^2`, each with
variance two.  Hence

\[
 \mathbb E C_{ij}^2=2(n-2).                                \tag{4}
\]

Summing (3)--(4) proves the exact identity

\[
 \boxed{
 \mathbb E_R\Sigma(A,R)
 =\sum_{i\ne j}(A^2_{ij})^2+3n(n-1)(n-2).
 }                                                           \tag{5}
\]

In particular some skew signing attains the right side or less.  It can be
found deterministically by exposing the tournament edges one at a time and
choosing the sign that does not increase the conditional expectation in
(5).

## 2. The codegree defect is paid for by `M-1/pi`

Apply Proposition 6.5e to the complete support of `A`, with `d=n-1`.  Put

\[
 L_n={n(n-1)\over\pi}\arcsin{1\over\sqrt{n-1}}.             \tag{6}
\]

The square correction there gives

\[
 \boxed{
 \sum_{i\ne j}(A^2_{ij})^2
 \le8\pi(n-1)(n-2)^{3/2}\bigl(\Phi(A)-L_n\bigr).
 }                                                           \tag{7}
\]

Combining (5)--(7), there is a skew signing with

\[
 \Sigma(A,R)\le
 8\pi(n-1)(n-2)^{3/2}(M-L_n)+3n(n-1)(n-2).                 \tag{8}
\]

Writing `M=alpha_n n^(3/2)` gives

\[
 \boxed{
 {\min_R\Sigma(A,R)\over n^4}
 \le8\pi(\alpha_n-1/\pi)+O(1/n).
 }                                                           \tag{9}
\]

Therefore, along any subsequence with `alpha_n->1/pi`, a skew mate with
`Sigma=o(n^4)` exists automatically.  Proposition 6.5e says that such a mate
is necessary for an outgoing-half solution near that floor; (9) shows this
necessary rigidity is consistent and cannot itself yield a contradiction.
There remains a factor-eight gap between the leading sufficient construction
bound here and the `64pi(alpha-1/pi)` necessary allowance in Proposition
6.5e, but neither direction controls the cube norm.

## 3. Including the missing-matching defect

For

\[
 K_0=\begin{pmatrix}A&R\\-R&-A\end{pmatrix},\qquad d=2n-2,
\]

the diagonal entries of `C=AR-RA` occupy the missing perfect matching in
`K_0^2`.  Put

\[
 T(R)=\sum_iC_{ii}^2.                                      \tag{10}
\]

Since `C_ii` is twice a random tournament row sum,

\[
 \mathbb ET(R)=4n(n-1).                                    \tag{11}
\]

Thus

\[
 \mathbb E(\Sigma+T)
 =\sum_{i\ne j}(A^2_{ij})^2+n(n-1)(3n-2).                  \tag{12}
\]

The same conditional-expectation argument produces a deterministic `R` no
worse than (12).

## 4. The strongest generic spectral bridge still loses `pi/2`

Direct multiplication and the diagonal identities give

\[
 K_0^2-dI=
 \begin{pmatrix}P-dI&C\\C&P-dI\end{pmatrix},\qquad
 \|K_0^2-dI\|_F^2=2(\Sigma+T),                             \tag{13}
\]

and `tr(K_0^2-dI)=0`.  For a traceless symmetric matrix of order `2n`, its
largest eigenvalue is at most
`sqrt((2n-1)/(2n))` times its Frobenius norm.  Therefore

\[
 \|K_0\|_{\rm op}^2
 \le d+\sqrt{2n-1\over2n}\sqrt{2(\Sigma+T)}.                \tag{14}
\]

Since `Phi(K_0)<=n||K_0||op` and `Phi(K_0)=4D_to`,

\[
 D_{\to}(A,S)\le{n\over4}
 \sqrt{d+\sqrt{2n-1\over2n}\sqrt{2(\Sigma+T)}}.            \tag{15}
\]

Even the fictitious ideal `Sigma+T=0` yields only

\[
 D_{\to}\le\left({1\over2\sqrt2}+o(1)\right)n^{3/2},       \tag{16}
\]

whereas at the universal lower floor the target is
`n^(3/2)/(pi sqrt(2))`.  The ratio is `pi/2`.

More generally, along a subsequence on which `alpha_n -> alpha` and
`(Sigma+T)/n^2 -> s`, (15) could imply the desired upper bound only if

\[
 2+\sqrt{2s}\le8\alpha^2.                                  \tag{17}
\]

The asymptotic upper bound gives `alpha<=1/2`, so (17) forces
`alpha=1/2` and `s=0`.  Thus this generic trace/Frobenius-to-operator
conversion cannot close the outgoing-half target in the unknown interior
range.  It could become sharp only on the value-specific upper edge.  This
does not exclude a special-structure estimate that uses more than the scalar
defect `Sigma+T`.

## 5. Verdict

Random orientation plus the signed-regular correction constructs the
approximate commuting/equal-square mate demanded near `alpha=1/pi`.
However, converting that second-moment fact to a cube-norm upper bound loses
the fixed factor `pi/2` even at zero defect.  The remaining theorem must use
the statewise energy-layer correlation or nonlinear cover, not merely the
spectral norm of an approximate mate.
