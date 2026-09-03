# Gaussian saturation forces the central two-half saddle

**Status:** proved all-orders necessary structure at the `1/pi` floor.  The
dual-Gaussian proof has an exact expectation-gap identity, and its positive
and negative outputs are Hamming-central up to a square-defect variance.
Consequently a sharp outgoing-half orientation near the floor necessarily
contains the balanced two-half saddle seen by the calculator.  This is not an
orientation construction and does not close the MathOverflow limit.

## 1. Exact saturation gap

Let `H` be a symmetric zero-diagonal `{0,+1,-1}` matrix of order `N` whose
nonzero support is `d`-regular, `d>=2`, and put `P=Phi(H)`.  For a standard
Gaussian vector `g`, define

\[
 Z^\pm=(I\pm H/\sqrt d)g,\qquad X^\pm=\operatorname{sgn}Z^\pm .       \tag{1}
\]

For a support edge `e={i,j}`, put

\[
 u_e={H_{ij}(H^2)_{ij}\over2d},\qquad v=d^{-1/2},\qquad
 \Delta(u,v)=\arcsin(u+v)-\arcsin(u-v),                              \tag{2}
\]

and define the exact Gaussian lower value

\[
 \mathcal L_G(H)={1\over\pi}\sum_{e\in E(H)}\Delta(u_e,v).           \tag{3}
\]

The covariance calculation in Proposition 6.5e gives the equality

\[
 {1\over2}\left(\mathbb E Q_H(X^+)-\mathbb E Q_H(X^-)\right)
 =\mathcal L_G(H).                                                     \tag{4}
\]

Since both terms below are nonnegative, (4) is equivalently the exact
saturation-gap identity

\[
 \boxed{
 \mathbb E\big[(P-Q_H(X^+))+(P+Q_H(X^-))\big]
 =2\big(P-\mathcal L_G(H)\big).}                                      \tag{5}
\]

In particular, for every `eta>0`,

\[
 \Pr\{Q_H(X^+)<P-\eta\ \hbox{or}\ Q_H(X^-)>-P+\eta\}
 \le {2(P-\mathcal L_G(H))\over\eta}.                               \tag{6}
\]

Thus near equality in the arcsine lower bound produces positive and
negative near-extremizers simultaneously; this is stronger than merely
forcing the square correction to be small.

## 2. The two outputs are Hamming-central

Put

\[
 \mathcal O=\sum_{i=1}^N X_i^+X_i^-,\qquad
 d_H(X^+,X^-)={N-\mathcal O\over2}.                                   \tag{7}
\]

For every coordinate `i`, `Z_i^+` and `Z_i^-` have variance two and
covariance zero, because `(H^2)_{ii}=d`.  Hence they are independent and
`E[X_i^+X_i^-]=0`.  In particular,

\[
 \mathbb E\mathcal O=0,\qquad \mathbb E d_H(X^+,X^-)=N/2.             \tag{8}
\]

There is also a quantitative variance bound:

\[
 \boxed{
 \operatorname{Var}\mathcal O
 \le 3N+{2\over d^2}\sum_{i<j}(H^2)_{ij}^2.}                          \tag{9}
\]

To prove it, normalize the two-dimensional Gaussian pair
`V_i=(Z_i^+,Z_i^-)/sqrt(2)`.  Its coordinates are independent standard
Gaussians.  For `i!=j`, the cross-correlation matrix of `V_i,V_j` is

\[
 C_{ij}=\begin{pmatrix}
  \widetilde u+\widetilde v&-\widetilde u\\
  -\widetilde u&\widetilde u-\widetilde v
 \end{pmatrix},\qquad
 \widetilde u={(H^2)_{ij}\over2d},\quad
 \widetilde v={H_{ij}\over\sqrt d}.                                 \tag{10}
\]

The function `sgn(s)sgn(t)` has Gaussian Hermite rank two.  The standard
Gaussian correlation-operator estimate therefore gives

\[
 |\operatorname{Cov}(X_i^+X_i^-,X_j^+X_j^-)|
 \le\|C_{ij}\|_{\rm op}^2
 \le\|C_{ij}\|_F^2
 ={2H_{ij}^2\over d}+{(H^2)_{ij}^2\over d^2}.                         \tag{11}
\]

For completeness, the correlation-operator estimate follows by taking a
singular-value decomposition of `C_ij` and expanding both functions in the
two-dimensional Hermite basis: all components of total degree below two
vanish, so every surviving multiplier is at most
`||C_ij||_op^2`.  Summing (11), using
`sum_(i<j) H_ij^2=Nd/2`, proves (9).  Chebyshev now gives

\[
 \Pr\{|d_H(X^+,X^-)-N/2|>r\}
 \le {3N+(2/d^2)\sum_{i<j}(H^2)_{ij}^2\over4r^2}.                     \tag{12}
\]

## 3. Consequence for an optimal signing at the universal floor

For a complete signing `A` of order `n`, let `M=Phi(A)` and

\[
 L_n={n(n-1)\over\pi}\arcsin{1\over\sqrt{n-1}}.                     \tag{13}
\]

Proposition 6.5e gives \(\mathcal L_G(A)\ge L_n\), while Proposition 6.5g
gives

\[
 \sum_{i\ne j}(A^2)_{ij}^2
 \le8\pi(n-1)(n-2)^{3/2}(M-L_n).                                    \tag{14}
\]

Equations (5) and (9) therefore imply the explicit bounds

\[
 \begin{aligned}
 \mathbb E\big[(M-Q_A(X^+))+(M+Q_A(X^-))\big]
   &\le2(M-L_n),\\
 \operatorname{Var}\mathcal O
   &\le3n+{8\pi(n-2)^{3/2}\over n-1}(M-L_n).
 \end{aligned}                                                        \tag{15}
\]

Consequently, along any sequence of optimal signings for which

\[
 {M\over n^{3/2}}\longrightarrow{1\over\pi},                         \tag{16}
\]

the coupled Gaussian outputs satisfy

\[
 \boxed{
 Q_A(X^+)=M-o_{\Pr}(n^{3/2}),\quad
 Q_A(X^-)=-M+o_{\Pr}(n^{3/2}),\quad
 d_H(X^+,X^-)=n/2+o_{\Pr}(n).}                                       \tag{17}
\]

Thus Hamming-central opposite-sign near-extremizers are forced if the
universal lower floor is approached.  A proof of the outgoing-half target
cannot dispose of its critical layer by showing that this geometry is
absent.

## 4. Consequence for a sharp outgoing-half orientation

Let `A` be optimal, let `R=A circ S`, and form

\[
 K_0=\begin{pmatrix}A&R\\-R&-A\end{pmatrix},\qquad d=2n-2.             \tag{18}
\]

Suppose along a sequence that

\[
 {M\over n^{3/2}}\to{1\over\pi},\qquad
 D_{\to}(A,S)\le {M\over\sqrt2}+o(n^{3/2}).                          \tag{19}
\]

Proposition 6.5e gives `Sigma(A,R)=o(n^4)`.  Moreover

\[
 \sum_{i<j}(K_0^2)_{ij}^2=\Sigma(A,R)+T(R),\qquad
 T(R)=\sum_i(AR-RA)_{ii}^2\le4n(n-1)^2.                              \tag{20}
\]

Since \(\Phi(K_0)=4D_{\to}\) and the base Gaussian value for \(K_0\) is

\[
 {4n(n-1)\over\pi}\arcsin{1\over\sqrt{2n-2}},                       \tag{21}
\]

(19) makes \(\Phi(K_0)-\mathcal L_G(K_0)=o(n^{3/2})\).  Applying (5) and
(9) to \(K_0\) shows that its two coupled rounding outputs \(Y^+,Y^-\) obey

\[
 Q_{K_0}(Y^+)=\Phi(K_0)-o_{\Pr}(n^{3/2}),\qquad
 Q_{K_0}(Y^-)=-\Phi(K_0)+o_{\Pr}(n^{3/2}),                            \tag{22}
\]

and

\[
 d_H(Y^+,Y^-)=n+o_{\Pr}(n).                                          \tag{23}
\]

This is exactly the calculator's balanced two-half saddle.  Put
`u=(Y^++Y^-)/2`, `v=(Y^+-Y^-)/2` and

\[
 O=Q_{K_0}(u),\qquad J=Q_{K_0}(v),\qquad X=u^TK_0v.                   \tag{24}
\]

Then \(u,v\) have complementary supports of sizes \(n+o_{\Pr}(n)\), and the two
corner energies in (22) give

\[
 \boxed{
 O+J=o_{\Pr}(n^{3/2}),\qquad
 X=\Phi(K_0)-o_{\Pr}(n^{3/2}),\qquad
 |O|,|J|\le\Phi(K_0)/2+o_{\Pr}(n^{3/2}).}                             \tag{25}
\]

The last estimate is the stable form of the exact square-calculus statement
`p(1,1)=P`, `p(1,-1)=-P` from (6.14j1): after division by `P`, compactness
and the full square bound leave `O+J=0`, `X=1`, and `|O|=|J|<=1/2` in every
subsequential limit.

Therefore any sharp orientation at the lower floor must make the signed
regular Gaussian lower bound itself asymptotically tight and must realize
the plotted central saddle.  This is a new necessary-structure theorem, not
an upper construction.

## 5. The canonical conference pair is not an obstruction by itself

There is one useful warning.  If `A^2=(n-1)I`, then the two Gaussian vectors
in (1) are independent as vectors.  With

\[
 c={2\over\pi}\arcsin{1\over\sqrt{n-1}},
\]

their sign covariance matrices are `I+cA` and `I-cA`.  Taking the positive
and negative roundings independently, every skew signing `R` satisfies

\[
 \mathbb E(X^{+T}RX^-)^2
 =\operatorname{tr}\big((I+cA)R(I-cA)R^T\big)
 \le(1+c\sqrt{n-1})^2n(n-1)=O(n^2).                                  \tag{26}
\]

Thus this canonical Hamming-central opposite-energy coupling has skew
interaction \(o_{\Pr}(n^{3/2})\) for every orientation.  The remaining
obstruction, if any, lies in uniform control of the whole diffuse layer, not
in this single Gaussian coupling.
