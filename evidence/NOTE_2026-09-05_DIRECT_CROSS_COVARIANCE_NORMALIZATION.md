# Direct cross normalization on actual conditional optimizers

2026-09-05. Analytic proof draft for independent complete review.
The result is a Gaussian reduction and a conditional-optimizer floor,
not an evaluated order upper bound.

## 1. Source and directly normalized covariance

Let \(n\ge2\), let \(A\) be a complete symmetric zero-diagonal signing
of order \(n\), and let \(B\) be any \(n\) by \(n\) sign matrix. Set
\[
 K=\begin{pmatrix}A&B\\B^T&-A\end{pmatrix},\qquad
 L_K=\|K\|_{\rm op},\qquad N=2n.
\]
On the \(n^2\) cross coordinates define the self-adjoint operator
\[
 {\cal S}_B(X)=BX^TB,\qquad
 H=A\otimes A-{\cal S}_B+I,\qquad
 \mu=\max(2,\|H\|_{\rm op}),\qquad \ell=\mu+1,
                                                        \tag{1}
\]
where the tensor has entries \(A_{ik}A_{jl}\) on coordinates \(ij,kl\).
The cross correlation matrix is
\[
 \boxed{R_\mu=I+H/\mu
       =\{\ell I+A\otimes A-{\cal S}_B\}/\mu.}             \tag{2}
\]
This normalization need not equal \(L_K^2-1\).

The diagonal of \(H\) is zero, since \(A_{ii}=0\) and \(B_{ij}^2=1\).
For a fixed row \(i\), its principal submatrix on \((i,j)\), \(1\le j\le n\),
is \(I-b_i b_i^T\), where \(b_i=(B_{i1},\ldots,B_{in})^T\).
It has the eigenvalue \(1-n\). Consequently
\[
 \mu\ge n-1,\quad \mu\ge2,\qquad
 \operatorname{diag}R_\mu=1,\qquad 0\preceq R_\mu\preceq2I.
                                                        \tag{3}
\]
Thus (2) is an admissible latent Gaussian covariance at every order.
For \(n\ge3\), the maximum with 2 in (1) is redundant.

Let \(G_\mu\) have covariance \(R_\mu\) on the actual matrix entries,
and let \(W\) have independent standard normal entries, independently.
For any deterministic real \(h\), set
\[
 B_{h,ij}=\operatorname{sign}(G_{\mu,ij}+hB_{ij}),\qquad
 s=2\Phi_{\rm Gauss}(h)-1,\quad k=4\phi(h)^2,\quad
 v=1-s^2-k\ge0,
\]
\[
 Z_h=sB+\sqrt{k}\,G_\mu+\sqrt v\,W.                        \tag{4}
\]
All cross entries, including those with \(i=j\), are treated identically.
For arbitrary fixed real \(I(x,y)\), and \(|\theta|\le1\), write
\[
 {\cal M}_I(D)=\max_{x,y\in\{-1,1\}^n}
                    |I(x,y)+\theta x^TDy|.
\]
There are absolute constants \(C_1,C_2\) such that
\[
 \boxed{
 |\mathbb E{\cal M}_I(B_h)-\mathbb E{\cal M}_I(Z_h)|
 \le C_1 n^{16/11}+C_2 n\sqrt{1+L_K}.}                    \tag{5}
\]
The constants are uniform over \(A,B,h,I,\theta\).

## 2. Hermite remainder after direct cross normalization

We reuse the complete Hermite and four-cycle proofs of
NOTE_2026-09-05_WHOLE_EDGE_SOURCE_PRESERVING_GAUSSIAN_REDUCTION.md,
but check their new normalization explicitly. No Gaussian covariance
on the unused unordered edges is assumed here.

For distinct cross coordinates, (2) gives
\[
 (R_\mu)_{ij,kl}=
 {A_{ik}A_{jl}-B_{il}B_{kj}\over\mu}.                      \tag{6}
\]
Coordinates sharing a row or column have correlations of magnitude
\(1/\mu\); all other pairs have correlations \(0,\pm2/\mu\).
Every argument below therefore belongs to \([-1,1]\), by (3).

For the centered threshold activation, let \(c_j(h)\) be its normalized
Hermite coefficients. Then \(c_1^2=k\), \(\sum_{j\ge2}c_j^2=v\), and
threshold orientation by \(B_{ij}\) changes the sign of exactly the
even coefficients. Write
\[
 o(t)=\sum_{\substack{j\ge3\\j\ {\rm odd}}}c_j^2t^j,
 \qquad e(t)=\sum_{\substack{j\ge2\\j\ {\rm even}}}c_j^2t^j,
 \qquad |o(t)|\le |t|^3,\quad 0\le e(t)\le t^2.
\]
Let \(b=\operatorname{vec}B\), let \(D_b=\operatorname{diag}b\),
and put \(u=e(1/\mu)\), \(w=e(2/\mu)/2\). Let \(C_h\) be the exact
centered covariance of \(B_h\), and \(C_0=kR_\mu+vI\).

For clarity one can index auxiliary matrices by ALL unordered edges
of the complete order-\(N\) signing \(K\). Let \(E_N\) be the line-graph
adjacency matrix, and let \(Q_N\) be the four-cycle matrix with entries
\[
 (Q_N)_{ij,kl}=K_{ik}K_{jk}K_{il}K_{jl}.
\]
Let \(\Pi\) select only cross edges. The purely algebraic bounds from
the cited complete proof are
\[
 \|E_N\|=2(N-2),\qquad \|Q_N\|\le L_K(N-1)/2.              \tag{7}
\]
In particular, the second bound follows by squaring the full ordered-pair
four-cycle matrix and using
\((K^2)\circ(K^2)\preceq L_K^2(N-1)I\). It does not depend on the
positivity of any auxiliary covariance on the unused edges.

The ENTIRE even-Hermite remainder on the cross coordinates is exactly
\[
 R_{\rm even}
 =wbb^T+
 D_b\Pi\{wQ_N+(u-w)E_N-wI\}\Pi^TD_b.                     \tag{8}
\]
Indeed, on disjoint edges the two numerator products in the compressed
whole-edge minus law have equal sign exactly when \(Q_N=1\); their
correlation magnitude is then \(2/\mu\), and otherwise zero.
Adjacent edges contribute \(u\), and the diagonal is zero. These
entry identities restrict to (6), regardless of the unused edges.
The covariance's diagonal residual Hermite mass is already in \(vI\).

The odd-Hermite remainder has operator norm bounded by the corresponding
full-edge row count:
\[
 \|R_{\rm odd}\|
 \le {2(N-2)(2N-5)\over\mu^3}.                            \tag{9}
\]
This uses only \(|o(t)|\le|t|^3\), not positivity of the formal full-edge
matrix. Combining (7)-(9),
\[
 C_h=(C_0+wbb^T)+{\cal E},\qquad
 \|{\cal E}\|\le
 {2(N-2)(2N-5)\over\mu^3}
 +{L_K(N-1)+6(N-2)+2\over\mu^2}
 \le {8L_K+180\over n}.                                  \tag{10}
\]
For the last loose absolute estimate use \(N=2n\), \(\mu\ge n/2\),
and \(n\ge2\). Both \(C_h\) and \(C_0+wbb^T\) are PSD; the error
need not have a sign.

The retained rank-one Gaussian is \(\sqrt w\,\xi B\), with a scalar
standard normal \(\xi\). Its actual bilinear-norm cost is at most
\[
 \sqrt w\,\mathbb E|\xi|\max_{x,y}|x^TBy|
 \le {2n^2\over\sqrt\pi\,\mu}=O(n).                       \tag{11}
\]
The Gaussian finite-maximum comparison under covariance operator error
\(\delta\), with \(m\) coordinates and \(J\) states of coefficient norm
at most \(\sqrt m\), costs \(\sqrt{2\delta m\log J}\).
It follows by Gaussian convex order and adding independent \(\delta I\)
noise. Here \(m=n^2\), \(J\le2^{2n+1}\), and (10) applies. Thus the
matched Gaussian with mean \(sB\) and covariance \(C_h\) differs from
(4) in expected \({\cal M}_I\) by at most
\[
 C n\sqrt{1+L_K}.                                       \tag{12}
\]
This accounts for the rank-one term through (11), not through its rank.

Finally the independently proved mean-preserving shifted-sign theorem
NOTE_2026-09-05_SHIFTED_SIGN_GAUSSIAN_UNIVERSALITY.md applies directly
to the \(n^2\) cross coordinates, with latent covariance constant two
by (3). Its observable family is the augmented set
\(\{\sigma\theta(x_i y_j)_{ij}\}\), with the fixed internal energy
absorbed into the prior. Its growing-temperature bound gives
absolute error \(C_1n^{16/11}\). This proves (5).

## 3. An elementary operator bound from the ORIGINAL norm

For any real matrix \(T\) with entries bounded in absolute value by one,
let
\[
 \beta_{\mathbb R}(T)=
 \max_{x,y\ {\rm real},\ \|x\|_\infty,\|y\|_\infty\le1}|x^TTy|.
\]
Then
\[
 \|T\|_{\rm op}^2\le2\beta_{\mathbb R}(T).                 \tag{13}
\]
To verify the real/complex constant, rotate a complex pairing \(u^*Tv\)
to the real axis. If \(u=a+ib\), \(v=c+id\) have complex coordinate
moduli at most one, its real part is the sum of the two real-cube forms
\(a^TTc+b^TTd\). Therefore the complex \(\ell^\infty\)-to-\(\ell^1\)
operator norm is at most \(2\beta_{\mathbb R}(T)\). The complex
\(\ell^1\)-to-\(\ell^\infty\) norm is at most one. Riesz--Thorin
interpolation halfway between these endpoints gives (13).

Equivalently, for real unit \(\ell^2\) vectors \(x,y\), omit zero
coordinates and apply the three-lines theorem to
\[
 f(z)=\sum_{ij}T_{ij}\operatorname{sign}(x_i)\operatorname{sign}(y_j)
                  |x_i|^{2z}|y_j|^{2z}.
\]
On the line \(\Re z=0\) its modulus is at most \(2\beta_{\mathbb R}(T)\);
on \(\Re z=1\) it is at most one. At \(z=1/2\) it equals \(x^TTy\).
This also proves the displayed finite-dimensional estimate directly.
For the classical interpolation theorem, see
[Tao, Theorem 21](https://terrytao.wordpress.com/2009/03/30/245c-notes-1-interpolation-of-lp-spaces/).

If \(T\) is symmetric with zero diagonal, its original quadratic norm
bounds the form throughout the real cube by independent sign rounding.
Polarization with \((x+y)/2\) and \((x-y)/2\) then gives
\(\beta_{\mathbb R}(T)\le4\Phi(T)\). Hence
\[
 \boxed{L_K^2\le8\Phi(K).}                                \tag{14}
\]
In particular \(\Phi(K)\le C_0 n^{3/2}\) implies
\[
 n\sqrt{1+L_K}=O_{C_0}(n^{11/8})=o(n^{16/11}).             \tag{15}
\]
This is not a conference-scale bound \(L_K=O(\sqrt n)\).

## 4. Actual conditional minima and the one-sided consequence

For a fixed \(A\), let
\[
 F_A^*=\min_{B\in\{-1,1\}^{n\times n}}
              \Phi\begin{pmatrix}A&B\\B^T&-A\end{pmatrix}.
\]
For every \(B\) the displayed original norm is exactly
\[
 \max_{x,y}\{|Q_A(x)-Q_A(y)|+|x^TBy|\}.                   \tag{16}
\]
Indeed, replacing \(x\) by \(-x\) preserves the internal difference
and changes the sign of the cross energy.

The independent-sign exponential-maximum estimate gives
\[
 F_A^*\le2\Phi(A)+2\sqrt{\log2}\,n^{3/2}.                 \tag{17}
\]
For completeness, with random independent cross signs each fixed
bilinear energy is a sum of \(n^2\) independent signs. The set of
distinct rank-one sign coefficient matrices has cardinality
\(2^{2n-1}\) and already contains both energy signs. Its expected
maximum is at most
\(\sqrt{2n^2(2n-1)\log2}\le2\sqrt{\log2}\,n^{3/2}\).
The minimum is at most the corresponding expectation, proving (17).

Choose ANY exact conditional minimizer \(B_*\), and use (1)-(4) with it.
If \(\Phi(A)\le C_A n^{3/2}\), then (17) supplies the required original
norm cap for \(K\), and (5), (14)-(15) give the uniform comparison
\[
 |\mathbb E{\cal M}_{Q_A(x)-Q_A(y)}(B_h)
       -\mathbb E{\cal M}_{Q_A(x)-Q_A(y)}(Z_h)|
                  \le C_{C_A} n^{16/11}.                 \tag{18}
\]
Every \(B_h\) is an admissible cross signing with the SAME internal
blocks. Conditional original norm optimality therefore proves
\[
 \boxed{F_A^*\le
   \mathbb E\Phi\begin{pmatrix}A&Z_h\\Z_h^T&-A\end{pmatrix}
          +C_{C_A}n^{16/11},\qquad h\in\mathbb R.}         \tag{19}
\]
For every exact order-\(n\) original minimizer \(A\), the elementary
cap \(\Phi(A)\le\sqrt{\log2}\,n^{3/2}\) makes this constant absolute.
Taking the infimum over deterministic thresholds in (19) is valid,
and \(m_{2n}\le F_A^*\). No minimizing threshold need exist.

Equation (19) is an actual conditional-optimizer floor for a new,
directly normalized covariance. It does not evaluate that Gaussian
right side against \(2\sqrt2\,\Phi(A)\), does not assert positivity
of uncorrected joint comparison fields at \(\ell<L_K^2\), and does
not prove original all-orders convergence.
