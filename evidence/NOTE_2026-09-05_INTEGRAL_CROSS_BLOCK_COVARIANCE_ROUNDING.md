# Integral cross-block covariance rounding and an actual-energy saving

## Scope

Let \(U,V\succeq0\) be \(n\times n\) matrices with diagonal one, and define
\[
\overline q(B)=\frac12\left[
 \operatorname{tr}(B^TUBV)+\operatorname{tr}(B^TVBU)\right],
 \qquad B\in\{-1,1\}^{n\times n}.                               \tag{1}
\]
All \(n^2\) entries of \(B\), including its diagonal, are free signs.
The bounds below optimize over all such integral cross blocks. They do
not identify \(B\) with the original signing \(A\).

This note proves a general covariance-rounding bound, its diagonal- and
scalar-shift refinements, ridge and spectral-tail consequences, an
actual-energy-dependent strict improvement over the independent-sign
value \(n^2\), and a sufficient low-effective-rank case of the local
endpoint comparison. The general covariance and spectral-tail bounds
use a precisely cited external theorem; the stronger complete-sign
Gaussian rounding bound is elementary and proved here.

No bound here proves, for actual global half-product minima in general,
\[
\min_B\overline q(B)\le2a_A'(\beta)/\beta+o(n^2).                \tag{2}
\]
Even (2) would only control the indicated endpoint derivative, not the
entire balanced interpolation or the original MO limit. The high-rank
algebraic examples in the last section are explicitly conditional
covariance models, not claimed actual Gibbs or global-minimizer families.

## 1. Imported covariance theorem and the general integral bound

The external input is the generalized covariance theorem of Harshaw,
Sävje, Spielman, and Zhang, *Balancing Covariates in Randomized Experiments
with the Gram–Schmidt Walk Design*, arXiv:1911.03071v8:
[Supplement S3.4, Theorem 6.3*, PDF p.41](https://arxiv.org/pdf/1911.03071v8#page=41).
Its hypotheses allow arbitrary input columns of Euclidean norm at most
one. With zero initial assignment, the output sign vector \(\sigma\)
has mean zero and satisfies
\[
\operatorname{Cov}(W\sigma)\preceq P_{\operatorname{range}W}.
                                                                    \tag{3}
\]
The mean-zero assertion is the zero-start consequence of Lemma 6.1 and
Lemma S3.3 in Supplement S3.3--S3.4, PDF pp.38--39. The generalized
statement, rather than a covariate-specific specialization, is used here.

Put \(m=n^2\), and use column vectorization:
\[
\mathcal K=\frac12(V\otimes U+U\otimes V),\qquad
\overline q(B)=\operatorname{vec}(B)^T\mathcal K\operatorname{vec}(B).
                                                                    \tag{4}
\]
Then \(\mathcal K\succeq0\), its diagonal is one, and
\(\operatorname{tr}\mathcal K=m\).

For any \(G\succ0\) with \(G_{ii}\le1\), the columns of \(W=G^{1/2}\)
have norm at most one. Since \(W\) is invertible, (3) gives
\[
\mathbb E[\sigma\sigma^T]\preceq G^{-1}.
\]
Taking the trace against \(\mathcal K\succeq0\) and then selecting one
sign outcome proves
\[
\boxed{\quad
\min_{B\in\{-1,1\}^{n\times n}}\overline q(B)
 \le\operatorname{tr}(\mathcal K G^{-1})
 \quad(G\succ0,\ \operatorname{diag}G\le1).
\quad}                                                        \tag{5}
\]
This is a bound for an actual integral \(B\), not only a fractional or
semidefinite relaxation. Independent signs also give
\(\mathbb E\overline q(B)=n^2\).

### 1.1. Retaining the fixed coordinate squares

The output vector has \(\sigma_i^2=1\) exactly. Therefore any real
diagonal matrix \(\Lambda\) with \(\mathcal K-\Lambda\succeq0\)
can be separated before applying the covariance domination:
\[
\mathbb E[\sigma^T\mathcal K\sigma]
=\operatorname{tr}\Lambda+
 \operatorname{tr}\bigl[(\mathcal K-\Lambda)
                       \mathbb E\sigma\sigma^T\bigr]
\le\operatorname{tr}\Lambda+
 \operatorname{tr}\bigl[(\mathcal K-\Lambda)G^{-1}\bigr].
\]
Selecting a sign outcome gives the stronger admissible family
\[
\boxed{\quad
\min_B\overline q(B)\le
 \operatorname{tr}\Lambda+
 \operatorname{tr}\bigl[(\mathcal K-\Lambda)G^{-1}\bigr],
\quad
G\succ0,\ \operatorname{diag}G\le1,\quad
\mathcal K-\Lambda\succeq0,\ \Lambda\text{ diagonal}.
\quad}                                                        \tag{5a}
\]
In particular, with \(\Lambda=\ell I\), this becomes
\[
\min_B\overline q(B)\le
 n^2\ell+\operatorname{tr}\bigl[(\mathcal K-\ell I)G^{-1}\bigr],
\qquad \ell\le\lambda_{\min}(\mathcal K).                      \tag{5b}
\]
The scalar \(\ell\) may be negative. The improvement uses the fixed
identity \(\|\operatorname{vec}B\|_2^2=n^2\), not any extra property
of the rounding distribution. Section 6 optimizes the scalar-shift
subclass in its conditional conference example; it does not optimize
all diagonal shifts in (5a).

## 2. Ridge and spectral-tail bounds

For \(\lambda>0\), take
\[
G=\frac{\mathcal K+\lambda I}{1+\lambda}.
\]
It is positive definite with diagonal one. Equation (5) gives
\[
\min_B\overline q(B)
 \le\inf_{\lambda>0}(1+\lambda)
       \operatorname{tr}\left[\mathcal K(\mathcal K+\lambda I)^{-1}\right].
                                                                    \tag{6}
\]

Suppose \(\mathcal K=M+E_0\), where \(M,E_0\succeq0\),
\(\operatorname{rank}M=d\), and \(\operatorname{tr}E_0=\eta\).
The sum of the eigenvalues of \(\mathcal K\) after its largest \(d\)
is at most \(\eta\): projection onto \(\operatorname{range}M\)
already captures at least \(\operatorname{tr}M\).
Thus
\[
\operatorname{tr}\left[\mathcal K(\mathcal K+\lambda I)^{-1}\right]
 \le d+\eta/\lambda.
\]
Optimizing \((1+\lambda)(d+\eta/\lambda)\), including limiting cases
when \(d\) or \(\eta\) is zero, yields
\[
\boxed{\quad
\min_B\overline q(B)\le
 \min\{n^2,(\sqrt d+\sqrt\eta)^2\}.
\quad}                                                        \tag{7}
\]

In particular, let \(U_r,V_s\) be spectral truncations of ranks \(r,s\),
and set
\[
u=\operatorname{tr}(U-U_r),\quad
v=\operatorname{tr}(V-V_s),\quad
h=\dim(\operatorname{range}U_r\cap\operatorname{range}V_s).
\]
Take
\[
M=\tfrac12(V_s\otimes U_r+U_r\otimes V_s).
\]
The remainder \(\mathcal K-M\) is positive semidefinite by expanding
each tensor product using its positive semidefinite spectral tail.
Moreover,
\[
d=2rs-h^2,\qquad \eta=n(u+v)-uv.                               \tag{8}
\]
For the rank formula, the two tensor-product ranges intersect in
\((\operatorname{range}U_r\cap\operatorname{range}V_s)^{\otimes2}\);
the range of a sum of positive semidefinite matrices is the sum of
their ranges. The trace identity follows from
\(\operatorname{tr}M=(n-u)(n-v)\).

Finite-temperature full-support Ising covariances have full rank.
A bare rank bound is therefore not automatically informative; the
quantitative tail traces in (7)--(8) are essential.

## 3. A centered bound using the actual Gibbs energy

Now let \(A\) be a nonzero symmetric zero-diagonal matrix, and let
\(U,V\) be its actual zero-field Gibbs covariances at \(+\beta,-\beta\),
with \(\beta>0\). Define
\[
Q_A(x)=\tfrac12x^TAx,\quad Z_A(t)=\mathbb E_xe^{tQ_A(x)},\quad
a_A(\beta)=\tfrac12(\log Z_A(\beta)+\log Z_A(-\beta)).
\]
The uniform mean of \(Q_A\) is zero. Convexity of \(\log Z_A\) gives
\[
p:=\operatorname{tr}(AU)\ge0,\qquad
q:=\operatorname{tr}(AV)\le0,\qquad
a_A'(\beta)=\frac{p-q}{4}>0.                                   \tag{9}
\]
Put
\[
L=\|A\|_{\rm op},\quad
\alpha=\frac{p+q}{2n},\quad
H=A\otimes A-\alpha(A\otimes I+I\otimes A),\quad
E=\frac{p^2+q^2}{2}.                                          \tag{10}
\]
The zero diagonal of \(A\) implies \(\operatorname{diag}H=0\).
Direct tensor trace identities give
\[
\operatorname{tr}(\mathcal K H)
 =pq-\alpha n(p+q)=-E,\qquad E\ge4(a_A')^2.                    \tag{11}
\]

Write \(r_2=\operatorname{tr}(A^2U)\) and
\(s_2=\operatorname{tr}(A^2V)\). Expanding \(H^2\) gives
\[
\begin{aligned}
F:=\operatorname{tr}(\mathcal K H^2)
 &=r_2s_2-2\alpha(qr_2+ps_2)
   +\alpha^2[n(r_2+s_2)+2pq].                                \tag{12}
\end{aligned}
\]
We claim
\[
0\le F\le n^2L^4,\qquad \|H\|_{\rm op}\le2L^2,\qquad
E\le n^2L^2.                                                   \tag{13}
\]
For the first upper bound, the derivative of the expression in (12)
with respect to \(r_2\), holding \(p,q,\alpha,s_2\) fixed, is
\[
s_2-2\alpha q+\alpha^2n
 =\operatorname{tr}((A-\alpha I)^2V)\ge0.
\]
The analogous derivative with respect to \(s_2\) is nonnegative.
First increase \(r_2\) to \(nL^2\); then increase \(s_2\) to \(nL^2\).
The latter derivative can only increase during the first step.
The result is
\[
F\le n^2L^4-2\alpha^2(n^2L^2-pq)\le n^2L^4.
\]
Nonnegativity follows from \(\mathcal K\succeq0\) and \(H^2\succeq0\).
Since \(0\le p\le nL\) and \(-nL\le q\le0\),
\(|\alpha|\le L/2\), which proves the operator bound in (13);
the last assertion is immediate.

Let \(h_0=2L^2\). For \(0<t<h_0^{-1}\), the matrix
\(G=I-tH\) is positive definite with diagonal one.
The exact resolvent identity and functional calculus give
\[
G^{-1}=I+tH+t^2H^2G^{-1},\qquad
H^2G^{-1}\preceq\frac{H^2}{1-th_0}.
\]
Equations (5) and (11) therefore imply
\[
\min_B\overline q(B)
 \le n^2-tE+\frac{t^2F}{1-th_0}.                               \tag{14}
\]
Here \(E>0\), and Cauchy--Schwarz gives \(E^2\le n^2F\), so \(F>0\).
Choose \(t=E/(2F+h_0E)\). The remainder in (14) is half of the
linear saving, and (13) now proves
\[
\boxed{\quad
\min_B\overline q(B)
 \le n^2-\frac{E^2}{4F+2h_0E}
 \le n^2-\frac{E^2}{8n^2L^4}
 \le n^2-\frac{2(a_A'(\beta))^4}{n^2L^4}.
\quad}                                                        \tag{15}
\]
For the degenerate \(E=0\) case in a more general covariance setting,
one simply uses independent signs.

This is a strict finite-temperature improvement over random cross-block
signs, tied to the actual Gibbs energy. It can be a leading \(n^2\)
saving if \(L=O(\sqrt n)\) and \(a_A'=\Omega(n^{3/2})\).
Neither of those two hypotheses is silently supplied for the desired
global-minimizer application, and (15) is not the comparison (2).

### 3.1. A stronger elementary Gaussian-sign bound for complete signings

Let \(\mathcal S_n\) be the set of symmetric zero-diagonal matrices
with every off-diagonal entry in \(\{-1,1\}\). Now specialize to
\(A\in\mathcal S_n\), \(n\ge2\), retaining the actual Gibbs
covariances and notation of Section 3. There is a stronger integral
bound in this case that does not require the imported theorem.

Only the negative spectrum of \(H\) matters for the construction.
Write the extremal eigenvalues of \(A\) as \(a>0\) and \(-b<0\),
so \(L=\max(a,b)\). Actual phase energies give
\(0\le p\le na\) and \(-nb\le q\le0\), hence
\(-b/2\le\alpha\le a/2\). The eigenvalues of \(H\) are
\(h(x,y)=xy-\alpha(x+y)\), for eigenvalues \(x,y\) of \(A\).
The bilinear function \(h\) is minimized on \([-b,a]^2\) at a
corner. Its two same-sign corner values satisfy
\[
h(a,a)=a^2-2\alpha a\ge0,\qquad
h(-b,-b)=b^2+2\alpha b\ge0,
\]
while its two mixed corner values equal \(-\mu\), where
\[
\mu=ab+\alpha(a-b)>0,\qquad
\mu\le\frac{\max(a,b)(a+b)}2\le L^2.
\]
Both extremal eigenvalues occur, so \(\lambda_{\min}(H)=-\mu\).
In particular, \(H\succeq-L^2I\). This is a bound on its negative
spectrum, not the assertion \(\|H\|_{\rm op}\le L^2\).

Set \(\tau=1/\mu\). The matrix
\(\Sigma=I+\tau H\) is positive semidefinite with diagonal one.
Let \(z\) be a centered Gaussian vector with covariance \(\Sigma\)
and define the integral cross block by
\(\operatorname{vec}B=\operatorname{sign}z\). Singular covariance
is allowed; every coordinate has variance one and is nonzero almost
surely. For standard Gaussian coordinates of correlation \(\rho\),
\[
\mathbb E[\operatorname{sign}z_i\operatorname{sign}z_j]
 =\frac2\pi\arcsin\rho.
\]
This identity follows by rotational invariance: the probability that
the two Gaussian linear forms have different signs is
\(\arccos(\rho)/\pi\), with endpoints obtained by continuity.

Under column vectorization, distinct tensor-coordinate pairs have
three possibilities: if both component indices differ, the entry of
\(H\) is the product of two off-diagonal signs; if exactly one index
is equal, the entry is \(-\alpha\) times one off-diagonal sign;
and if both are equal, the diagonal entry is zero. The supports of
these terms are disjoint. The oddness of scalar arcsine therefore
gives the exact sign-covariance identity
\[
\mathbb E[\operatorname{vec}B\operatorname{vec}B^T]
=I+\frac2\pi\left[
 \arcsin(\tau)(A\otimes A)
 -\arcsin(\tau\alpha)(A\otimes I+I\otimes A)\right].
\]
This is entrywise scalar arcsine, not matrix functional calculus.
Taking the trace against \(\mathcal K\) gives
\[
\mathbb E\overline q(B)
=n^2+\frac2\pi\left[
 \arcsin(\tau)pq
 -\arcsin(\tau\alpha)n(p+q)\right]
\le n^2-\frac{2\tau E}{\pi}.
\]
For the last inequality, \(pq\le0\),
\(\arcsin(\tau)\ge\tau\), and \(\alpha\) has the same sign as
\(p+q\). Thus replacing each arcsine by its argument has the stated
upper direction, using \(|\arcsin u|\ge|u|\).
All arcsine arguments belong to \([-1,1]\) because they are entries
of the correlation matrix \(\Sigma\).
Selecting an outcome proves
\[
\boxed{\quad
\min_B\overline q(B)
 \le n^2-\frac{2E}{\pi\mu}
 \le n^2-\frac{2E}{\pi L^2}
 \le n^2-\frac{8(a_A'(\beta))^2}{\pi L^2}.
\quad}                                                        \tag{15a}
\]
One may alternatively take \(\tau=1/L^2\), since \(\mu\le L^2\).
The displayed estimate is stronger than the centered bound (15) for
complete signings, but is not asserted to improve
every possible covariance-rounding certificate. The signing is fixed
throughout the Gaussian construction; no independence assumption for
a noise-dependent optimizer is used. Neither (15) nor (15a) alone
establishes the desired radial comparison (2).

## 4. Elementary skew rounding and the endpoint normalization

A second direct bound requires no imported rounding theorem. Choose
independent fair signs \(B_{ij}\) for \(i<j\), put \(B_{ji}=-B_{ij}\),
and choose the diagonal signs independently. The only nontrivial
off-coordinate covariances are
\(\mathbb E B_{ij}B_{ji}=-1\). Expanding (1) gives
\[
\mathbb E\overline q(B)=n^2+n-\operatorname{tr}(UV),
\quad
\min_B\overline q(B)\le
 \min\{n^2,n^2+n-\operatorname{tr}(UV)\}.                        \tag{16}
\]

For the remaining statements involving signing optimality, retain
\(A\in\mathcal S_n\) as defined in Section 3.1.
At an actual edge-local half-product minimum in this set, the exact
single-edge flip inequality is
\[
A_e(U_e-V_e)\le\tanh(2\beta)(1-U_eV_e).
\]
Summing over unordered edges gives
\[
2a_A'(\beta)\le
 \frac{\tanh(2\beta)}2\bigl(n^2-\operatorname{tr}(UV)\bigr).
                                                                    \tag{17}
\]
The direction is important: (17) bounds the desired radial quantity
from above by approximately the overlap expression in (16). It does
not turn the skew bound into (2).

For clarity, take \(\beta=c/\sqrt n\), write
\(R_n(\beta)=\min_{A\in\mathcal S_n}a_A(\beta)\), and define the
optimized balanced cosh path exactly by
\[
\begin{split}
f_{2n}(t)=\min_{\substack{A_1,A_2\in\mathcal S_n\\
                         B\in\{-1,1\}^{n\times n}}}
\log\mathbb E_{x,y}\cosh\Bigl[
 &\beta\sqrt{1-t/2}\bigl(Q_{A_1}(x)+Q_{A_2}(y)\bigr)\\
 &+\beta\sqrt{t/2}\,x^TBy\Bigr],\qquad 0\le t\le1.
\end{split}                                                    \tag{18a}
\]
Here \(x,y\) are independent uniform signs. This is a minimum of a
log-cosh expectation, not a half-product pressure at order \(2n\).
At \(t=0\), the expression inside the logarithm is
\[
\tfrac12\bigl[Z_{A_1}(\beta)Z_{A_2}(\beta)
             +Z_{A_1}(-\beta)Z_{A_2}(-\beta)\bigr]
\ge \exp\bigl(a_{A_1}(\beta)+a_{A_2}(\beta)\bigr)
\ge e^{2R_n(\beta)}.
\]
For any minimizing \(A\), the choice \((A_1,A_2)=(A,-A)\) attains
equality. Thus \(f_{2n}(0)=2R_n(\beta)\). In this active paired
branch with cross block \(B\), the first right derivative is
\[
-\frac{\beta}{2}a_A'(\beta)+\frac{\beta^2}{4}\overline q(B).
                                                                    \tag{18}
\]
Indeed, the internal coefficient derivative is \(-\beta/4\), and
the signed internal energy at the paired endpoint is \(2a_A'\).
The cross mean vanishes; its second-moment coefficient is
\((\beta^2/2)/2=\beta^2/4\), with the two opposite-phase orderings
averaged exactly as in (1).

The fully optimized order-\(2n\) path is no larger than this branch and
agrees with it at the active endpoint. Consequently,
\[
f_{2n}'(0+)\le
-\frac{\beta}{2}a_A'(\beta)
+\frac{\beta^2}{4}\min_B\overline q(B).                         \tag{19}
\]
This is only an endpoint upper bound; the exact derivative minimizes
over all active block pairs as well. No assertion about the derivative
away from zero, an integrated path error, or convergence follows from
(19) alone. In particular, even proving (2) at every order would leave
the integrated comparison \(f_{2n}(1)\le f_{2n}(0)+o(n)\)
unjustified by this endpoint argument.

## 5. A sufficient low-effective-rank case

Let \(S=(U+V)/2\). Suppose a rank-\(r\) spectral truncation of \(S\)
has tail trace \(\tau\). Since \(U,V\preceq2S\), each of \(U,V\)
has a rank-\(r\) spectral truncation with tail trace at most \(2\tau\).
To verify this, let \(P\) project onto the chosen top subspace of \(S\);
then \(\operatorname{tr}((I-P)U)\le2\tau\), and the best rank-\(r\)
spectral subspace of \(U\) captures at least as much trace as \(P\).
The same applies to \(V\).

Using (7)--(8), we obtain the explicit sufficient estimate
\[
\min_B\overline q(B)
 \le\left(\sqrt2\,r+2\sqrt{n\tau}\right)^2.                     \tag{20}
\]
Thus \(r=o(n)\) and \(\tau=o(n)\) imply
\(\min_B\overline q(B)=o(n^2)\), which suffices for the endpoint
inequality (2), since \(a_A'\ge0\). No entropy-to-low-rank statement
or implication from signed correlation diffuseness is being assumed.

## 6. What remains missing in the high-rank regime

The signed edge-local estimate
\(\sum_{i<j}|U_{ij}-V_{ij}|=O_c(n^{3/2})\) does not itself establish
the low-effective-rank hypothesis of Section 5. A useful algebraic
check is the conditional conference-form model with
\(A\in\mathcal S_n\) a symmetric conference sign matrix:
\[
A^2=(n-1)I,\qquad
U=I+tA,\quad V=I-tA,\quad s=t\sqrt{n-1},\quad 0\le s\le1.
                                                                    \tag{21}
\]
This paragraph makes no claim that arbitrary chosen \(t\) gives actual
Gibbs covariances, or that \(A\) is globally minimizing.

In (21), \(S=I\), signed diffuseness holds for bounded \(s\), and
\[
\overline q(B)=n^2-t^2\operatorname{tr}(B^TAB A).
\]
As \(A^2=(n-1)I\),
\(\operatorname{tr}(B^TAB A)\le(n-1)n^2\). Equality is attained by
the integral matrix \(B=A+I\). Hence
\[
\min_B\overline q(B)=n^2(1-s^2).                               \tag{22}
\]
If this covariance model is actual Gibbs data, its radial derivative is
\(a_A'=tn(n-1)/2\). At \(\beta=c/\sqrt n\), comparison (2) would require
asymptotically
\[
1-s^2\le s/c.
\]
The scalar edge-local inequality reduces only to \(t\le\tanh\beta\),
or \(s\le c+o(1)\), which is an upper response bound, not the needed
lower response estimate.

One can compare the unshifted bound (5) and its scalar-shift refinement
(5b) exactly in this model. These calculations are not a limitation
theorem for every use of Gram--Schmidt rounding. Put
\[
J_0=(A\otimes A)/(n-1),\qquad
\mathcal K=I-s^2J_0.
\]
Then \(J_0^2=I\), its diagonal and trace are zero, and its two
eigenspaces have equal dimension. For every admissible \(G\),
Cauchy--Schwarz gives
\[
\operatorname{tr}(\mathcal K G^{-1})
 \ge\frac{(\operatorname{tr}\sqrt{\mathcal K})^2}{n^2}.
\]
For \(s<1\), equality is attained by
\(G=\sqrt{\mathcal K}/\gamma\), where
\(\gamma=\operatorname{tr}\sqrt{\mathcal K}/n^2\):
the diagonal of \(\sqrt{\mathcal K}\) is the constant \(\gamma\).
Taking a limit covers \(s=1\). Therefore
\[
\inf_{\substack{G\succ0\\\operatorname{diag}G\le1}}
 \operatorname{tr}(\mathcal K G^{-1})
 =\frac{n^2}{2}\left(1+\sqrt{1-s^4}\right).                     \tag{23}
\]
Equation (23) describes only the unshifted admissible family (5).
Using the fixed-square refinement makes its value strictly smaller
when \(0<s<1\). To see the full scalar-shift optimum, let
\(a=1-\ell\ge s^2\). The residual matrix
\(\mathcal K-\ell I=aI-s^2J_0\) has eigenvalues \(a-s^2,a+s^2\),
each with multiplicity \(n^2/2\). The same trace inequality and
normalized-square-root construction, with limits at zero eigenvalues,
give
\[
\begin{split}
\inf_{\substack{G\succ0\\\operatorname{diag}G\le1}}
\left[n^2\ell+
 \operatorname{tr}\bigl((\mathcal K-\ell I)G^{-1}\bigr)\right]
 &=n^2\left[1-\frac a2+\frac12\sqrt{a^2-s^4}\right].
\end{split}
\]
For \(s>0\), this expression is increasing in \(a>s^2\), so its
minimum occurs at \(a=s^2\), or \(\ell=1-s^2\). For \(s=0\),
it is identically \(n^2\). Consequently
\[
\boxed{\quad
\inf_{\substack{\ell\le1-s^2,\ G\succ0\\\operatorname{diag}G\le1}}
\left[n^2\ell+
 \operatorname{tr}\bigl((\mathcal K-\ell I)G^{-1}\bigr)\right]
 =n^2\left(1-\frac{s^2}{2}\right).
\quad}                                                        \tag{24}
\]
For \(s>0\) the boundary infimum over positive-definite \(G\) is
not attained: at \(\ell=1-s^2\), one may take
\(G_\varepsilon=I-(1-\varepsilon)J_0\), \(0<\varepsilon<1\).
Its certificate value is
\(n^2(1-s^2)+n^2s^2/(2-\varepsilon)\), tending to (24).
Passing to an infimum is valid for the finite integral minimum.
For \(s=0\), \(G=I\) attains the value \(n^2\).
No optimization over nonconstant diagonal shifts \(\Lambda\) is
asserted here.

For \(0<s\le1\), (24) is still strictly above the exact integral
minimum \(n^2(1-s^2)\) from (22); both (23) and (24) have lower
endpoint value \(n^2/2\) at \(s=1\). This is a comparison of these
two specific bound families with an explicit integral construction,
not a no-go theorem for the full diagonal-shift class (5a), other
rounding arguments, or actual Gibbs covariances.

Indeed, the Gaussian construction of Section 3.1 also applies to
these conditional covariance matrices: here \(\alpha=0\) and
\(\mu=L^2=n-1\). It gives
\(\min_B\overline q(B)\le n^2(1-2s^2/\pi)\), already strictly
below (24) for \(s>0\), though still above the exact value (22).
This further distinguishes the trace certificates (23)--(24) from
all integral rounding constructions.

For the actual global-minimizer application, the unresolved endpoint
step is still (2): a useful integral cross block must be tied to the
actual radial energy. The conditional model shows why a lower response
estimate would matter even when the cross-block minimum is known
exactly; it does not supply an actual-Gibbs counterexample. Beyond
that endpoint step, the integrated order comparison after (19)
also remains open for this route.
