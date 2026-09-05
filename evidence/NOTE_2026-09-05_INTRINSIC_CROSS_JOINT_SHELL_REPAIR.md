# Intrinsic cross normalization and a rank-four joint-shell repair

2026-09-05. This extends the upper in
original_mo_conditional_cross_joint_shell_upper.md to the intrinsic
cross covariance. The independent Gaussian cushion is retained. The
upper is not evaluated sharply enough to prove original convergence.

## 1. Intrinsic covariance

Let n>=2, let A be a complete symmetric zero-diagonal order-n signing,
and let B be an n by n sign matrix. Put

\[
 K=\begin{pmatrix}A&B\\B^T&-A\end{pmatrix},\quad
 H=A\otimes A-\mathcal S_B+I,\quad
 (\mathcal S_B)_{ij,kl}=B_{il}B_{kj},
\]
\[
 \mu=\max(2,\|H\|_{\rm op}),\quad \ell=\mu+1,\quad
 R=I+H/\mu={\ell I+A\otimes A-\mathcal S_B\over\mu}.       \tag{1}
\]

Then \(0\le R\le2I\), diag R=1, and \(\mu\ge n-1\).
The principal H submatrix on a fixed row is \(I-b_i b_i^T\), proving
the last assertion. We do not assume \(\ell\ge\|K\|_{\rm op}^2\).

Let Z be a centered Gaussian cross matrix with covariance \(C=kR+vI\),
where k,v>=0. For the actual shifted-sign proposal,
\(s=2\Phi_G(h)-1,\ k=4\phi(h)^2,\ v=1-s^2-k\), so \(k+v\le1\).

Fix a nonempty shell
\[
 \mathcal T_{p,q,c}=\{(x,y)\in\{-1,1\}^{2n}:
           x^TAx=p,\ y^TAy=q,\ x^TBy=c\},                \tag{2}
\]
and choose a deterministic representative \((\bar x,\bar y)\).
Write \(L_A=\|A\|_{\rm op},L_B=\|B\|_{\rm op}\), and
\[
 a_0=n(k\ell/\mu+v),\quad t=k/\mu,\quad
 M=\begin{pmatrix}a_0I+tqA&-tcB\\-tcB^T&a_0I+tpA\end{pmatrix}.
                                                               \tag{3}
\]
M is not assumed positive semidefinite.

## 2. A direct positive rank-four repair

Define
\[
 P={k\over\mu}\begin{pmatrix}
 \ell\bar x\bar x^T+(A\bar x)(A\bar x)^T&0\\
 0&\ell\bar y\bar y^T+(A\bar y)(A\bar y)^T
 \end{pmatrix},\qquad \widetilde M=M+P.                  \tag{4}
\]
Then P>=0, rank P<=4, and \(\widetilde M\ge0\).
For arbitrary real u,v set \(U=u\bar y^T+\bar x v^T\).
Expanding the quadratic form of R gives
\[
\begin{aligned}
 \mu\langle U,RU\rangle
 ={}&\ell n(\|u\|^2+\|v\|^2)
       +2\ell(u\cdot\bar x)(v\cdot\bar y)\\
 &+q\,u^TAu+p\,v^TAv+2(u^TA\bar x)(v^TA\bar y)\\
 &-(u^TB\bar y)^2-(\bar x^TBv)^2-2c\,u^TBv.
\end{aligned}
\]
Consequently
\[
\begin{aligned}
 \mu (u,v)^T\widetilde M(u,v)
 ={}&k\mu\langle U,RU\rangle
       +k(u^TB\bar y)^2+k(\bar x^TBv)^2\\
 &+k\ell(u\cdot\bar x-v\cdot\bar y)^2
       +k(u^TA\bar x-v^TA\bar y)^2\\
 &+\mu vn(\|u\|^2+\|v\|^2)\ \ge0.                       \tag{5}
\end{aligned}
\]
This proof uses R>=0 directly. It does not use a whole-source
operator-norm denominator.

Let \(g=(\xi,\eta)\) be Gaussian with covariance \(\widetilde M\).
Because M may be indefinite, g must NOT be described as a Gaussian
with covariance M plus an independent Gaussian with covariance P.

## 3. Actual Gaussian increment comparison

For two states in one shell write
\[
 r_x=x^Tx',\quad r_y=y^Ty',\quad
 d=x^TBy',\quad e=x'^TBy.
\]
The algebraic increment excess of M over the covariance of
\(X_{x,y}=x^TZy\), divided by two, is exactly
\[
 {k\over4\mu}\left[
 \mu\langle(x-x')\otimes(y-y'),
 R((x-x')\otimes(y-y'))\rangle+(d-e)^2\right]
       +v(n-r_x)(n-r_y)\ \ge0.                          \tag{6}
\]
This is the previously checked full-exchange identity with
\(L^2=\ell,D=\mu\); alternatively direct expansion proves it.
Equation (6) is algebraic and makes no covariance claim about M.
Adding P adds a nonnegative increment. Since \(\widetilde M\) is PSD,
Gaussian increment comparison now proves
\[
 \boxed{\mathbb E\sup_{\mathcal T_{p,q,c}}x^TZy
 \le\widetilde w_{p,q,c}:=
       \mathbb E\sup_{\mathcal T_{p,q,c}}(\xi^Tx+\eta^Ty).} \tag{7}
\]

## 4. Exact nonuniform-marginal mismatch bound

If a_0=0 all noise and widths are zero. Otherwise let
\[
 a_i=a_0+P_{ii},\quad a'_j=a_0+P_{n+j,n+j},\quad
 E_x={\operatorname{tr}P_x\over a_0},\quad
 E_y={\operatorname{tr}P_y\over a_0},\quad
 E_*=1+{L_A^2\over\ell}.
\]
The ratio \(k/(\mu a_0)\le1/(n\ell)\) proves
\[
                 0\le E_x,E_y\le E_*.                  \tag{8}
\]

Let \(\widetilde m_p,\widetilde m_q,\widetilde m_c\) be the exact
expected three energies of \(z^0=(\operatorname{sign}\xi,
\operatorname{sign}\eta)\), computed by the Gaussian arcsine identity.
Set
\[
 \widetilde b=((p-\widetilde m_p)/2,
              (q-\widetilde m_q)/2,c-\widetilde m_c),
\]
\[
 H_\lambda=\begin{pmatrix}\lambda_1A&\lambda_3B\\
                   \lambda_3B^T&\lambda_2A\end{pmatrix},\quad
 R_\lambda={\lambda_{\max}(H_\lambda)-\lambda_{\min}(H_\lambda)\over2},
 \quad
 \widetilde r=\sup_{\lambda\ne0}
                  {|\lambda\cdot\widetilde b|\over2nR_\lambda}.
                                                               \tag{9}
\]
For n>=2, \(R_\lambda>0\) when lambda is nonzero. Define
\[
 \kappa=2/\pi,\quad \delta(r)={1-\sqrt{1-r^2}\over2},\quad
 z_d=\Phi_G^{-1}((1+d)/2),\quad
 \mathcal L(d)=1-e^{-z_d^2/2}.
\]
Then \(0\le\widetilde r\le1\), and the exact finite-n upper is
\[
 \boxed{\widetilde w_{p,q,c}\le
 2n\sqrt{\kappa a_0}[1-2\mathcal L(\delta(\widetilde r))]
             +{\sqrt{\kappa a_0}\over2}(E_x+E_y).}       \tag{10}
\]

For completeness, let z* maximize the field on the finite shell and
let \(\bar d=\mathbb E\operatorname{Ham}(z^*,z^0)/(2n)\).
Subtracting the spectral midpoint of H_lambda and polarizing gives
\[
 |\lambda\cdot\widetilde b|
       \le4nR_\lambda\sqrt{\bar d(1-\bar d)}.
\]
Thus \(\bar d\ge\delta(\widetilde r)\).
Every field coordinate has variance at least a_0, so
\[
 \mathbb E(u-|g_i|)_+
       \le\mathbb E(u-\sqrt{a_0}|G|)_+,\quad u\ge0.
\]
The threshold argument for arbitrary jointly selected flip indicators,
optimized at \(u=\sqrt{a_0}z_{\bar d}\), therefore gives
\[
 \sum_i\mathbb E[|g_i|1_{z_i^*\ne z_i^0}]
             \ge2n\sqrt{\kappa a_0}\mathcal L(\bar d).
\]
Meanwhile, concavity of square root gives
\[
 \sum_i\mathbb E|g_i|
 \le2n\sqrt{\kappa a_0}
                  +{\sqrt{\kappa a_0}\over2}(E_x+E_y).
\]
Subtracting twice the flip cost proves (10). No coordinate or
selection independence has been used.

## 5. Reference-energy stability, uniform even as the noise vanishes

Assume, for fixed constants \(C_A,C_K\),
\[
 \Phi(A)\le C_A n^{3/2},\qquad \Phi(K)\le C_K n^{3/2}.     \tag{11}
\]
Every shell then has \(|p|,|q|\le2C_A n^{3/2}\) and
\(|c|\le C_K n^{3/2}\). We use
\[
                         L_A^2\le8\Phi(A).             \tag{12}
\]
Indeed real cube polarization gives
\(\beta_{\mathbb R}(A)\le4\Phi(A)\). The complex infinity-to-one norm
is at most \(2\beta_{\mathbb R}(A)\): rotate a bilinear value to the
real axis and split real and imaginary parts. The complex one-to-
infinity norm is one, and finite-dimensional Riesz--Thorin interpolation
gives (12).

Put
\[
 \rho_p=tq/a_0,\quad \rho_q=tp/a_0,\quad \rho_c=tc/a_0,
\]
\[
 m_p^0=\kappa n(n-1)\arcsin\rho_p,\quad
 m_q^0=\kappa n(n-1)\arcsin\rho_q,\quad
 m_c^0=-\kappa n^2\arcsin\rho_c.                         \tag{13}
\]
For all sufficiently large n these correlations have absolute value
at most one half, because each is bounded by its corresponding
energy divided by \(n\ell\). Uniformly over shells, representatives,
and all k,v>=0 with a_0>0,
\[
 \boxed{
 |\widetilde m_p-m_p^0|+|\widetilde m_q-m_q^0|
                \le C L_A\sqrt n=O(n^{5/4}),\qquad
 |\widetilde m_c-m_c^0|\le Cn.}                          \tag{14}
\]

Here is a quantitative proof that retains variance normalization.
In the first block set
\[
 \alpha_i=P_{ii}/a_0,\quad w_i=(1+\alpha_i)^{-1/2},\quad
 W=\operatorname{diag}(w_i),\quad N=W(P_x/a_0)W.
\]
Then N>=0, \(\operatorname{tr}N\le E_x\), and
\(\sum_i(1-w_i)\le E_x/2\). Its off-diagonal correlation matrix is
\((\rho_p WAW+N)_{\rm off}\).
The linear arcsine term differs from \(\rho_p n(n-1)\) by at most
\[
                         (|\rho_p|n+L_A)E_x.            \tag{15}
\]
For the first summand use
\(\sum_{i\ne j}(1-w_iw_j)\le(n-1)E_x\); for the other use
\(|\operatorname{tr}(AN)|\le L_A\operatorname{tr}N\).
The universal inequality \(|\arcsin u-u|\le C u^2\) on [-1,1] and
\[
 \|(\rho_p WAW+N)_{\rm off}\|_F^2
                   \le2\rho_p^2n(n-1)+2E_x^2
\]
show that the full error against the FIRST formula in (13) is at most
\[
 C[(|\rho_p|n+L_A)E_x+\rho_p^2n^2+E_x^2].                \tag{16}
\]
The baseline arcsine-minus-linear error is included in this bound.
The second block has the same estimate. By (11)--(12),
\(E_x,E_y=O(\sqrt n)\), \(|\rho_p|,|\rho_q|=O(n^{-1/2})\), and
\(L_A\ge\sqrt{n-1}\), so (16) is at most \(C L_A\sqrt n\).

The cross block has no repair entry. If w'_j are the second-block
weights, its exact mean is
\[
       \widetilde m_c=-\kappa\sum_{i,j}
                            \arcsin(\rho_c w_iw'_j).
\]
The derivative of arcsin on [-1/2,1/2] is at most two, and hence
\[
 |\widetilde m_c-m_c^0|
                \le\kappa|\rho_c|n(E_x+E_y)=O(n).        \tag{17}
\]
This proves (14). All constants depend only on \(C_A,C_K\), not on
k,v or h. In particular no lower bound on a_0 or k has been used:
the ratios are bounded uniformly as either tail \(h\to\pm\infty\)
makes the actual noise vanish.

## 6. The same leading mismatch formula at intrinsic mu

Let
\[
 b^0=((p-m_p^0)/2,(q-m_q^0)/2,c-m_c^0),\quad
 r_0=\sup_{\lambda\ne0}{|\lambda\cdot b^0|\over2nR_\lambda},
 \qquad \bar r_0=\min(1,r_0).
\]
Uniformly under (11),
\[
                    |r_0-\widetilde r|\le C n^{-1/2}.  \tag{18}
\]
Indeed H_lambda has trace zero, so
\(\|H_\lambda\|_{\rm op}\le2R_\lambda\). Its three blocks give
\[
 |\lambda_1|L_A,\ |\lambda_2|L_A,\ |\lambda_3|L_B
                                      \le2R_\lambda.
\]
Use (14), \(L_B\ge\sqrt n\), and divide the numerator difference
by \(2nR_\lambda\). Clipping r_0 to [0,1] cannot increase its
distance from \(\widetilde r\).

The function \(\mathcal L(\delta(r))\) is uniformly Holder continuous
of exponent one half on [0,1]. This follows from
\[
 |\delta(r)-\delta(r')|\le\sqrt{|r-r'|/2}
\]
and the bounded derivative of \(\mathcal L\) on [0,1/2].
Since \(a_0\le2n\) when \(k+v\le1\), (7), (10), and (18) prove
\[
 \boxed{\mathbb E\sup_{\mathcal T_{p,q,c}}x^TZy
 \le2n\sqrt{\kappa a_0}
       [1-2\mathcal L(\delta(\bar r_0))]+C n^{5/4}.}      \tag{19}
\]
The direct marginal inflation in (10) is only O(n). The larger
O(n^(5/4)) error allows the square-root endpoint loss in converting
the exact repaired mismatch to its baseline value.

Thus the opposite cross reference mean is retained, now with
\[
             \rho_c={kc\over n(k(\mu+1)+v\mu)},           \tag{20}
\]
rather than a whole-source spectral-square denominator.

## 7. The actual Gaussian maximum and the unevaluated step

There are at most \(J=(2n^2+1)^3\) nonempty shells. Gaussian
concentration, \(\|C\|_{\rm op}\le2k+v\), and (7) give for real s
\[
 \mathbb E\Phi\begin{pmatrix}A&sB+Z\\sB^T+Z^T&-A\end{pmatrix}
 \le\max_{\mathcal T_{p,q,c}\ne\varnothing}
       \left[\left|{p-q\over2}+sc\right|
                          +\widetilde w_{p,q,c}\right]
                +n\sqrt{2(2k+v)\log(2J)}.               \tag{21}
\]
The right side of (10) or (19) may replace the shell width.
A shell maximum in either absolute-value phase is n-Lipschitz in
the cross Gaussian vector, giving the displayed concentration
remainder. Representatives may differ between shells.

If B is an actual conditional cross optimizer, every shell also
satisfies \(|(p-q)/2|+|c|\le F_A^*\). The separate intrinsic-
covariance sign-to-Gaussian theorem supplies the conditional floor.
This note does not substitute the whole-source covariance for that
distinct theorem.

The verified new implication is that intrinsic cross positivity
suffices for a repaired, leading-order-equivalent joint-shell upper.
What remains is a sharp evaluation of the upper using attainable
shells and conditional optimality. Neither this repair nor the
opposite mismatch alone proves
\(F_A^*\le2\sqrt2\,\Phi(A)+o(n^{3/2})\).

