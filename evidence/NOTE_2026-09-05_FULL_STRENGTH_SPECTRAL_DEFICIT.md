# A spectral-deficit criterion for the full-strength Gaussian-sign law

2026-09-05. **All-orders conditional obstruction for the singular
canonical law; exact-minimizer spectral hypothesis OPEN.**

Fix (c>0), (0<t\le1), and a small (0\le r<1\), independently
of (n). Put (d=n^2), \(\gamma=c\sqrt{t/(2n)}\), and

\[
 K_0=\frac4{3\sqrt\pi},\qquad
 g_r(c,t)=c(\sqrt{2t}K_0-1)-2\log2
                  -\frac{c^2t}{2\pi}\arcsin r.
\]

There are fixed positive (r,c,t) with (g_r(c,t)>0), by the
reviewed conditional-noise theorem. Throughout, the covariance-generating
host is fixed before the Gaussian and its signs are drawn.

Let \(\Sigma\succeq0\) be any (d\times d) correlation matrix,
including a singular one. Define its deficit below the nearly independent
Gaussian variance level (1-r\) by

\[
 D_r=[(1-r)I_d-\Sigma]_+,\qquad V_r=\operatorname{tr}D_r.
                                                               \tag{1}
\]

Here the positive part is spectral functional calculus. Generate
\(G\sim N(0,\Sigma)\) and (B=\operatorname{sign}G\), reshaped as
an (n\times n\) cross block. For the actual paired pressures and
endpoint notation of the quenched note,

\[
\boxed{\begin{aligned}
 \mathbb E\min_A F_{A,B}(t)
 \ge{}&[c\sqrt{2t}K_0-2\log2
                   -c^2t\arcsin r/(2\pi)]n\\
 &-\frac{c\sqrt{2t}}\pi\sqrt{nV_r}
  -\frac{c^2t}{2\pi n}V_r-o_{c,t}(n).
\end{aligned}}                                               \tag{2}
\]

The error is uniform over \(\Sigma\) and fixed (r\). In particular,
if (g_r(c,t)>0\) and (V_r=o(n)\), then this law has a positive
linear mean excess over \(2R_n(c/\sqrt n)\). For the actual canonical
construction \(\Sigma=I+H/\mu\), including its singular endpoint,
the sufficient spectral hypothesis is exactly

\[
 \operatorname{tr}[-rI-H/\mu]_+=o(n).                         \tag{3}
\]

This condition is not supplied for all actual global minimizers.
It is a new full-strength conditional criterion, not an application
of a vanishing-strength or (o(n)\)-relative-entropy assumption.

## 1. Repairing only the missing Gaussian variance

Write (D=D_r\), (V=V_r\), and (q=1-r>0\). Couple

\[
 \widetilde G=G+W,\qquad W\sim N(0,D)\text{ independent of }G,
 \qquad \widetilde\Sigma=\Sigma+D\succeq qI_d.
\]

The diagonal of \(\widetilde\Sigma\) is (1+d_e\), with
\(d_e=D_{ee}\ge0\) and \(\sum_e d_e=V\). Each coordinate pair
\((G_e,\widetilde G_e)\) has correlation \((1+d_e)^{-1/2}\).
The Gaussian angular sign formula therefore gives

\[
 \Pr(\operatorname{sign}G_e\ne\operatorname{sign}\widetilde G_e)
 =\frac1\pi\arctan\sqrt{d_e}
 \le\frac{\sqrt{d_e}}\pi.
\]

Let (L(B)=\log\mathbb E_{x,y}e^{\gamma x^TBy}\), with uniform
spins. Changing one cross sign changes (L\) by at most (2\gamma\).
Thus Cauchy--Schwarz proves the finite coupling bound

\[
 \mathbb E|L(B)-L(\widetilde B)|
 \le\frac{2\gamma}\pi\sum_e\sqrt{d_e}
 \le\frac{2\gamma}\pi\sqrt{dV},
 \qquad\widetilde B=\operatorname{sign}\widetilde G.          \tag{4}
\]

No Gaussian density ratio, invertibility, or independence between
different coordinate pairs is required for this coupling.

## 2. Conditional independent signs with nonconstant total variances

In distribution write
\(\widetilde G=Y+\sqrt q\,z\), where
\(Y\sim N(0,\widetilde\Sigma-qI_d)\) and (z\) is an independent
standard Gaussian vector. Conditional on (Y\), the signs are independent
with means (m_e=2\Phi_{\rm std}(Y_e/\sqrt q)-1\).

The proof in the fixed-strength conditional-noise note uses only these
conditional means and variances. In detail, coordinatewise Taylor
replacement by independent Gaussians of means (m_e\) and variances

\[
 w_e=1-m_e^2
\]

has total error at most
\(C\gamma^3d=O_{c,t}(\sqrt n)\), because
\(|\partial_e^3 L|\le2\gamma^3\) and the centered third absolute
moments of biased signs are at most one. Convexity and global evenness
of (L\) imply that removing the entire deterministic mean matrix
from that centered Gaussian expectation cannot increase its lower
bound. Finally the independent-Gaussian variance derivative lies in
\([0,\gamma^2/2]\). Hence, uniformly in (Y\),

\[
 \mathbb E[L(\widetilde B)\mid Y]
 \ge \mathbb E L(z)-\frac{\gamma^2}{2}\sum_e m_e^2
                                      -O_{c,t}(\sqrt n).       \tag{5}
\]

The latent variance of (Y_e\) is (r+d_e\). Two independent noise
copies conditional on (Y_e\) have total variance (1+d_e\) and
correlation \((r+d_e)/(1+d_e)\). Therefore

\[
 \mathbb E m_e^2
 =\frac2\pi\arcsin\frac{r+d_e}{1+d_e}
 \le\frac2\pi(\arcsin r+d_e).                               \tag{6}
\]

For the last inequality the derivative, for (s\ge0\), is exactly

\[
 \frac d{ds}\arcsin\frac{r+s}{1+s}
 =\frac{\sqrt{1-r}}{(1+s)\sqrt{1+r+2s}}\le1.
\]

Consequently (5), averaged over (Y\), gives

\[
 \mathbb E L(\widetilde B)
 \ge \mathbb E L(z)-\frac{\gamma^2}\pi(d\arcsin r+V)
                                      -O_{c,t}(\sqrt n).       \tag{7}
\]

The already proved Gaussian pure-cross lower bound is
\(\mathbb E L(z)\ge[c\sqrt{2t}K_0-2\log2]n-o_{c,t}(n)\).
Combine it with (4), (7), and the pointwise host-free inequality
\(\min_A F_{A,B}(t)\ge L(B)\). This proves (2).

## 3. Necessary deficit for success, and probability scope

At parameters with (g_r(c,t)>0\), a sequence satisfying the desired
mean comparison with \(2R_n(c/\sqrt n)\) must obey

\[
 \boxed{\quad
 \liminf_{n\to\infty}\frac{V_r}{n}
 \ge \frac{\pi^2 g_r(c,t)^2}{2c^2t}>0.
 \quad}                                                       \tag{8}
\]

Indeed \(2R_n(c/\sqrt n)\le cn+o(n)\). If (8) failed, take a
subsequence along which (V_r/n\) stays bounded strictly below its
displayed threshold. Divide (2) by (n\); its last finite-(V_r\)
term is then (o(1)\), and its square-root loss is strictly smaller
than (g_r(c,t)\), a contradiction. The statement also covers
sequences with unbounded (V_r/n\), by the same subsequence argument.

For the actual canonical covariances, (\|\Sigma\|_{\rm op}=O(n)\)
uniformly by the reviewed bound (4n-3\). If (V_r=o(n)\), one may
also transfer the conditional-noise lower-tail argument to obtain that
the success probability tends to zero: in (5), the function
\(Y\mapsto\sum_e m_e^2\) has a uniformly bounded coordinate derivative
at fixed (q\), covariance operator norm (O(n)\), and mean at most
\((2/\pi)(d\arcsin r+V_r)\). Its order-(d\) upper fluctuations
are exponentially unlikely by Gaussian concentration; the remaining
conditional sign fluctuations have the same bounded-differences bound
as in the fixed-strength note. The repaired law therefore has an
exponentially small lower tail at a level a fixed positive multiple
of (n\) below the floor in (2) without its coupling loss.

Equation (4) and Markov's inequality transfer this only as an (o(1)\)
success probability for the original law. No exponential transfer,
no exclusion of \(\exp(o(n))\) original-law proposals, and no statement
that all actual minimizers satisfy (3) is made. The unrestricted
selected-outcome comparison and original convergence remain open.

No computation, signing census, or simulation was run for this proof.
