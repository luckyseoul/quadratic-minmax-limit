# Corollary: actual-Gibbs covariance tails and a quadratic cross-block floor

2026-09-05. A short corollary of the reviewed field-response proof,
not a radial upper comparison. No optimality or computation is used.

Let J be a real symmetric zero-diagonal n by n interaction with
`Phi(J)<=b n`, and let U be its actual zero-field Ising covariance.
Define

\[
 \kappa=\log(1+\sqrt2),\qquad
 B_b=1+\frac{4\pi b}{\kappa},\qquad
 \chi_b=\frac12\exp\left[-2\sqrt{8B_b(1+4B_b)}\right]>0.
                                                               \tag{1}
\]

The completely proved augmentation in
`evidence/NOTE_2026-09-05_NORM_CAP_FIELD_RESPONSE.md`
(SHA256 `46f6465c9a889dc485b9c24dac6f7fef8849d27271cc86df11b94ab732ed52dd`)
supplies an actual joint spin/Gaussian law. Set

\[
 R=\operatorname{diag}(r_i),\qquad
 r_i=\mathbb E\operatorname{sech}^2(g_i).
\]

Total covariance gives

\[
 U=R+\operatorname{Cov}(\tanh g)\succeq R,\qquad 0\preceq R\preceq I.
                                                               \tag{2}
\]

Section 4 of that proof, at zero field, bounds the average posterior
second moment on at least n/2 coordinates by `8B_b(1+4B_b)`.
Jensen/Cauchy--Schwarz and `sech^2 u>=exp(-2|u|)` give

\[
 \boxed{\qquad \operatorname{tr}R\ge\chi_b n.\qquad}          \tag{3}
\]

All expectations are under the actual augmented posterior, not its
Gaussian prior. For every rank-r orthogonal projection P, (2)--(3) imply

\[
 \operatorname{tr}((I-P)U)
 \ge\operatorname{tr}R-\operatorname{tr}(PR)
 \ge\chi_b n-r.                                             \tag{4}
\]

Thus every rank `r=o(n)` covariance truncation has an extensive tail,
also for the covariance at -J and for the average of the two phases.

Let U and V be the actual covariances at J and -J, respectively, and
use their diagonal lower bounds R_+ and R_-. For any real n by n
matrix C, successive PSD trace comparisons give

\[
 \operatorname{tr}(C^T U C V)
 \ge\operatorname{tr}(C^T R_+ C R_-)
 =\sum_{i,j} C_{ij}^2(r_+)_i(r_-)_j.                         \tag{6}
\]

For integral full cross blocks, (6) is at least
`(tr R_+)(tr R_-)>=chi_b^2 n^2`. Exchanging U,V gives

\[
 \boxed{\quad
 \min_{C\in\{-1,1\}^{n\times n}}
 \frac12\left[\operatorname{tr}(C^T U C V)
              +\operatorname{tr}(C^T V C U)\right]
 \ge\chi_b^2 n^2.
 \quad}                                                       \tag{7}
\]

If instead C=A is a complete symmetric zero-diagonal signing with
J=beta A, then its diagonal entries vanish. Equation (6) gives

\[
 \boxed{\quad
 T_A=\operatorname{tr}(A U A V)
 \ge\chi_b^2 n^2-\sum_i(r_+)_i(r_-)_i
 \ge\chi_b^2 n^2-n.
 \quad}                                                       \tag{8}
\]

At fixed `beta=c/sqrt(n)`, a Boolean cap `Phi(A)<=C_0 n^(3/2)`
gives fixed `b=cC_0`. The conditional sublinear-rank/sublinear-tail case
in Section 5 of the integral covariance-rounding note is therefore
unavailable under this cap. Independent signs give the universal upper
bound `min_C qbar(C)<=n^2`, so this minimum has quadratic order here.
The desired radial upper bound `2a_A'(beta)/beta+o(n^2)` can also be
quadratic, and chi_b may be tiny: it is not refuted by this corollary.
No zero-temperature-uniform or integrated order comparison is claimed.
