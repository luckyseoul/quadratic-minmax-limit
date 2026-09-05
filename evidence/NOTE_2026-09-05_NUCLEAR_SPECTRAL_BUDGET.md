# A nuclear spectral budget for complete signing minima

2026-09-05. **All-orders spectral inequality; no spectral flatness or
cross-order transport conclusion.** No computation was run.

For a complete symmetric zero-diagonal signing A of order n>=2, put

\[
 Q_A(x)=\tfrac12x^TAx,\qquad
 \Phi(A)=\max_{x\in\{-1,1\}^n}|Q_A(x)|,
 \qquad |A|=(A^2)^{1/2},
\]
\[
 d_i=|A|_{ii},\qquad D=\operatorname{tr}|A|=\sum_j|\lambda_j(A)|.
\]

Then

\[
 \boxed{\quad
 \Phi(A)\ge\frac2\pi\sum_{i<j}\frac1{\sqrt{d_i d_j}}
       \ge\frac{n^2(n-1)}{\pi D}.
 \quad}                                                     \tag{1}
\]

For the actual half-product and symmetric pressures

\[
 a_A(\beta)=\tfrac12\left(\log\mathbb E_x e^{\beta Q_A(x)}
                         +\log\mathbb E_x e^{-\beta Q_A(x)}\right),
 \qquad F_A(\beta)=\log\mathbb E_x\cosh(\beta Q_A(x)),
\]

one has, for every beta>0,

\[
 \boxed{\quad
 F_A(\beta)\ge a_A(\beta)
 \ge\frac{2\beta}\pi\sum_{i<j}\frac1{\sqrt{d_i d_j}}-n\log2
 \ge\frac{\beta n^2(n-1)}{\pi D}-n\log2.
 \quad}                                                     \tag{2}
\]

## Proof of the finite inequalities

Both matrices |A|+A and |A|-A are positive semidefinite and have the
same diagonal d. Each d_i is positive: a zero diagonal entry of the PSD
matrix |A| would make its ith row zero, contradicting
`(A^2)_ii=n-1>0`.

Let g^+ and g^- be centered Gaussian vectors of covariances |A|+A and
|A|-A, respectively, and let x^+=sign(g^+), x^-=sign(g^-). No joint
coupling between these two vectors is needed. Their coordinate-pair
correlations are

\[
 r_{ij}^{\pm}=\frac{|A|_{ij}\pm A_{ij}}{\sqrt{d_i d_j}}\in[-1,1].
\]

The elementary planar Gaussian sign identity gives
`E x_i^+ x_j^+=(2/pi)arcsin(r_ij^+)`, and likewise for the minus law.
Because arcsine has derivative at least one on (-1,1), with its
continuous endpoint extension, and A_ij is a sign,

\[
 A_{ij}\bigl(\arcsin r_{ij}^+-\arcsin r_{ij}^-\bigr)
 \ge\frac2{\sqrt{d_i d_j}}.
\]

Consequently

\[
 \mathbb E Q_A(x^+)-\mathbb E Q_A(x^-)
 \ge\frac4\pi\sum_{i<j}\frac1{\sqrt{d_i d_j}}.              \tag{3}
\]

The left side is at most 2 Phi(A). This proves the first inequality
in (1). With k=binom(n,2), AM--GM and Cauchy--Schwarz give

\[
 \sum_{i<j}\frac1{\sqrt{d_i d_j}}
 \ge2\sum_{i<j}\frac1{d_i+d_j}
 \ge\frac{2k^2}{(n-1)D}
 =\frac{n^2(n-1)}{2D},                                     \tag{4}
\]

proving its second inequality.

Each of the two Gaussian-sign laws is a probability measure on the
2^n spin states, so its relative entropy from the uniform spin law is
at most n log 2, including when its Gaussian covariance is singular.
The finite Gibbs variational principle applied separately at +beta and
-beta therefore gives

\[
 a_A(\beta)\ge\frac\beta2
       \bigl(\mathbb E Q_A(x^+)-\mathbb E Q_A(x^-)\bigr)-n\log2.
\]

Combine (3)--(4). Finally `F_A>=a_A` is the arithmetic--geometric mean
inequality for the two one-sided partition functions. This proves (2).

## Consequences for the actual minima

For every complete signing, `tr A^2=n(n-1)` and
`D<=n sqrt(n-1)` by Cauchy--Schwarz. Thus (1) includes the existing
uniform lower bound `Phi(A)>=n sqrt(n-1)/pi`, and strengthens it when
the nuclear norm is below its maximum.

If A_n is an exact global norm minimizer, the reviewed conference
construction upper bound `m_n<=(1/2+o(1))n^(3/2)` and (1) give

\[
 \frac{D}{n^{3/2}}\ge\frac2\pi-o(1),\qquad
 \frac{D^2}{\operatorname{tr}A^2}
       \ge\left(\frac4{\pi^2}-o(1)\right)n.                 \tag{5}
\]

For a fixed c>0, if A_n minimizes either a_A(c/sqrt(n)) or
F_A(c/sqrt(n)), compare with the same admissible conference-based norm
construction. The optimized pressure is at most `(c/2+o(1))n`.
Equation (2), applied to the actual pressure minimizer, consequently
gives

\[
 \frac{D}{n^{3/2}}
 \ge\frac{c}{\pi(c/2+\log2)}-o_c(1),\qquad
 \frac{D^2}{\operatorname{tr}A^2}
 \ge\left(\frac{c^2}{\pi^2(c/2+\log2)^2}-o_c(1)\right)n.
                                                               \tag{6}
\]

The quantity D^2/tr(A^2) is defined here as the nuclear effective rank;
it is at most the ordinary matrix rank. The proof uses only the stated
objective cap, not a new global-minimizer variation. In particular,
the finite inequalities remain valid for every complete signing.

These bounds enforce a linear-size bulk in this effective-rank sense.
They do not imply `||A||_op=O(sqrt(n))`, rule out isolated eigenvalues
of size n^(3/4), or establish a small or large fixed-threshold canonical
tensor deficit `tr[-rI-H/mu]_+` without additional hypotheses. They do
not compare different optimized orders or prove the original limit.

The proof reuses only the elementary Gaussian angular sign identity
already proved in CORE.md, the finite Gibbs variational identity, and
CORE.md's reviewed conference upper bound. No numerical spectral or
Gaussian calculation is used.
