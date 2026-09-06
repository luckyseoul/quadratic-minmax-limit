# Root contributing review: actual complete-cross spectral gain

2026-09-05. Full-source analytic PASS; no corrections requested.

## Source and contribution disclosure

Root read the COMPLETE 411-line frozen source
`/tmp/original_mo_complete_cross_flat_spectral_gain.md`, SHA256
`b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`.
The source hash was reverified before preparing this receipt. Root also
read the complete author, exact-worker and independent docs-worker
receipts. Their scope and contributor disclosures match the source.

Root independently checked the proposed contraction formulas, Gaussian
approximation bridge and robust variance/chord extension before they
were integrated. Root is therefore a contributing reviewer of this
proof, not its authorship-independent reviewer. The separate docs
worker had no involvement in this derivation and supplied a complete
independent final-source review.

## Uniform Gaussian-sign lemma

The Schur-product induction gives R^{circ q}<=C I for every positive
integer q. The normalized Hermite coefficients of sign have total
squared mass one and first coefficient sqrt(kappa). These statements
give the claimed dimension-independent L2 bounds, including the tails.

For mixed contractions with r<min(p,q), root expanded the four indices
and verified the trace tr(M N_q M N_p) and its Frobenius-square form.
All Schur powers have positive exponent, and the coefficient magnitudes
1/sqrt(n) give the O(C^4/n) squared norm. For r=p<q, the full
contraction collapses to b=R^{circ p}a and the stated quadratic form,
bounded by C^3/n. Reversal handles r=q<p, including mixed first-chaos
terms. The p=q=r scalar case is removed as the mean variance term.

The Wick product and differentiation identities have the correct
factorials. Root checked the coefficient p(r-1)! binom(p-1,r-1)
binom(q-1,r-1) in Gamma. Its constant terms equal p!||f_p||^2,
so their sum is exactly sigma_Q^2. Every other term has one of the
bounded mixed contractions, proving Var(Gamma_Q)<=K_(Q,C)/n.

Gaussian integration by parts yields the characteristic equation with
the displayed sign and factor t. The integrating-factor kernel is
bounded by one; thus the t^2 sqrt(Var Gamma_Q)/2 estimate is valid
without an inverse-variance assumption. It covers degenerating variances.

At fixed Q, compactness of the variance interval and the characteristic
continuity theorem give subsequential normal convergence. Uniform
second moments imply uniform integrability of absolute values. The
contradicting-sequence argument establishes uniformity over all inputs.
The tail estimate then permits n to infinity first and Q to infinity
second. This proves uniform absolute-moment convergence, not merely
variance control or a growing-dimensional joint column limit.

## Actual sign entries and the spectral penalty

The actual R=BB^T/n has diagonal one. Each actual column gives the
required coefficients of magnitude 1/sqrt(n). The assumption d>=||B||
and m=n/d^2>=m_0 is precisely the bounded-covariance hypothesis.
Values m_0>1 leave an empty family and cause no substantive exception.

Root verified ell_j, t_j and their averages with the normalization nd.
The first Hermite term in the average variance is kappa(1-epsilon).
For each higher odd order the entry exponent is even, making every
summand nonnegative; its diagonal contribution supplies (1-kappa)m.
The exact best-response inequality beta(B)>=E||B^T sign(G)||_1
then uses the uniform separate column marginals with the stated error.

The piecewise lower chord equals t below kappa and is the secant
between kappa and one above kappa. The averaged deficit bound
at most kappa epsilon gives the exact penalty -kappa epsilon.
The gain coefficient is sqrt(kappa)-kappa. Spectral trace inequalities
give 0<=epsilon<=1-m and characterize epsilon=0 as flatness of the
nonzero singular values at d. No finite-order rate is asserted.

## Scope and logical dependencies

The proof is self-contained apart from standard elementary Gaussian,
Hermite and characteristic-function identities, which it derives or
states explicitly. Root did not need the optional multiple-chaos PDF
theorems to validate its logical core and does not claim to have read
those PDFs. The author and independent reviewer record their separate
primary-source statement checks and hashes.

The resulting lower bound belongs to ACTUAL unweighted B. An active
ratio substitution requires the original p=q_A=0 and c=Phi(K), which
imply c=beta(B), together with the stated scalar operator bound.
Neither a trace cap nor near-scalarity alone supplies that operator
bound. Transfer to the full weighted cross law and evaluation over all
actual cells remain separate work. The theorem excludes the flat
kappa endpoint in its stated actual setting, not the convergence problem.

No mathematical computation, checker, scan, numerical integration or
optimization was run for this review. File reads and hashes establish
provenance only. Original convergence and its possible value remain OPEN.
