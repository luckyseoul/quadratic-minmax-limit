# Independent review: near-scalar internal flat-law transfer

2026-09-06. Reviewer: optimized_profile_docs_gate.

## Frozen sources and complete direct reads

I directly read the complete 141-line source
`/tmp/original_mo_near_scalar_internal_flat_law_transfer.md`, SHA256
`f65ce2200fd926ba969c9bc5bbaf8ecec8a79b8d228e0f17865fc56c9d9775a8`.

I also directly read its complete 50-line author receipt,
`/tmp/original_mo_near_scalar_internal_flat_law_transfer_author_receipt.md`,
SHA256
`3cd3e469d7669c72922746f8238fd9f819fd8964a2792d8bc5e28e73af112631`,
and the complete named 280-line prerequisite
`/tmp/original_mo_near_scalar_diagonal_spectral_normalization.md`, SHA256
`c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`.
All three hashes were checked directly in this review task.

The new source uses the prerequisite's good-coordinate, interlacing,
and congruence mechanism. Those facts are also derived in the new source
and checked below; its conclusion does not import the prerequisite's
separate Gaussian phase, nuclear, or active-field inequalities.

I did not derive, propose, amend, or supply any step of this 141-line
transfer theorem. I previously contributed to a separate prospective
local-update scalar argument, but that argument is not a premise of this
lemma. This receipt is therefore an independent review of the complete
new transfer theorem, not of the separate prospective gain theorem.

## 1. Hypotheses and common original labels

The matrix K is an actual complete signing of the literal paired form
with the same internal A on both diagonal halves, with opposite signs.
The diagonal D is positive and feasible for both signs of K. Optimality,
a small canonical-primal gap, and a separately imposed trace cap are not
required hypotheses. Dispersion tends to zero and the FULL actual left
weighted internal law tends to the specified chi_*.

With t_i=d_i/dbar, the mean of t_i is one and the mean of 1/t_i is
1+delta. Expansion of (t_i-1)^2/t_i gives exactly the stated dispersion
identity. Outside [1-epsilon,1+epsilon], its lower bound is
epsilon^2/(1+epsilon): on the lower side the denominator at the endpoint
is smaller, giving an even stronger bound.

The number of removed original index labels is at most the number of
bad coordinates among all 2n coordinates. Consequently
b=(n-q)/n <= 2delta(1+epsilon)/epsilon^2. The choice
epsilon=delta^(1/3) makes b tend to zero. When delta=0, every t_i=1,
and taking all labels and epsilon=0 proves the stated zero-dispersion
convention, including sequences alternating between the two cases.
In particular q/n tends to one and q tends to infinity.

Intersecting the two half-good sets is essential and is actually done:
the same literal principal source A_J is available in the left and right
halves. No separately chosen internal source is inserted into the pair.

## 2. Feasibility and congruence bounds

Principal feasibility gives D_L plus or minus A positive semidefinite,
so H_L and its principal compression H_J are contractions. On J,
A_J/dbar=Q H_J Q exactly, where Q has entries sqrt(t_i).
The good interval implies ||Q|| <= sqrt(1+epsilon) and
||Q-I|| <= epsilon. Hence

    ||A_J|| <= (1+epsilon)dbar,
    ||Q H_J Q-H_J||
       <= epsilon(sqrt(1+epsilon)+1) <= 3epsilon.

This bounds only the retained actual complete principal source. There
is no unsupported operator bound on the full untrimmed A/dbar.

## 3. Full-law transfer and exact scale recovery

Deleting n-q coordinates from a symmetric matrix gives the stated
unnormalized counting-function inequalities. After normalization,
a F_J <= F_L <= a F_J+b, whence the uniform CDF difference is at
most b. Thus the empirical law of H_J has the same weak limit as that
of the FULL H_L.

The finite normalized second-moment estimate is also valid. Writing
m_L=tr(H_L^2)/n and m_J=tr(H_J^2)/q, the positive and negative squared
parts each lose between zero and n-q under principal compression.
Thus 0 <= m_L-a m_J <= 2b. Since 0 <= m_J <= 1,

    -b <= m_L-m_J <= 2b,

so |m_L-m_J| <= 2b as claimed. This explicitly accounts for the
different n and q normalizations rather than treating them as equal.

Weyl comparison pairs the ordered eigenvalues of H_J and A_J/dbar
within 3epsilon. Both supports lie in
[-(1+epsilon),1+epsilon]. Eventually epsilon<=1/8, so x^2 is
Lipschitz with constant at most 3 on that common interval. The normalized
second moments therefore differ by at most 9epsilon. Completeness and
the zero diagonal give tr(A_J^2)=q(q-1) exactly, proving (3.2).

The full H_L supports lie in [-1,1], so the assumed weak convergence
also gives second-moment convergence. The limit is
(16/25)(3/4)^2=9/25. Accordingly (q-1)/dbar^2 tends to 9/25;
since q tends to infinity, dbar/sqrt(q) tends to 5/3.

The same congruence comparison transfers the entire empirical law,
not just its second moment. Rescaling by the now convergent positive
factor dbar/sqrt(q) sends the nonzero atoms plus or minus 3/4 to
plus or minus 5/4 and retains the zero mass 9/25. The previously
proved norm bound gives limsup ||A_J||/sqrt(q) <= 5/3.

Thus the specified nonzero full internal second moment supplies the
source scale without an external trace cap. No equality of the operator
norm limit to the empirical nonzero atom 5/4 is asserted or inferred;
possible spectral outliers are covered by the weaker 5/3 cap.

## 4. Retaining the original quadratic norm

Fixing any Boolean state on J and extending outside J by independent
unbiased signs leaves its expected quadratic energy equal to its A_J
energy. All crossing and removed internal edges have zero expectation.
The original maximum absolute energy therefore dominates each such
fixed energy in absolute value, proving Phi(A_J)<=Phi(A) exactly.

After normalization this gives
Phi(A_J)/q^(3/2) <= a^(-3/2)alpha_n. To replace the right side by
alpha_n+o(1), the source correctly proves alpha_n is bounded:
feasibility gives Phi(A)<=tr(D_L)/2<=S/2=n dbar, while the established
scale and a tending to one give dbar/sqrt(n) tending to 5/3.
This completes every conclusion of (1.1).

## Verdict and scope

PASS, with no required correction. The theorem yields one large complete
principal submatrix of the SAME original A, its stated empirical flat
law, a bounded normalized operator norm, and an original quadratic norm
upper comparison. The source expressly keeps the original paired field,
cross covariance, and full weighted input law unchanged.

A later independent strict source lower bound could combine with this
lemma to exclude its particular actual near-scalar internal-law regime.
That lower bound is neither proved nor assumed in this transfer lemma.
No current profile exclusion, all-profile result, all-cell width bound,
or original MO closure follows from this receipt alone.

No mathematical calculation program, solver, construction, numerical
evaluation, or search was run. I used tools only for complete source
reads, line counts and hashes, and writing this review in /tmp. No
canonical repository file was edited and no publication was performed.
