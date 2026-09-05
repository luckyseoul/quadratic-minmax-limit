# Independent proof review: fixed-internal-block cross regularization

2026-09-05. Reviewed by the proof agent.

Final source: `/tmp/original_mo_cross_only_operator_regularization.md`.
The entire 264-line source was read. Its SHA-256 is
`27d9ab77768e8b7afa2d48d041cf3fe6bf3b66e8b16e481ca12abcf906a28d4f`.

Result: **PASS within the stated scope.**

The proof correctly obtains one diagonal SDP majorizer for both signs
of the bipartite dilation, thresholds its diagonal to bound the retained
cross rectangle, and chooses one filler satisfying both the operator and
Boolean-norm estimates. The two displayed failure probabilities sum to
less than one already at n=2. No independence of the success events is
used.

The two blockwise phase flips delete the selected cross strips but also
delete the two internal boundary cuts. Restoring those cuts incurs exactly
the two explicit beta terms in (3). This essential correction is retained:
the argument never assumes monotonicity of cross-strip deletion with the
internal blocks fixed. Cauchy--Schwarz then gives the printed
`(2+sqrt(2) K_A) sqrt(Lambda C/K)` bound.

The near-conditional corollary correctly assumes `epsilon_B>0`; its
regularized cross block is not asserted to be an exact minimizer. The
source-regularization corollary keeps the source perturbation and
`epsilon_A^(-4) epsilon_B^(-2)` cross-cap dependence explicit. The final
section correctly leaves arbitrary exact-source regularization and
absorption of cap-dependent errors unresolved.

The imported finite tensor-rounding inequality and same-order source
regularization statements were checked directly against their cited
repository sources during the prior full review. The final source differs
from that reviewed 264-line version only by adding `For epsilon_B>0` in
Section 5. Reversing exactly that text change reproduces its SHA-256
`fcb20c448f022caf374c43ae1fca610e62164862a6d63d75b25c305bf8c34a71`.

This receipt certifies the theorem and its carefully limited corollaries;
it does not certify the still-open sharp Gaussian upper comparison or
original convergence.
