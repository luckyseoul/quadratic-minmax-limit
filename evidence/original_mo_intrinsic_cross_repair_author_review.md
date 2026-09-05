# Author receipt: intrinsic cross normalization and rank-four shell repair

2026-09-05. Author: optimized_profile_docs_gate.
This is an AUTHOR receipt, not an independent review of my own proof.

## Frozen artifact and complete-read provenance

Source: /tmp/original_mo_mu_joint_shell_extension.md, 338 lines.
SHA-256:
1dcd9b1e76b00887e406e505113c854b80f0661bb3bd69283f6486fb59fa2d53.

I authored and checked the complete proof, equations (1)--(21).
The source has remained byte-identical since its first complete draft.
Both optimized_profile_proof and root independently read all 338 lines
at this hash and returned mathematical PASS with no corrections.
Those two independent reads, not this author receipt, provide the
independent-review provenance.

## Author checks and exact result

The intrinsic covariance is normalized by
mu=max(2,||A tensor A-S_B+I||op). Its diagonal is one, its operator
norm is at most two, and mu>=n-1. No whole-source spectral-square
bound is assumed for this normalization.

The unrepaired joint-shell matrix M is explicitly allowed to be
indefinite. The four-rank PSD addition P and the exact sum-of-squares
identity (5) prove M+P is PSD directly from R>=0. The independent
cushion is added with its correct vn I covariance. The proof makes
no invalid Gaussian decomposition involving an indefinite covariance.

The full exchange increment identity (6) retains its mixed terms.
Adding P preserves increment domination, so (7) is a genuine upper
for the actual cross Gaussian process on each nonempty joint shell.

The nonuniform marginal bound (10) uses only marginal variances
at least a_0, the exact jointly selected Hamming mismatch, and
concavity of square root for marginal inflation. It does not assume
independent field coordinates or independent flip selection.

Under the original norm caps, the normalized arcsine calculation
controls internal reference errors by O(L_A sqrt(n)) and the cross
reference error by O(n). The proof includes diagonal normalization,
the signed linear trace estimate, and the nonlinear Frobenius
remainder. The bound L_A^2<=8 Phi(A) is justified by real cube
polarization and real-to-complex interpolation.

The ratios k/(mu a_0)<=1/(n ell) make all constants uniform as the
actual threshold tends to either infinity and the noise vanishes.
The exact dual mismatch changes by O(n^(-1/2)); the possible
square-root endpoint loss gives the final O(n^(5/4)) shell-width
error in (19). The actual Gaussian maximum (21) includes both phases
and the explicit polynomial-shell concentration remainder.

The result preserves the leading opposite-cross mismatch at the
intrinsic mu denominator. It does not evaluate the shell upper
against 2 sqrt(2) Phi(A), reverse the conditional floor, assume full
order-2n optimality, or prove original all-orders convergence.

## Publication scope

This receipt concerns only the frozen mathematical source identified
above. No repository source or tests were changed, and no unchanged
mathematical computation or experiment was rerun for this receipt.
Canonical destination paths and cross-reference aliases may be
recorded in the import manifest without changing these reviewed bytes.

