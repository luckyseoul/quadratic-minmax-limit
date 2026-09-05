# Independent full-source proof review: source/cross trace boundary

Date: 2026-09-05. Reviewer: optimized-profile proof worker.

Verdict: PASS for the finite actual source/cross nuclear inequality and
the explicitly stated formal-relaxation obstruction. This is not a
complete-signing counterexample and does not close the original target.

## Exact source and independence

I read the complete final 444-line source
`/tmp/original_mo_source_cross_nuclear_trace_boundary.md`, SHA256
`106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`.
This review covers its whole argument, not just the previous 370-line
version or the final all-threshold delta. I did not author its finite
transfer, AM-GM refinement, formal model, or all-threshold extension.
The docs-gate reviewer contributed the AM-GM refinement; the root
authored Section 9; neither contribution is treated here as independent
authorship-free review. This receipt supplies a separate full-source
mathematical review of the combined final artifact.

The previously fully read original nuclear, near-scalar normalization,
cross singular-moment, and pure-cross actual-measure prerequisites have
the hashes recorded in the source. I also rechecked the threshold
coefficients directly in Section 4 of
`original_mo_diagonal_majorizer_cross_covariance.md`, hash
`0b3921d43d88424457ad2ed777ee158e8ac34c6751c995f0c3b86aee870e95ff`,
and the original drift and completion-square constant in Sections 1
and 4 of `original_mo_diagonal_majorizer_weighted_shell_upper.md`, hash
`9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.

## Finite actual inequality

The Frobenius comparison adds only nonnegative diagonal summands.
Jensen applies to the positive numbers t_i through t^(-1/2), so the
cross term has the claimed sign. The trace-norm transfer bound is
conservative but valid for each n-dimensional internal block.

The exact complete cross sign squares give m=m_0 ell h, with
ell+h=2(1+delta); AM-GM therefore gives the lower bound for sqrt(m_0)
in the needed direction. The complete source nuclear bound has the
factor (1-1/n) and normalization 2alpha exactly as displayed.

T is the actual weighted contraction. Its block row squares imply
A_L^2<=I-WW^T and A_R^2<=I-W^TW. Operator-monotonicity of the square
root legitimately gives the two nuclear trace upper bounds without
assuming that either internal block commutes with W. Taking traces
finishes (1.1), with the same actual cross measure on both sides.

The contraction row-square inequality and Cauchy--Schwarz imply
dbar^2>=N-1. The source nuclear lower bound and the Frobenius upper
bound on tr|A| give the stated positive lower bound on alpha. These
facts make the asymptotic error uniform without assuming a trace cap.
No use is made of a false contraction claim for K/dbar.

## Formal retained data and their exact scope

I checked the spectral laws, the commuting block realization, all
three first absolute moments, r=2m, u=f sqrt(m), and positive gamma.
The balanced internal eigenvalues give zero diagonal after Hadamard
conjugation; constant squared diagonals follow from flat basis entries.
For irrational m the limiting/weighted-trace qualification is necessary
and is present. The text does not assert complete-signing entries.

Both full moment inequalities, the source nuclear and common-zero-odd-
diagonal cubic inequalities, the cross cubic inequality, and the new
source/cross inequality are satisfied. In particular the source cubic
normalization really gives the same kappa/[2sqrt(1-m)] bound as the
source nuclear normalization. The cross cubic inequality is retained
and is saturated at u=kappa; it has not been omitted from the model.

The formal alpha, f, and u are not assigned to Boolean norms or active
states of the displayed real matrices. Thus further restrictions on
those particular matrices cannot be silently read into the specified
trace relaxation. Conversely the formal example supplies no actual
complete-signing counterexample and no realized active optimizer.

## Centered and shifted functional checks

I checked the endpoint-law substitution in both terms of (6.1).
The same nu is used in both expectations; no separate envelope is
substituted. The Minkowski lower bound and the affine-ratio endpoint
minimum have the correct directions. The original target is
sqrt(2)alpha, not f/2. The rational squared margin is 501/125000.
The previously verified pi enclosure is reused, not recomputed.

For shifted thresholds, the square-to-disk Gaussian probability bound
proves a_h<=kappa. This leaves s in [0,kappa^2], while u and the
Boolean completion-square constant stay fixed at kappa. Factoring out
sqrt(w) from the covariance is legitimate. The derivative of the
combined C(s) is negative throughout this interval; the proof does
not incorrectly reuse the centered lower bound on A_* separately.

Negative metrics increase both endpoint-law terms: the unit-atom first
numerator difference is 2(kappa-s)>0, and the second difference is
4st y/(1-t^2 y)^2>=0. The endpoint limits exist after cancellation.
The drift is the ORIGINAL absolute drift z f/2. Adding it to the
noise lower bound gives (9.2) by sqrt(1-z^2)>=1-z. Infinite-threshold
limits are included and cannot evade the positive target gap.

## Execution and stopping point

No numerical calculation, solver, signing search, new rational checker,
or metric optimization was run for this review. It is an analytic
full-source audit; file reads and hashes are provenance checks only.
The actual finite inequality survives independently of the formal
obstruction. The missing implication is actual entry/frame/Boolean
information beyond the explicitly retained spectral/block constraints.
No stronger impossibility theorem or original-inequality closure is
asserted by this PASS.
