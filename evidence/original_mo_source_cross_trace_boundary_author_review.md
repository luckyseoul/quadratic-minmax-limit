# Author review: actual source/cross coupling and trace-relaxation boundary

2026-09-05. Final-source end-to-end author review; not an independent review.

## 1. Frozen source and review extent

Source: `/tmp/original_mo_source_cross_nuclear_trace_boundary.md`

SHA256:
`106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`

The author read all 444 lines of this exact final source end to end after
the all-threshold extension and scope changes were integrated. Verdict:
PASS for the stated finite theorem and explicitly delimited formal
relaxation calculation. No original-problem closure is asserted.

The threshold covariance source was refreshed in full, all 384 lines:
`original_mo_diagonal_majorizer_cross_covariance.md`, SHA256
`0b3921d43d88424457ad2ed777ee158e8ac34c6751c995f0c3b86aee870e95ff`.
The weighted-shell source was refreshed in full, all 381 lines:
`original_mo_diagonal_majorizer_weighted_shell_upper.md`, SHA256
`9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.
Their threshold coefficients and original drift were checked directly,
not inferred from a centered-law analogy.

The nuclear, near-scalar, cross-moment, and pure-cross prerequisites
are separately identified by final source hashes inside the proof.
They had been completely read during the source's preparation and
previous bounded reviews. The finite transfer and formal calculations
in the new source are themselves given explicitly.

## 2. Finite actual-signing theorem checked

The literal off-diagonal sign squares and Jensen's inequality give the
full Frobenius comparison (2.1). The trace-norm comparison on each
internal block has precisely the error in (2.2); no operator-norm
approximation of the whole unweighted matrix is assumed.

The original source nuclear lower bound supplies (2.3). The exact
cross second-moment identity and AM-GM give the directional estimate
sqrt(m_0)>=sqrt(m)/(1+delta), not its reverse. The two actual block
inequalities from T^2<=I, followed by operator monotonicity of square
root, give (2.5) without requiring commutation. These facts prove
(1.1) for every feasible positive D, without trace optimality or a cap.

The row-square estimate dbar^2>=N-1, the universal source lower bound,
and m<=1 make the error O(sqrt(delta)+1/n) absolute for delta<=1.
The matrices, original energies, final cells, and actual cross measure
remain unchanged throughout. There is no trimmed covariance or
replacement of the original source by a contractive scalar proxy.

## 3. Formal relaxation and all-metric calculation checked

The chosen mu, nu, and chi laws have the stated first three moments.
They are mutually consistent with the explicit block contraction built
from a projection P and balanced internal eigenvalues on its complement.
Rational-rank Hadamard-basis approximations realize only the stated
real-matrix diagonal and trace data, not complete-signing magnitudes.

The full nuclear and cubic inequalities, source nuclear and cubic
inequalities, cross cubic inequality, and new source/cross inequality
all hold at the displayed formal parameters. Their satisfaction does
not create actual Boolean norms, active states, or optimal frames.

For the centered law, the two expectations use the same endpoint law.
Minkowski's inequality and the exact affine-ratio endpoint minimum give
(7.1). The coarse rational comparisons prove the uniform squared lower
bound 40501/125000, exceeding the actual target squared 8/25 by
501/125000. Negative metrics have the stated termwise comparison.

## 4. All-threshold extension independently checked by the author

The Gaussian square/disk containment yields w>=exp(-h^2), hence
k/w<=kappa. It does not require a numerical Mills-ratio approximation.
For spectral s in [0,kappa^2], C(s) is strictly decreasing, so its
centered-law lower bound is valid throughout. The individual centered
bound A_*>7/10 is deliberately not extended to smaller s.

The negative-metric unit-atom difference is 2(kappa-s)>0. The Gaussian
noise contribution scales by sqrt(w), while the Boolean completion
constant sqrt(kappa) remains unchanged. On the formal active face the
original drift divided by 2n^(3/2) is exactly zf/2. Its combination with
the noise lower bound gives (9.2) for all thresholds and signed metrics,
including all continuous endpoint limits.

The target is sqrt(2)alpha, not f/2. The latter appears only as the
normalized original drift at infinite threshold. Every stopping claim
remains restricted to the explicitly retained formal relaxation.

## 5. Authorship, computations, and status

The exact worker authored the finite Frobenius/nuclear transfer and the
formal trace-relaxation example and centered-law comparison. The docs
reviewer independently checked the transfer and contributed its AM-GM
refinement before the note was written. The root authored the
all-threshold extension; the exact worker independently checked that
extension and incorporated it in the final source. This author receipt
is not a claim of independent authorship or independent review.

No mathematical computation, checker, optimization, signing search,
metric scan, simulation, local job, or remote job was run for this note.
The already verified rational enclosure of pi was reused; its original
checker and result hashes appear in Section 7 and were not rerun.
File reads, scoped edits, line counts, and SHA256 hashes are the only
execution used in finalizing this source and receipt.

Actual source/cross coupling: proved. Explicit formal trace/block
relaxation: insufficient for the stated all-threshold ellipsoid
certificate. Actual complete-signing realizability and active Boolean
frame compatibility: not proved. Original MO inequality: still open.
