# Independent new-source review: cross-endpoint to source nuclear transfer

2026-09-06. Reviewer: optimized_profile_docs_gate.

## Frozen source and complete reads

I directly read the ENTIRE 230-line source
`/tmp/original_mo_near_scalar_cross_endpoint_source_nuclear_transfer.md`,
SHA256
`6a486df0fd46aa76259e3f02e3734eb2529162500f98f89af58e90562e6a2187`.
I also directly read its complete 80-line author receipt,
`/tmp/original_mo_near_scalar_cross_endpoint_source_nuclear_transfer_author_receipt.md`,
SHA256
`345c96871c15a798d23803af4c00bb618607e6a770fe65823a2fde7ab3ed11c9`.
Both hashes and line counts were checked directly.

Both named prerequisites were directly read completely earlier in this
same analytic review sequence, and their unchanged hashes were refreshed
for this task:

- `original_mo_near_scalar_diagonal_spectral_normalization.md`, 280 lines,
  `c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`;
- `original_mo_near_scalar_internal_flat_law_transfer.md`, 141 lines,
  `f65ce2200fd926ba969c9bc5bbaf8ecec8a79b8d228e0f17865fc56c9d9775a8`.

The required good-coordinate and congruence mechanisms are restated in
the new source. Its conclusion does not import the earlier 141-line
theorem's prescribed internal spectral law.

I supplied no proposal, derivation, correction, or proof step to this
new 230-line transfer. I previously contributed an elementary directional
AM-GM scale comparison to a different older source/cross note and authored
a separate local-update scalar lemma. Those historical roles are disclosed;
this receipt does not claim independence for every older mechanism in the
research chain. I independently checked all steps of the newly written
transfer, including its stronger inverse-mean lower, full internal second
moment recovery, and noncommuting nuclear estimate. The separate scalar
and prospective all-law gain arguments are not premises here.

## 1. Actual inputs and exact scale bounds

The theorem assumes an actual complete paired signing, a positive
DIAGONAL D feasible for both signs, dispersion tending to zero, and
weak endpoint-law convergence for ALL n actual squared singular values
of W, including zeros. It imposes no internal spectral shape.

Feasibility makes the full weighted matrix a contraction and therefore
W a contraction. Its squared singular values have common support [0,1],
so endpoint-law convergence implies m_D=tr(WW^T)/n tends to m.

Literal cross sign squares give m_D=m_0 ell h exactly. The normalization
of t gives ell+h=2(1+delta) and L+R=2. The two directional
Cauchy--Schwarz bounds imply 1/ell<=L and 1/h<=R. Hence

    (ell+h)/(ell h)<=2,
    ell h >= (ell+h)/2 = 1+delta.

The upper ell h<=(1+delta)^2 is ordinary AM-GM. Both quantities
are positive, so division gives precisely (2.3), with its printed
directions. Since delta tends to zero and m_D tends to positive m,
m_0 tends to m and dbar/sqrt(n) tends to 1/sqrt(m).
Thus the nonzero endpoint mass supplies the scale without a trace cap.

## 2. Common original labels and full internal second moment

The dispersion identity is exact. Each excluded index label accounts
for at least one bad coordinate among the two halves, giving the factor
two in (3.1). With epsilon=delta^(1/3), the removed fraction tends
to zero. The explicit delta=0 convention covers mixed zero/nonzero
dispersion sequences. In particular q/n tends to one and q grows.

The same literal A_J is used on both halves. Principal feasibility
gives the contraction bounds for H and H_J. Good-coordinate congruence
gives A_J/dbar=Q H_J Q, the actual retained-source operator bound,
and the at-most-3epsilon ordered-eigenvalue comparison. It does not
give an operator bound on the full untrimmed A/dbar.

For the second moments, unnormalized positive and negative squared
eigenvalue sums each lose between zero and n-q under compression.
Thus M_H-a M_J is in [0,2b]; since M_J is in [0,1],
M_H-M_J is in [-b,2b], proving the asserted absolute 2b estimate.
This accounts for the distinct n and q normalizations.

The square function has Lipschitz constant at most 3 on the common
eventual spectral interval. Weyl comparison contributes 9epsilon,
and completeness gives tr(A_J^2)=q(q-1) exactly. This proves (3.3).
The scale already recovered from the FULL cross moment then gives
M_H->m. No assumed internal-law convergence was used to obtain it.

Taking traces in the ACTUAL block inequality H^2+WW^T<=I yields
M_H+m_D<=1. Both moments tend to m, hence m<=1/2. This restriction
is derived, not an unstated extra premise on the endpoint measure.

## 3. Noncommuting full internal nuclear bound

For fixed tau in (0,1), the threshold 1-tau carries no limiting
endpoint-law mass, so the rank fraction p_tau tends to m.
Trace pairing of H^2+WW^T<=I with the PSD projection P_tau
gives the first inequality of (4.1). The second follows from the
definition of this spectral projection of WW^T. It does not require
H to commute with W, WW^T, or P_tau.

The decomposition H=H(I-P_tau)+H P_tau is exact. Its summands
need not be symmetric; the nuclear-norm triangle inequality and the
rank--Frobenius bound apply to general matrices. In detail,

    rank H(I-P_tau)<=n-r_tau,
    ||H(I-P_tau)||_F^2=tr[H^2(I-P_tau)]<=tr H^2,
    rank H P_tau<=r_tau,
    ||H P_tau||_F^2=tr(H^2 P_tau)<=tau r_tau.

After division by n, these give exactly
sqrt[(1-p_tau)M_H]+p_tau sqrt(tau), including the second coefficient.
The first term deliberately retains the FULL internal second moment;
no unjustified block diagonalization or internal spectral atom is used.

First taking n to infinity at fixed tau gives the upper
sqrt[m(1-m)]+m sqrt(tau). Taking tau to zero afterward proves
the second assertion of (1.4). There is no moving-cutoff rate claim.
This is solely a matrix nuclear-moment estimate, not a replacement
of an original Boolean objective by a rectangular norm.

## 4. Original-source nuclear and quadratic comparisons

Unnormalized nuclear norm decreases under principal compression of a
symmetric matrix. Combining tr|H_J|<=tr|H| with the 1-Lipschitz
absolute-eigenvalue comparison gives the exact inequality in Section 5.
Multiplication by dbar/sqrt(q)->1/sqrt(m), with a tending to one,
transfers the full weighted bound sqrt[m(1-m)] to the original retained
source cap sqrt(1-m).

Unbiased independent extension of a Boolean state on J preserves its
expected original quadratic energy and proves Phi(A_J)<=Phi(A).
Feasibility gives Phi(A)<=tr(D_L)/2<=n dbar, so the recovered scale
bounds the full normalized Phi(A). Consequently the factor a^(-3/2)
in (5.1) differs from one by an admissible o(1) error.

For m>=9/25, the retained source has limsup operator ratio at most
1/sqrt(m)<=5/3 and limsup nuclear ratio at most sqrt(1-m)<=4/5.
This proves both caps without the earlier internal flat-law hypothesis.

## Verdict and stopping scope

PASS for the complete frozen 230-line source, with no required correction.
Every conclusion in (1.2)--(1.4) follows from the stated actual full
cross-endpoint law and near-scalar feasible diagonal assumptions.
The common original principal source is used only for source moments
and the original Phi comparison. The original paired covariance,
cross block, active field, and optimizer are not replaced.

The theorem itself provides no additional source-norm lower bound or
profile exclusion. A subsequent gain theorem using the two source caps
still requires its separate complete proof and review. Arbitrary active
sequences are not shown to have small dispersion or an endpoint cross
law; the all-profile implication and global original MO target stay OPEN.

No mathematical program, checker, solver, numerical evaluation,
construction, or search was run. Tools were used only for complete
source reads, line counts and hashes, and this /tmp review. I made no
canonical edit and performed no publication or backup operation.
