# Author receipt: full-SDP gap and two actual original quadratic phases

2026-09-05. Author: optimized_profile_exact.

## Frozen source and direct reread

The final source is all 274 lines of
`/tmp/original_mo_full_sdp_gap_original_phase_bound.md`, SHA256
`1d36878bdd157be36b1e935f0e92a0e977cbbabb1bbf23784a645860ac1142c0`.
I read that entire frozen source directly after writing it and again
after receipt of its independent reviews. I confirm the finite inequality,
its original-norm normalization, its asymptotic interpretation, and its
explicit restrictions. No source correction was requested by the reviewers.

## Independently reviewed objects

I directly read the entire proof-worker review at
`/tmp/original_mo_full_sdp_gap_original_phase_bound_proof_review.md`, SHA256
`18bc090624c2453cc76ee500669ce7a1b3aa54441896dc971ce499117895755f`.
It independently checks the full source and records analytic PASS.

I also directly read the entire 129-line independent docs-gate review at
`/tmp/original_mo_full_sdp_gap_original_phase_bound_docs_review.md`, SHA256
`e52ad4b2811640495793639b3e19510e1e2a86af22594e308b4321bf37d21474`.
It independently checks the full source and records analytic PASS.

The proof-worker supplied the residual/compatibility prerequisite and
discussed the desired small-gap consequence and the need to avoid a
maximum-coordinate phase normalizer. That participation is disclosed in
its receipt. The docs-gate reviewer records no theorem-development
participation. Neither reviewer authored or edited this phase construction
or its masked-residual proof. The final source hash above is unchanged.

## Confirmed result

For an actual complete symmetric zero-diagonal signing K of order N>=3,
an actual trace-optimal same-diagonal majorizer D, S=tr D, q=N-1,
L=||K||op, kappa=2/pi, and rho=1-kappa, define

    gamma=(S-tr|K|^3/q)/S,
    eta=S/[N sqrt(q)],
    c_*=1-1/sqrt(2),
    b_*=min(1/2, eta(1+N/q+c_*)sqrt(gamma)).

The finite original quadratic-norm bound is

    Phi(K)/S >= (kappa/2)[1-gamma-b_*-sqrt(2gamma b_*)]
                           -rho(L/q)(1-gamma).

The two PSD sign-phase correlations use the same coordinatewise
normalization v_i=q+|(K|K|)_ii|. Their original expected energies
lie in [-Phi(K),Phi(K)], yielding the exact coefficient kappa/2
after cancellation. No rectangular-objective transfer supplies that
coefficient. Only off-diagonal sign squares are used in the normalization
loss. The complete weighted residual bound, arithmetic-harmonic
dispersion, and mask m_ij<=1/2 give the retained b_* and mixed loss.

Under a fixed original norm cap Phi(K)<=C N^(3/2), the consequence is

    Phi(K) >= (kappa/2)S
                -O_C(N^(3/2)sqrt(gamma)+N^(5/4)).

Thus gamma tending to zero gives the actual source-scale conclusion
Phi(K)>=(kappa/2-o_C(1))S. Small gamma remains an additional hypothesis;
it is not inferred from a norm cap or from optimizer status.

## Paired use and remaining open ranges

For the literal paired source of order N=2n, general original shells
retain the ratio c/Phi(K). The conclusion c_D/n>=kappa-o_C(1) needs
the separate active original assumptions p=q_A=0 and c=Phi(K), with
the positive sign chosen. The source does not evaluate the resulting
Gaussian field width or maximize over the other original cells.

The displayed finite lower estimate is provably vacuous for every
gamma>=1/4: then b_*=1/2 and its main bracket is at most -1/4.
This is a limitation of this particular fully evaluated formula, not
an impossibility theorem for all optimized original-phase arguments.
The complementary positive-gap sources remain unresolved. Original
MO convergence is not claimed.

## Execution status

This package is analytic. I ran no mathematical computation, numerical
phase evaluation, signing search, census, optimization, SDP, simulation,
or test, on the coordinator or on any offload host. Tool use for this
package was limited to reading, writing, listing, and hashing documents.
There is no new machine-result artifact or execution claim.
