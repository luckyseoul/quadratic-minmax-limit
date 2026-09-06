# Author checkpoint: actual complete-cross spectral gain

2026-09-05. Author: optimized-profile proof worker.

Frozen source: `/tmp/original_mo_complete_cross_flat_spectral_gain.md`,
complete 411 lines, SHA256
`b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`.

I reread the complete written source after freezing it and rechecked
the derivation end to end. Author verdict: PASS for its actual-cross
theorem and the scope explicitly stated there. This is an author
checkpoint, not an independent review or a statement of publication.

## New actual theorem

For complete sign B, d>=||B||op, m=n/d^2>=m_0>0, and
epsilon=1-tr[(B^T B)^2]/(n^2 d^2), the theorem gives uniformly

    beta(B)/(nd)>=kappa+(sqrt(kappa)-kappa)m
                                      -kappa epsilon-o_(m_0)(1).

Here kappa=2/pi. Exact flat nonzero singular values at d give
epsilon=0, hence a positive leading-order gain above kappa.
The asymptotic error is uniform but has no asserted finite-n rate.

The substantive bridge is a uniform scalar absolute-moment Gaussian
limit for sums n^(-1/2)sum_i sign_i sign(G_i), with unit marginal
variances and bounded covariance operator. It is proved directly,
not inferred from a variance lower bound alone. The argument checks
all mixed Hermite contractions, including order-one interactions;
subtracts the constant variance terms in Gamma; controls its remaining
variance by O_(Q,C)(1/n); and solves the characteristic equation.
Variance compactness, second-moment uniform integrability, and a
uniform L2 Hermite tail bound establish the required absolute moment.
Only uniform separate marginals are used for the n columns.

The literal complete entries are used twice: every column has the
delocalized coefficients needed by that lemma, and the Gaussian
correlation matrix has unit diagonal. Higher Hermite covariance terms
retain positive diagonal contributions in the averaged variance.
The final piecewise chord converts these into the displayed gain
and its exact kappa epsilon penalty.

## Provenance and contribution limits

The exact worker contributed the mixed-contraction/characteristic-
equation route and the robust epsilon extension. Root independently
checked those contributions and requested integration. The source
discloses those facts. Their full-source checks are contributing-author
audits and are not counted as the independent review.

The docs-gate worker had no involvement in this derivation and is
performing the independent full-source audit of this frozen artifact.
Its receipt is a separate artifact; this author checkpoint does not
pretend to supply or supersede that review.

The two primary multiple-chaos theorem statements were read directly
from the downloaded PDFs and their hashes verified. They are optional
cross-checks, not logical prerequisites of the self-contained proof.
No mathematical computation, new checker, solver, scan, signing
search, or numerical experiment was run. Downloading source PDFs,
reading files, and hashing artifacts were provenance operations only.

## Remaining boundary

The operator bound is on the ACTUAL unweighted B, at scale sqrt(n).
It is not inferred from arbitrary near-scalar diagonal majorizers.
The gain transfers to a pure-cross ratio u only when an actual active
state gives c=beta(B)=Phi(K) and the scale d satisfies the stated
hypotheses. It excludes the kappa-floor endpoint in that restricted
actual setting, not every weighted trace limit or every active cell.
The original inequality and full weighted/frame transfer remain open.
