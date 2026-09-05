# Author receipt: near-scalar actual full-spectrum normalization

2026-09-05. Author: optimized_profile_exact.

## Frozen source and complete direct reads

The final source is all 280 lines of
`/tmp/original_mo_near_scalar_diagonal_spectral_normalization.md`, SHA256
`c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`.
I directly read that complete frozen source after its final pre-review
edits. I confirm the exact finite phase and nuclear bounds, the full
weighted moment transfer, the asymptotic error, and the scope limits.

I also directly read the entire 175-line independent review at
`/tmp/original_mo_near_scalar_diagonal_spectral_normalization_docs_review.md`,
SHA256
`5e97d723cbcb05628ee2f6bbba591250ec14e8fbc37e784e0471b2e3fb13c1f2`.
It records full analytic PASS, including complete reads of the 280-line
source and its 262-, 303-, and 252-line prerequisites. No mathematical
or documentary correction was requested. The frozen source is unchanged.
Root separately reported its complete 280-line read and checked the
finite phase/nuclear constants, interlacing, and fixed-cap asymptotics.

## Development and review provenance

Root posed the bounded problem of separating near-scalar diagonal
dispersion from a small canonical gap, and required that any trimming
remain an auxiliary original-norm LOWER bound. I derived and authored
the new moment-transfer and normalization argument. The independent
docs-gate reviewer had no theorem-development, source-writing, or
pre-derivation involvement, as recorded in its receipt.

I previously authored the metric-stability prerequisite and independently
reviewed the gap compatibility and original phase/nuclear prerequisites.
Those dependencies have their own independent review records. This
receipt does not advertise my author reread as an independent review.

## Confirmed theorem and constants

For an actual complete signing K of order N>=4 and ANY positive
diagonal D with D+-K>=0, let S=tr D, dbar=S/N,
delta=S tr(D^(-1))/N^2-1, r=(N-1)/dbar^2, and
T=D^(-1/2)KD^(-1/2). Let mu_j be the empirical average of
|lambda|^j over the FULL actual weighted spectrum of T.

Under the separate cap S<=C N^(3/2), the proved normalization is

    2Phi(K)/S >= max{kappa(1+r)mu_3/(2r), kappa r/mu_1}
                               -O_C(delta^(1/3)+N^(-1/2)),
    mu_2=r+O(delta^(1/3)).

No optimality assumption on D or small-gap assumption is required.
The proof removes at most bN auxiliary coordinates, with
b<=delta(1+epsilon)/epsilon^2, and compares the principal signing's
moments with the full actual T by interlacing loss at most 2b and
congruence moment errors C_j epsilon, C_j=3,9,81/4.

The exact finite phase and nuclear bounds are (4.3) and (5.2) of the
source. The phase's common variance and retained fourth-moment term
are explicit; its error is O_C(N^(-1/2)). The nuclear estimate uses
an UPPER estimate for its denominator. Taking epsilon=delta^(1/3)
gives the asserted uniform error. The scalar delta=0 case needs no
trimming. The alternative fixed-epsilon then epsilon-to-zero limit
is also controlled by the finite estimates.

Only in the exactly scalar case does the source use
mu_3=r(1-gamma), giving the positive-gap phase consequence
kappa(1+r)(1-gamma)/2-O_C(N^(-1/2)) for 2Phi(K)/S.
The equality is not inferred for near-scalar diagonals with exceptional
coordinates. The actual nuclear moment remains present.

## Source preservation and open endpoint

The auxiliary principal signing only proves Phi(K)>=Phi(K_I), by
random extension of an extremizing state. It never replaces the
original K, its covariance, its weighted full spectral measure, or its
actual cross matrix W_D in the upper-bound construction.

For the original paired source and the explicit active final-cell
conditions p=q_A=0 and positive c=Phi(K), the delta-only uniform
energy estimate transfers the theorem to the actual representative
u_D=c_D/n, at extra cost at most 2sqrt(delta). General cells retain
the ratio c/Phi(K). Small delta is not inferred from source nearminimality
or conditional optimality, and arbitrary scalar enlargement is not
advertised as a useful sharp-scale substitute.

The actual full/cross spectral coupling, all-cell field-width evaluation,
and desired original comparison remain open. This theorem is a new
normalization for the near-scalar positive-gap branch, not its closure
or a proof of original MO convergence.

## Execution status

No mathematical computation, signing construction, numerical evaluation,
simulation, optimization, solver, census, or test was run by the author
for this package, locally or remotely. Document operations only were
used. No machine-result artifact is claimed.
