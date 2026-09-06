# Contributing-author audit: actual complete-cross flat-spectrum gain

2026-09-05. Reviewer: optimized_profile_exact. Complete analytic PASS,
with contribution disclosed. This is NOT an independent full-source review.

## 1. Exact source and review extent

I read all 411 lines of the frozen source end to end:
`/tmp/original_mo_complete_cross_flat_spectral_gain.md`, SHA256
`b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`.
No correction is requested. The verdict covers the complete scalar
Gaussian lemma, robust actual-cross theorem, flat specialization, and
precisely limited active-cell consequence.

Before this review I completely reread the existing 168-line cross
singular-moment source, SHA256
`6d5129a1572842c76c8f11a008b0093cb3c340684a40219b7db8828fdeeaf756`.
That source gives only the earlier cubic floor and does not contain
the present delocalized-sign central-limit gain.

I directly read the optional primary theorem statements, not merely
their titles or summaries: Nualart--Peccati, Theorem 1, PDF page 3;
Noreddine--Nourdin, Theorem 1.1, PDF pages 1--2. Their downloaded
PDF hashes match the source's stated values:

    b8c80cb9b03b016da921316239d0acaed6594fc9303a6176c860fe25db173ea9
    bb07986ac769c63ae02596b470903ba3c4f89c5d83d494c4d3eff2c98d6cb259

The first theorem assumes a normalized fixed chaos of order at least
two; the second permits fixed orders including one and does not require
an invertible limiting covariance. The source accurately describes them
as optional cross-checks. Its logical proof does not import either.

## 2. Uniform Gaussian-sign lemma

For every positive integer q, R^{circ q} is positive semidefinite,
has diagonal one, and is bounded above by C I. The Schur-product
order argument is valid even when R has negative entries or is singular.

The Hermite expansion of sign is centered and odd; its coefficient
of the normalized first Hermite polynomial is sqrt(kappa). Distinct
chaoses remain orthogonal across correlated Gaussian coordinates.
Consequently the tail bound C tau_Q and the full second-moment bound
C are uniform in n, the correlation matrix, and the sign pattern a.

I expanded the mixed contractions directly. When r<min(p,q), their
squared norm is exactly tr(M N_q M N_p), a nonnegative Frobenius
square. The two operator bounds and ||M||_F^2<=C^2/n give C^4/n.
When r=p<q, the collapsed vector b=R^{circ p}a gives the quadratic
form in (3.2), bounded by C^3/n. Reversing factors handles r=q<p.
This includes all interactions with the first chaos. The scalar
p=q=r contraction is the only unbounded-in-this-sense case, and it
is precisely the mean variance term rather than a fluctuation.

The Gaussian polynomial product formula and derivative formula follow
from the displayed generating functions. Applying them to
grad F_Q dot grad U_Q gives exactly the coefficient
p(r-1)! binom(p-1,r-1) binom(q-1,r-1) in (4.4).
The constants are p!||f_p||^2, summing to sigma_Q^2. Symmetrization
does not increase tensor norm. Every nonconstant term is controlled
by (3.1) or (3.2), proving Var(Gamma_Q)<=K_(Q,C)/n with no hidden
dimension dependence for fixed Q.

Gaussian integration by parts yields the exact characteristic-function
equation, including its sign and t factor. Its integrating-factor kernel
is at most one on the integration interval, giving the printed
(t^2/2)sqrt(Var Gamma_Q) error without an inverse-variance bound.

For fixed Q, a compact variance subsequence and the characteristic-
function continuity theorem give distributional convergence to N(0,v),
including v=0. The uniform second-moment bound gives uniform
integrability of absolute values. A contradicting sequence of inputs
would violate this argument, so the absolute-moment convergence is
indeed uniform over all admissible R and a.

The orthogonal tail variance also bounds |sigma-sigma_Q|. Thus the
factor (1+sqrt(kappa))sqrt(C tau_Q) in (5.2) is correct. Taking the
limits in the stated order proves the full sign lemma. No growing-
dimensional joint column limit, independence of chaos components,
or lower bound inferred from variance alone is used.

## 3. Actual columns and the robust flatness penalty

For complete B, the actual row Gram R=BB^T/n has unit diagonal,
and each column supplies coefficients exactly in {+1/sqrt(n),
-1/sqrt(n)}. The actual operator bound gives R<=I/m. These are
the precise hypotheses needed for the uniform lemma; they would
not follow merely by supplying an arbitrary spectral measure.

The full sign covariance is the convergent positive Hermite series.
Its first term gives t_j>=kappa ell_j, while its operator bound gives
t_j<=1. The bound ell_j<=1 and its average 1-epsilon have the
normalizations stated in (6.2).

I checked the averaged variance trace in (6.3): its q=1 term is
exactly kappa(1-epsilon). Every higher odd Hermite order has an even
power of R_ij in that trace, so the n diagonal entries alone give
the retained (1-kappa)m contribution. No off-diagonal cancellation
is assumed to have a favorable sign before the exponent is even.

The uniform marginal absolute-moment estimates may be summed over
columns with error sqrt(m)e_(1/m_0)(n) after division by nd.
This step needs no column independence. The best-response inequality
beta(B)>=E||B^T sign(G)||_1 is an actual Boolean norm inequality.

The corrected chord (7.1) equals t on [0,kappa] and is the ordinary
square-root chord on [kappa,1]. Since the mean positive-part deficit
is at most kappa epsilon, its two epsilon contributions sum to exactly
kappa epsilon. Thus (7.3) proves the announced robust bound with its
printed coefficient, not only an unspecified O(epsilon) loss.

The spectral upper bound and trace Cauchy--Schwarz give
0<=epsilon<=1-m. Equality epsilon=0 forces every nonzero squared
singular value to equal d^2. The flat specialization and its positive
gain when m is bounded below are therefore correct.

## 4. Scope, contribution, and disposition

The conclusion concerns actual unweighted beta(B). Substitution of an
active cross value requires the separately stated actual assumptions
p=q_A=0 and c=Phi(K), which imply c=beta(B). The operator bound on
the scalar scale d remains explicit. It is not inferred from a general
diagonal majorizer, and no transfer across exceptional weighted scales
or evaluation of all original active cells is claimed.

The proof worker originated the flat-cross variance/chord gain and
the delocalized-sign Gaussian-limit bridge. I independently checked
that proposal, supplied the mixed-contraction and characteristic-
equation implementation and the robust epsilon extension, and then
reviewed their complete integration here. The root independently
checked those contributions. Accordingly this is a contributing-author
audit; a genuinely independent final-source review is still separate.

No mathematical computation, checker, optimization, signing search,
simulation, local mathematical job, or remote mathematical job was run
by me. Final review used source reads, PDF text extraction, and SHA256
checks only. No source or canonical repository file was edited.

Verdict: complete analytic PASS for (1.1)--(1.2) and their stated
scope. Original all-orders inequality and convergence remain open.
