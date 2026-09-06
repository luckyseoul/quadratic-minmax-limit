# Independent full-source review: actual complete-cross spectral gain

2026-09-05. Verdict: PASS. No mathematical correction requested.

## 1. Reviewed artifact and independence

I directly read the complete 411-line source
`/tmp/original_mo_complete_cross_flat_spectral_gain.md`, SHA256

    b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20

This review covers its self-contained Gaussianization proof, robust
complete-cross gain, and actual-versus-formal scope. I had no role in
deriving this theorem, its contraction/characteristic-function argument,
or its robust chord extension. My earlier contribution to a different
source/cross nuclear-coupling argument is not a prerequisite or a step
of this source. The present full-source mathematical review is independent
of the proofer, exact worker, and root contributions disclosed in Section 9.

No mathematical computation, solver, signing search, numerical integral,
spectral scan, test, or documentation gate was run for this review. No
source or canonical repository file was edited. Reads, PDF text extraction,
source-page access, and hash checks were the only diagnostic operations;
this review receipt is the sole new artifact written by me.

## 2. Uniform Gaussian-sign absolute-moment lemma

The Schur-product induction is valid for arbitrary real correlation
matrices, including negative correlations and singular matrices:
`(C I-R^{circ q}) circ R` is positive semidefinite and
`(C I) circ R=C I`. Consequently every positive entrywise power has
operator norm at most C. Each power also remains a correlation matrix.

The odd Hermite expansion of sign has squared coefficients summing to
one and first coefficient `sqrt(kappa)`. The correlated Hermite identity
and orthogonality across orders give the stated covariance and tail
formulas even when R is singular. Since `sum_i a_i^2=1`, both retained
and full variances are at most C, and the omitted-chaos variance is at
most `C tau_Q`, uniformly in all admissible inputs.

I checked the complete mixed-contraction expansion, rather than only
the same-order contractions. For `r<min(p,q)`, expanding the squared norm
does give `tr(M N_q M N_p)`, equal to the squared Frobenius norm of
`N_q^(1/2) M N_p^(1/2)`. All powers have positive exponents here.
The coefficient signs cause no problem: `||diag(a)||=1/sqrt(n)` gives
`||M||_F^2<=C^2/n`, and the two remaining operator factors give `C^4/n`.

For `r=p<q`, the collapsed coefficient is exactly
`a_j (R^{circ p}a)_j`. Its tensor norm squared is bounded by
`(C/n)||R^{circ p}a||^2<=C^3/n`. Reversing the roles treats `r=q<p`.
In particular this proves the mixed first/higher-chaos estimates. The
only unbounded-by-these-estimates contraction, `r=p=q`, is scalar and
becomes the constant variance term. Symmetrization is norm decreasing.

## 3. Finite-chaos approximation and tail passage

The Wick generating-function argument gives the product identity with
the stated factorials. Polarization extends it from tensor powers to
general symmetric tensors. Its derivative identity and Gaussian
integration by parts give `-L U_Q=F_Q` and
`E[F_Q h(F_Q)]=E[h'(F_Q) Gamma_Q]` without an independence assumption.

The coefficient in the displayed Gamma expansion is correctly

    p (r-1)! binom(p-1,r-1) binom(q-1,r-1).

For `p=q=r` this is `p!`, so those terms sum precisely to
`sigma_Q^2`, not to an additional random component. Every remaining
term is controlled by one of the mixed estimates. There are finitely
many for fixed Q; polynomial isometry and a finite triangle inequality
therefore prove `Var(Gamma_Q)<=K_(Q,C)/n`. Cross-covariances between
different displayed terms do not need to vanish for this bound.

Differentiating the characteristic function gives the exact sign and
normalization

    phi_Q'(t)+sigma_Q^2 t phi_Q(t)
      =-t E[(Gamma_Q-sigma_Q^2) exp(it F_Q)].

For nonnegative t, the integrating-factor kernel is at most one on
`0<=s<=t`, which yields `t^2 sqrt(Var Gamma_Q)/2`. Negative t follows
by conjugation. This estimate uses no positive lower variance bound,
so it remains valid for a zero or degenerating variance.

For fixed Q, every potentially nonuniform sequence has a subsequence
whose variances converge in the compact interval `[0,C]`. The displayed
characteristic-function estimate then forces convergence to the normal
law with that variance, including the degenerate law at zero. The
uniform second-moment bound makes the absolute values uniformly
integrable, with tail expectation at most `C/M`. This proves convergence
of absolute first moments and, by the same subsequence contradiction,
uniformity over R and the coefficient signs at each n.

The final truncation step is also uniform: orthogonality gives
`|sigma-sigma_Q|<=sqrt(C tau_Q)`, and the absolute-value Lipschitz
inequality controls the change of the first absolute moment by the
same quantity. Taking n to infinity first with Q fixed, then Q to
infinity, proves (2.1). A variance lower bound alone is never substituted
for this Gaussianization argument. No growing-dimensional joint limit
of the n cross columns is asserted or needed.

## 4. Actual complete-column normalization and retained gain

Every row and column of the actual sign matrix has squared norm n.
Thus `R=BB^T/n` has diagonal one, and `d>=||B||op` gives
`R<=I/m`. Feasibility implies `m<=1`; cases with the stipulated
`m_0>1` are simply vacuous. The covariance series of `X=sign(G)`
has the claimed positive-semidefinite upper bound `I/m`.

For each actual column, the definitions of ell_j and t_j give

    0<=ell_j<=1,   average ell_j=1-epsilon,
    kappa ell_j<=t_j<=1.

I checked both upper bounds with `||b_j||^2=n` and `m d^2=n`.
The average of the t_j is exactly `tr(E C_X)/n`, with `E=mR`.
Its first-chaos summand is `kappa(1-epsilon)`. For every higher odd
order q, the exponent q+1 is even, so all entry summands are nonnegative;
the n diagonal entries already contribute m. This proves the additional
`(1-kappa)m` variance term. Passing the infinite series is justified by
its nonnegative scalar summands, as well as the established L2 limit.

For each j the coefficient vector of `h_j/sqrt(n)` has the required
exact magnitudes `1/sqrt(n)`. Therefore the uniform lemma supplies an
error at most `sqrt(n)e_(1/m_0)(n)` in E|h_j|. Summing these marginal
errors and dividing by nd gives exactly `sqrt(m)e_(1/m_0)(n)`.
The deterministic inequality `beta(B)>=E||B^T X||_1` is valid for this
random sign vector; no independence between column sums is used.

For `b=sqrt(kappa)/(1+sqrt(kappa))`, the robust inequality (7.1)
has right side t below kappa and the secant through kappa and one
above kappa. Both inequalities have the correct direction. Moreover
`(kappa-t_j)_+<=kappa(1-ell_j)`. Averaging gives the penalty
`-[b+(1-b)]kappa epsilon=-kappa epsilon`, while
`b(1-kappa)=sqrt(kappa)-kappa`. This verifies precisely

    beta(B)/(nd)
      >=kappa+(sqrt(kappa)-kappa)m-kappa epsilon-o_(m_0)(1).

The trace inequalities establish `0<=epsilon<=1-m`. Equality
`epsilon=0` says that the sum of the nonnegative terms
`lambda(d^2-lambda)` is zero for the eigenvalues lambda of `B^T B`.
It is therefore equivalent to every nonzero singular value being d.
Since `sqrt(kappa)>kappa`, exact or asymptotic flatness and fixed
positive lower m force a strict leading-order gain over kappa.
The source appropriately gives no finite-n rate for its final error.

## 5. Active-cell interpretation and remaining scope

For clarity, I checked the established definitions of the paired matrix,
Phi(K), and original p,q,c in the defining excerpts of
`original_mo_diagonal_majorizer_cross_covariance.md` and
`original_mo_diagonal_majorizer_weighted_shell_upper.md`.
No result from those notes is used in Sections 2--7.

For any paired signing, flipping the right signs changes its cross
energy from c to -c while fixing the internal energy. Thus the maximum
absolute full energy is at least |c| for every cross state, and
`beta(B)<=Phi(K)`. At a stated actual state with `p=q_A=0` and
`c=Phi(K)`, the reverse inequality follows from the definition of beta,
so `c=beta(B)` as claimed. At scalar scale d, `u=c/(nd)` then has
the bound in (8.1), provided the actual operator bound and lower m hold.

The formal flat cross law with `u=kappa` is excluded from the specified
actual scalar, operator-bounded pure-cross setting. This is an
actual-entry restriction beyond the previously retained cubic cross
inequality, not a claim that the complete formal relaxation is realized.
The source does not derive the required actual unweighted operator
bound from near-scalarity or trace control of an arbitrary diagonal
majorizer. It makes no transfer to all weighted cross laws, no assertion
that every active cell is pure-cross, no all-cell ellipsoid evaluation,
and no closure of the original inequality. These exclusions are material
and correctly stated.

## 6. Optional primary-source checks

I read Theorem 1 on page 3 of the supplied Nualart--Peccati PDF and
Theorem 1.1 on pages 1--2 of the supplied Noreddine--Nourdin PDF, and
opened their arXiv source pages. The first theorem is the normalized
fixed-chaos contraction/fourth-moment/Gaussian-limit equivalence for
order at least two. The second allows fixed orders at least one and
a general Gaussian covariance matrix; its theorem statement does not
require invertibility. The later quantitative bound that does use an
inverse covariance is a different theorem and is not being applied.

The independently checked local PDF hashes match Section 9:

    original_mo_flat_cross_clt_primary_nualart_peccati.pdf
    b8c80cb9b03b016da921316239d0acaed6594fc9303a6176c860fe25db173ea9

    original_mo_flat_cross_clt_primary_noreddine_nourdin.pdf
    bb07986ac769c63ae02596b470903ba3c4f89c5d83d494c4d3eff2c98d6cb259

Their source links are https://arxiv.org/pdf/math/0503598 and
https://arxiv.org/pdf/1009.1310 . The descriptions in the reviewed
source are accurate. These theorem checks are optional corroboration;
the self-contained finite-chaos and tail proof above supplies the
logical bridge actually used by this artifact.

Final disposition: the complete reviewed source passes independent
analytic and scope review at the stated hash, with no requested changes.
