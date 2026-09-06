# Full-source review: strict original-source near-flat gain

2026-09-06. Reviewer: optimized_profile_exact.

## 1. Frozen source, full reads, and precise independence

I directly read the ENTIRE final 612-line source
`/tmp/original_mo_original_source_near_flat_strict_gain.md`, SHA256
`7726b89e1c39429cde75ff887b981cbd3cf831adb17b04f20193a3c6dbb35298`.
I previously read its complete initial 611-line version. The author's
sole final change explicitly identifies D as diagonal in Section 10.
I reread the entire final source, not only that changed passage.

I fully refreshed and directly read all five named prerequisites in the
present review sequence, and checked each exact hash:

- `original_mo_complete_cross_flat_spectral_gain.md`, 411 lines,
  `b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`.
- `original_mo_original_source_local_update_scalar_gain.md`, 209 lines,
  `7de99c4bbf997fc25eafa2742cb55c220dc13fdf29d0b1ae535358ea8c73f155`.
- `original_mo_near_scalar_internal_flat_law_transfer.md`, 141 lines,
  `f65ce2200fd926ba969c9bc5bbaf8ecec8a79b8d228e0f17865fc56c9d9775a8`.
- `original_mo_original_phase_spectral_moment.md`, 262 lines,
  `7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`.
- `original_mo_source_cross_nuclear_trace_boundary.md`, 444 lines,
  `106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`.

These names refer to their frozen /tmp copies. The scalar proof also has
my separate complete independent 104-line review, SHA256
`8f01834cc576419c61f51ac790426e95b3ce1f4997bc390658f1bd42381f56e0`,
at `/tmp/original_mo_original_source_local_update_scalar_gain_exact_review.md`.

I did not propose, derive, correct, or supply any NEW step of the present
Sections 1--9. I independently checked that new source-gain chain after
the full written source became available. However, I supplied the older
411-line prerequisite's mixed-contraction/characteristic-equation route
and robust cross extension, and authored the separate 141-line transfer
used in Section 10. I also contributed to the older phase/trace sources.
This review therefore does NOT claim independence for that entire
prerequisite provenance chain or for my own transfer theorem. It records
an independent check of the new robust-source construction, its new
distinguished-coordinate extension, and the new application of the
separately disclosed prerequisites.

## 2. Spectral projections and actual correlation normalization

The operator cap supplies a common compact spectral interval. The fixed
cutoffs at plus or minus lambda/2 avoid all limiting atoms. The bounded
piecewise continuous squared-error functions in Section 2 therefore have
empirical integrals tending to zero, giving both Frobenius estimates and
the two rank limits. No exact finite-n kernel or spectral symmetry is
silently assumed.

Completeness gives A_ii=0 and (A^2)_ii=n-1. Taking diagonals in those
Frobenius estimates yields mean-square convergence T_ii to zero and
P_ii to rho. Since 2P_+=P+T, the raw R0=(2/rho)P_+ has diagonals
converging to one in mean square, not necessarily uniformly.

The printed eta choice indeed tends to zero and has h_n/eta^2 tending
to zero, even if h_n=0 on some orders. Hence only o(n) coordinates are
bad. On good coordinates the positive diagonal normalization is defined;
the independent identity block on bad coordinates preserves an actual
PSD unit-diagonal covariance on ALL n original coordinates. Its operator
norm is at most 2/rho+o(1).

Removing o(n) rows/columns of bounded-operator R0 costs o(sqrt(n))
in Frobenius norm, as does inserting the identity there. Good diagonal
normalization costs O(eta sqrt(n)). Thus R-R0=o_F(sqrt(n)) is proved
with all exceptional coordinates accounted for. No part of A is deleted.

## 3. Schur powers, baseline, and full higher-chaos mean

The positive Schur multiplier argument correctly gives
R^{circ q}<=||R||op I for every positive integer q, including negative
entries of R. Since tr R=n, tr R^2<=C_R n. The Hermite tail covariance
series converges in operator norm, with total coefficient mass 1-kappa.

For every odd q>=3, actual unit-magnitude off-diagonal entries and zero
diagonal give |tr(A R^{circ q})|<=tr R^2=O(n), uniformly in q.
Summing this bound controls the entire original-energy nonlinear tail.
The R-R0 trace pairing is o(n^(3/2)); the positive spectral rank and
flat-atom limits give tr(AR0)=n lambda+o(n^(3/2)). Consequently the
actual positive phase baseline tends to 5kappa/8 with the correct
Q_A=x^TAx/2 normalization.

The identity

    A^2/lambda^2-rho R+A/lambda=o_F(sqrt(n))

is valid by the two projector approximations and P+T=rho R0. Pairing
it with each odd Schur power gives the claimed trace lower: the leading
term is at least rho n because every R_ij^(q+1) is nonnegative and all
diagonal terms equal one. The complete-entry term is O(sqrt(n)), and
the Frobenius error is o(n), UNIFORMLY over all odd q>=3.

Summation with the nonnegative Hermite coefficients is therefore valid
without an uncontrolled infinite-series error. The exact identity
lambda^2 rho=n yields mu_n>=1-kappa-o(1), with no missing factor of
two. Each actual row has squared length n-1, so the tail covariance
operator bound also gives v_i<=2(1-kappa)/rho+o(1) uniformly. This is
the strong mean and heterogeneous variance cap actually needed later.

## 4. First-chaos alignment and row parameters

The average squared first-chaos alignment error equals exactly

    (kappa/n^2)tr[(A-lambda I)^2 R].

For R0 the trace is (2/rho)||(A-lambda I)P_+||_F^2=o(n^2).
The remaining trace is o(n^2) by
||(A-lambda I)^2||_F=O(n^(3/2)) and R-R0=o_F(sqrt(n)). This uses
the actual operator cap and does not infer covariance control from a
matrix-product Frobenius bound or invert the possibly singular R.

Higher chaoses are orthogonal to every first-chaos coordinate, giving
c_i=E[G_i F_(i,1)] and sigma_i^2=b_i+v_i. Row lengths and the bounded
covariance operator supply uniform second-moment bounds. Cauchy--Schwarz
then proves both empirical L1 convergences b_i to a^2 and c_i to a.
In particular all but o(n) rows have c_i>=a/2.

## 5. The new joint marginal Gaussian lemma

I checked the generalized coefficient bounds, not only the previous
all-nonzero coefficient case. For M=diag(d)R^{circ r}diag(d),
||M||_F<=C/sqrt(n) follows from max|d_j|<=1/sqrt(n) and ||d||<=1.
Thus the exact mixed-contraction trace formula gives C^4/n. If one
tensor order is fully contracted into a larger order, the displayed
z=R^{circ p}d formula gives C^3/n. Equal-order full contractions are
precisely constant covariances and are not erroneously bounded as small.

The added distinguished-coordinate contraction has squared norm

    (d circ R_i)^T R^{circ(q-1)}(d circ R_i)
       <=(C/n)sum_j R_ij^2<=C^2/n,

because R^2<=C R and R_ii=1. Together these estimates cover every
nonconstant contraction for sG_i+tF_Q, including its first/higher-chaos
terms. The first/first contribution is a constant variance term.

The finite Gaussian integration-by-parts and product identities are
used with the correct Hermite normalization. They yield
Var Gamma_Q<=K_(Q,C,s,t)/n. The characteristic-function equation has
the correct sign and coefficient, and its integrating-factor kernel
has modulus at most one on the relevant integration interval. Therefore
the error in (6.7) is at most z^2 sqrt(K/n)/2 even for zero or vanishing
variance; no inverse-variance bound has been hidden.

The omitted Hermite tail is L2-bounded by C tau_Q and orthogonal both
to G_i and to the retained chaoses. Taking n first and then Q gives
the Gaussian characteristic approximation for every fixed linear
combination. Compact covariance subsequences and Cramer--Wold prove
the required TWO-dimensional joint limit. Uniform second moments give
uniform integrability of |F|, and a contradiction subsequence proves
uniformity over all admissible inputs.

For c>=c_*>0 both limiting marginals have positive variance. The two
coordinate axes thus have zero limiting Gaussian probability, even
when the pair is degenerate along a line. This validates convergence
of the sign-disagreement probability, including correlation one.
The angular formula arccos(c/sigma)/pi has the correct normalization.
This proof does not require a growing-dimensional CLT, independence
of different local fields, or any external chaos-limit theorem.

## 6. Actual gain, disagreement, and quadratic update

Each actual row has coefficients d_j=A_ij/sqrt(n), with its one zero
diagonal entry. It meets both generalized coefficient conditions.
Uniform absolute-moment Gaussianization, row parameter convergence,
and the baseline therefore give precisely (7.1).

On the all-but-o(n) rows with c_i>=a/2, the joint lemma applies
uniformly. The discarded fraction costs o(1) in mismatch frequency.
The Gaussian angle is uniformly continuous on the compact domain
where its marginal variance is bounded below. Restricting first to
small row-parameter errors and then shrinking the tolerance correctly
handles its non-Lipschitz correlation-one endpoint. At b=a^2, c=a the
angle is arctan(sqrt(v)/a), including v=0. Thus (7.2) follows without
assuming a positive lower bound or homogeneity of v_i.

The Boolean update X'=X+xi circ Delta remains on the same n coordinates.
Conditional independent Bernoulli masks have distinct-coordinate product
mean epsilon^2; the diagonal terms vanish because A_ii=0. Expanding
Q_A therefore gives the exact linear epsilon term and quadratic
epsilon^2/2 term in (8.1), not an independence assumption on the spins.

Each mismatch contributes exactly four to ||Delta||^2, so the lower
penalty after normalization is 2epsilon^2 C_n p_n, with
C_n=||A||op/sqrt(n). The cap 5/3+o(1), not the limiting atom 5/4,
is retained. Taking the expectation of actual original Boolean energies
then proves (8.2); there is no rectangular-norm or polarization factor.

## 7. Scalar conclusion, transfer application, and verdict

The separate scalar prerequisite was fully and independently checked.
Section 9 correctly extends its chord bound to v_i<=V+o(1), with an
o(1) average error on a fixed compact interval. The probability bound
and fixed epsilon=1/10 yield

    liminf[epsilon G_n-2epsilon^2 C_n p_n]>=16/3125.

The strong mean gives liminf sqrt(mu_n)>3/5; monotonicity of 3r^2-r
on that range justifies the displayed endpoint lower. Adding the
baseline proves the first bound of (1.3). The rational comparisons
16/3125>1/200 and 1/200-1/440=3/1100, together with
5kappa/8>35/88, prove its stated strict inequality.

Section 10 accurately applies the separate same-original-source transfer
to A_J with q/n tending to one. That lemma supplies the actual cap 5/3,
the required empirical law, and the normalized Phi upper comparison;
therefore the strict lower transfers back to the original A. The explicit
diagonal hypothesis on D is present in the final source. No diagonal
optimality or separate trace cap is required. I checked this application,
but my authorship of the transfer itself is disclosed above.

Verdict: PASS for the ENTIRE final 612-line source, with no required
correction. It proves the stated bounded near-flat actual source theorem
and excludes its specified actual near-scalar internal-law regime at
alpha tending to 2/5. It does not force arbitrary extremizing sequences
into that regime, exclude all other profiles, close all active cells,
or prove the global original MO target. The scope limitations in
Sections 1, 10, and 11 are correct and must be retained.

No mathematical computation, solver, checker, numerical integration,
optimization, signing construction, or search was run on any host.
Tools were used only for complete reads, hashes/line counts, and writing
this review in /tmp. No canonical repository file was edited and no
publication was performed during this review task.
