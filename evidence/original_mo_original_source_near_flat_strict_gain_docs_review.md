# Complete contributing review: original-source near-flat strict gain

2026-09-06. Reviewer: optimized_profile_docs_gate.

## Frozen source, direct reads, and contribution disclosure

I directly read the complete initial 611-line source and then the
complete final 612-line source:

`/tmp/original_mo_original_source_near_flat_strict_gain.md`

SHA256
`7726b89e1c39429cde75ff887b981cbd3cf831adb17b04f20193a3c6dbb35298`.

The final change explicitly makes D diagonal in Section 10. I checked
the final hash and line count and reviewed the entire final version,
not merely that change or a selected excerpt. I also fully read the
111-line author receipt, SHA256
`dc9fd2ef2b174c8c994a6475d26787e902aa7ca67e3ab8617eabbdc461f03f9f`.

All five named supporting sources were directly read in full during
this review sequence, with their hashes checked against Section 11:

- `original_mo_complete_cross_flat_spectral_gain.md`, 411 lines,
  `b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`;
- `original_mo_original_source_local_update_scalar_gain.md`, 209 lines,
  `7de99c4bbf997fc25eafa2742cb55c220dc13fdf29d0b1ae535358ea8c73f155`;
- `original_mo_near_scalar_internal_flat_law_transfer.md`, 141 lines,
  `f65ce2200fd926ba969c9bc5bbaf8ecec8a79b8d228e0f17865fc56c9d9775a8`;
- `original_mo_original_phase_spectral_moment.md`, 262 lines,
  `7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`;
- `original_mo_source_cross_nuclear_trace_boundary.md`, 444 lines,
  `106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`.

The complete 141-line transfer and its 280-line normalization
prerequisite were read and independently reviewed immediately before
this source-gain task. That separate transfer review is
`original_mo_near_scalar_internal_flat_law_transfer_docs_review.md`,
146 lines, SHA256
`bbade6ad8cd6539970f8e73075a22fe460d3e4b06fb7ceb13b5d24417f500465`.

This is a CONTRIBUTING review of the combined result. I authored the
distribution-free scalar support used in Section 9, with root supplying
the corrected operator cap and the live fixed-probability specialization.
I also previously supplied an AM-GM step in the older 444-line source;
that step is unrelated to its reused pi enclosure. Neither role is
misrepresented as independent review of my own derivation.

I supplied no derivation or amendment to the new projector construction,
stronger variance mean, first-chaos alignment, joint Gaussianization,
or actual Boolean-update identity. The checks of those links below are
independent of their derivation, but do not make the combined package
a wholly independent review. I requested no source correction.

## 1. Actual hypotheses and robust spectral projectors

The theorem assumes an actual complete symmetric zero-diagonal A,
the specified weak empirical law of A/sqrt(n), and the actual bound
limsup ||A||/sqrt(n)<=5/3. These are asymptotic hypotheses, not exact
finite-order flatness or an exact kernel assumption.

With lambda=sqrt(n/rho), the normalized limiting atoms are 0 and
plus or minus one. The two projector cutoffs plus or minus 1/2 carry
no limiting mass. The operator hypothesis gives a common bounded
spectral support, so the bounded piecewise continuous squared errors
indeed integrate to zero. This proves both Frobenius approximations
in (2.1), as well as the limiting positive and negative ranks.

The exact identities A_ii=0 and (A^2)_ii=n-1 give the claimed diagonal
L2 estimates for T and P. Since 2P_+=P+T, the raw positive-projector
correlation R0 has diagonal approaching one in empirical L2. No
uniform positive lower bound on that raw diagonal is inferred.

For the stated eta, the bad fraction is at most h_n/eta^2, which
tends to zero. On good coordinates the normalization is well defined;
the independent identity block on the remaining coordinates makes R
PSD with diagonal exactly one. Its operator norm is at most
2/[rho(1-eta)], up to the harmless maximum with one.

The row/column removal difference has rank at most twice the bad
count and bounded operator norm, hence Frobenius size O(sqrt(bad
count)). Identity filling costs the same order. The good-coordinate
diagonal normalization costs O(eta sqrt(n)). Thus (2.3) follows.
All original A coordinates remain in the energy and in the Boolean
vector; this correlation repair does not delete source edges.

## 2. Schur powers and the positive original-energy baseline

Schur multiplication by a correlation matrix preserves PSD order
and fixes scalar diagonal matrices. This proves the uniform Schur-
power operator bound, including negative entries of R. In particular
tr R^2 <= ||R|| tr R = ||R|| n.

The sign Hermite expansion has first coefficient sqrt(kappa), total
squared coefficient mass one, and higher odd mass 1-kappa. Its matrix
covariance series converges in operator norm under the uniform Schur
bound. Consequently the PSD tail bound in (3.3) is justified.

For every odd order at least three, the complete off-diagonal entries
and zero diagonal give |tr(A R^{circ q})| <= sum R_ij^2 = O(n),
uniformly in the order. Pairing with the Hermite series therefore
leaves only the first chaos on the n^(3/2) baseline scale.

The Frobenius error R-R0 pairs with ||A||_F=O(n) to give
o(n^(3/2)). The positive spectral atom and rank imply
tr(AR0)=n lambda+o(n^(3/2)). Hence the expectation of Q_A itself,
with its factor 1/2 retained, tends to 5kappa/8 after normalization.
This is neither a bilinear objective nor an absolute-moment substitute.

## 3. Stronger higher-chaos mean and individual variance bounds

The error in (4.1) follows by adding A^2/lambda^2-P,
A/lambda-T, and rho(R0-R). Each is o_F(sqrt(n)). For every higher
odd q, the trace pairing with R^{circ q} has error o(n) uniformly:
the latter Frobenius norm is O(sqrt(n)) uniformly in q.

Moreover tr(R R^{circ q}) is the sum of even powers R_ij^(q+1),
so it is at least the n unit diagonal terms. The subtracted trace
tr(A R^{circ q})/lambda is only O(sqrt(n)). This proves (4.3)
with a uniform error for all higher odd orders.

The tail coefficients are nonnegative with total mass 1-kappa.
Multiplication by lambda^2 and division by n^2, using
lambda^2 rho=n, gives exactly mu_n>=1-kappa-o(1). There is no
missing factor of two and no interchange of a nonuniform error
with the infinite Hermite sum.

The tail covariance cap and the exact row-square sum n-1 give the
uniform individual bound 0<=v_i<=2(1-kappa)/rho+o(1). The proof
does not require homogeneous variances or convergence of their
empirical distribution.

## 4. First-chaos alignment uses a trace of the square

The averaged alignment error is precisely
kappa tr[(A-lambda I)^2 R]/n^2. For R0, spectral commutation makes
this trace (2/rho)||(A-lambda I)P_+||_F^2=o(n^2).

For the R-R0 difference, the actual operator cap gives
||(A-lambda I)^2||_F=O(n^(3/2)); pairing with o_F(sqrt(n))
therefore gives o(n^2). This proves the required variance statement
without using an inverse correlation matrix or mistaking an AR
Frobenius approximation for a covariance estimate.

Orthogonality to the higher chaoses gives c_i=E[G_i F_(i,1)].
The row lengths and covariance cap bound all relevant second moments.
Cauchy--Schwarz applied to the established alignment proves both
averaged errors in (5.3), hence all but o(n) rows have c_i>=a/2.

## 5. Complete joint Gaussianization check

The coefficient assumptions in the extended lemma allow the actual
zero diagonal coefficient. They imply ||diag(d)||<=n^(-1/2)
and ||d||<=1. In the partial-contraction case,
||M||_F^2<=C^2/n and the two remaining Schur powers have norm at
most C, giving C^4/n. In the unequal full-contraction case,
||R^{circ p}d||<=C and the remaining diagonal-weighted Schur power
has norm at most C/n, giving C^3/n.

The new distinguished-coordinate contraction has squared norm at
most (C/n)sum_j R_ij^2 <= C^2/n, because R^2<=C R and R_ii=1.
This covers the interaction of u_i with each higher kernel. The
interaction of g_1 with those kernels is already in the unequal
full-contraction bound. First-order full contractions contribute
constants, just as equal higher-order full contractions do.

Thus every nonconstant term in the gradient product for
sG_i+tF_Q is controlled. The finite product expansion, polynomial
isometry, and triangle inequality give the stated O(1/n) variance
bound for Gamma_Q; independence between chaoses is not presumed.

Gaussian integration by parts gives the printed characteristic-
function differential equation with the correct sign. Its integrating-
factor kernel has modulus at most one along the interval from zero
to the argument, giving the z^2/2 bound without division by variance.
Singular R is covered by its finite-dimensional Gram representation.

The omitted sign tail has variance at most C tau_Q and is orthogonal
to both G_i and all retained chaoses. Thus the characteristic-function
and variance errors of each fixed linear combination tend uniformly
to zero by first increasing n and then Q. Compactness of the two-by-
two covariance matrices and Cramer--Wold give the required Gaussian
pair limit on every covariance-convergent subsequence.

The uniform second-moment bound gives uniform integrability of |F|.
For the disagreement event, c>=c_*>0 forces both limiting marginal
variances to be positive. Its boundary is contained in the coordinate
axes, whose limiting probability is zero even for correlation one.
The Gaussian angular formula and continuity therefore give (6.2),
including its degenerate endpoint. Compactness contradiction supplies
uniformity of both assertions over the full admissible data.

This checks the probability statement as well as the absolute moment.
A variance-only argument, or the older scalar absolute-moment lemma
without the distinguished-coordinate extension, would not suffice.

## 6. Averaging the local formulas

Each actual row has coefficient squared norm (n-1)/n and the required
maximum coefficient. Uniformity permits summing the separate marginal
errors and dividing by n. The average error in b_i-a^2 and continuity
of square root give the absolute-field formula. Subtracting twice the
baseline uses sqrt(kappa)a=5kappa/4 and proves (7.1).

For the sign frequencies, the exceptional c_i<a/2 rows have vanishing
fraction. On the remaining rows the Gaussian pair variance is bounded
away from zero. Restricting first to rows with small b_i and c_i
errors, then using uniform continuity of the angle on its compact
covariance domain, proves (7.2). This handles v_i=0 without an
unjustified Lipschitz assertion at perfect correlation.

No growing-dimensional joint local-field limit is needed, and no
covariance operator bound for sign(AX) is assumed.

## 7. Same-source Boolean update and exact penalty

Conditioning on X makes Y and Delta fixed. The independent Bernoulli
variables select between X_i and Y_i, so every X' is Boolean on the
same original source. The exact quadratic expansion has a first-order
term epsilon(AX)^T Delta and second-order term
epsilon^2 Delta^T A Delta/2. The only potential diagonal exceptions
to the product expectation vanish because A_ii=0.

The first term is epsilon(||AX||_1-X^TAX), including zero fields
with the chosen sign convention. The identity ||Delta||^2=4 times
the disagreement count gives the lower penalty
-2epsilon^2 ||A|| times that count. Division by n^(3/2) proves
(8.2) with the exact coefficient 2epsilon^2 C_n p_n.

Every Q_A(X') is at most Phi(A), so averaging gives a valid original
quadratic norm lower bound. No polarization, doubling of order, or
norm transfer is used. The actual cap 5/3 is retained; weak spectral
convergence cannot replace it by the empirical atom 5/4.

## 8. Contributed scalar step and strict margin

I rechecked this portion as a contributor. The reused pi enclosure
gives 7/11<kappa<16/25. The chord constant exceeds 8/25 because
the positive-square comparison is equivalent to 13kappa>8.
The endpoint enlargement V+o(1) changes the continuous chord constant
by o(1), and the variance means remain bounded. Thus the displayed
asymptotic lower for G_n is valid. The arctangent inequality and
Cauchy--Schwarz give the displayed upper for p_n.

With the admissible fixed epsilon=1/10 and cap 5/3+o(1), the update
lower is (4/125)mu_n-(4/375)sqrt(mu_n)-o(1). The stronger mean
puts liminf sqrt(mu_n) above 3/5. The polynomial 3r^2-r increases
there and has endpoint value 12/25, giving liminf gain at least
16/3125. All vanishing errors are harmless for this fixed probability.

The rational comparisons 16/3125>1/200,
2/5-35/88=1/440, and 1/200-1/440=3/1100 give exactly (1.3).
The unconstrained quadratic optimizer from the optional scalar
discussion is not used and needs no admissibility claim here.

## 9. Transfer composition and final verdict

The separately reviewed 141-line lemma supplies one common original
principal source A_J with q/n tending to one, the precise limiting
law at 0 and plus or minus 5/4, limsup operator cap 5/3, and
Phi(A_J)/q^(3/2)<=Phi(A)/n^(3/2)+o(1). Since q tends to infinity,
the main theorem applies directly and proves the stated paired-source
corollary. Its positive diagonal D is explicitly feasible for both
signs; neither optimality nor a separate trace cap is imposed.

PASS on the complete final 612-line source, with no required correction.
The conclusion is an actual original-source strict lower and excludes
the specified near-scalar internal-law regime at objective 2/5. It does
not retract the earlier FORMAL trace-relaxation calculation: that
calculation never asserted actual complete-signing realizability.

The paired covariance, cross block, active field, and other profiles
are not replaced or optimized here. The all-profile/all-active-cell
implication and the global original MO target remain OPEN.

No mathematical program, solver, checker, numerical integral, scan,
construction, or search was run in this review. Tools were used for
complete reads, hashes and line counts, and this /tmp receipt only.
I made no source changes or canonical repository edits and performed
no gate, publication, commit, or backup operation.
