# Independent new-source review: all-law adaptive nuclear gain

2026-09-06. Reviewer: optimized_profile_docs_gate.

## Frozen complete reads and exact independence

I directly read the ENTIRE 553-line frozen source
`/tmp/original_mo_all_law_adaptive_nuclear_gain.md`, SHA256
`0a7c553e29d4e3ac1572edb0e3fc795bc4d252d090061181365f01764c500a51`.
I also directly read its complete 109-line author receipt,
`/tmp/original_mo_all_law_adaptive_nuclear_gain_author_receipt.md`,
SHA256
`cc47e39b7eead99cbec74d6c684e5046d9217e350b8b7cc7f36a36ba3c097983`.
The source and author receipt hashes and line counts were checked.

All five named prerequisites were directly read completely in this
analytic review sequence, and their hashes were refreshed for this task:

- `original_mo_original_source_near_flat_strict_gain.md`, 612 lines,
  `7726b89e1c39429cde75ff887b981cbd3cf831adb17b04f20193a3c6dbb35298`;
- `original_mo_nuclear_spectral_budget.md`, 147 lines,
  `ee8ad5ff3dbf9aa9e251c4190e98ee1671c9a2140c759ba6f768f8c9c03ef13d`;
- `original_mo_complete_cross_flat_spectral_gain.md`, 411 lines,
  `b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`;
- `original_mo_original_phase_spectral_moment.md`, 262 lines,
  `7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`;
- `original_mo_source_cross_nuclear_trace_boundary.md`, 444 lines,
  `106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`.

The entire 147-line source was newly read in this all-law task; the
other four complete reads occurred in the preceding source-gain review
sequence. In particular the whole final612, not merely its changed
corollary, was read. The reused Gaussian argument is Section 6 of that
source, not its near-flat projector or alignment conclusion.

I supplied no proposal, derivation, parameter selection, correction,
or proof step to this NEW 553-line source. I independently checked the
written proof after its freeze. I previously authored scalar209 support
for the older fixed-probability argument, so my review of the combined
older612 theorem was contribution-disclosed. Scalar209 is NOT a logical
prerequisite of this adaptive argument. I also previously contributed
an AM-GM comparison to the older444 source; that step is unrelated to
the reused pi enclosure. These historical roles are disclosed, without
claiming independence for every earlier item in the provenance chain.
The present complete review is independent of the new all-law derivation.

## 1. Actual normalized phase frame and all-law hypotheses

The input is any actual complete symmetric zero-diagonal A with
||A||/sqrt(n)<=C, C fixed. No weak spectral limit, flatness, symmetry
of the eigenvalue law, diagonal homogeneity, or feasible SDP diagonal
is assumed. With M=A/sqrt(n) and q=(n-1)/n, exact complete row
squares give (M^2)_ii=q.

Spectral calculus gives M^2=|M|^2<=C|M|. Taking diagonals supplies
h_i>=q/C; Cauchy--Schwarz in the coordinate spectral measure gives
h_i<=sqrt(q). These prove every lower and upper bound in (3.1),
including positivity of all phase normalization diagonals.

Both |M|+sM are PSD and have exactly the same diagonal h. Therefore
their coordinate-normalized R_s are genuine unit-diagonal correlations.
The operator cap is at most (C/q)(2C)=K_n=2C^2/q<=4C^2.
No replacement of H by its mean ell I is made or needed. Singular
correlations and independently sampled positive/negative phases are
permitted throughout.

## 2. The exact baseline belongs to these same phases

The phase-correlation difference is 2M_ij/sqrt(h_i h_j). On each
complete edge, the oriented arcsine difference multiplied by A_ij is
at least 2/(sqrt(n)sqrt(h_i h_j)). Normalizing the energy difference
by n^(3/2) and averaging the two oriented phases gives exactly the
first inequality in (4.1), including its factor kappa/n^2.

The unordered-pair AM-GM and Cauchy--Schwarz denominator is
(n-1)sum_i h_i. With k=n(n-1)/2 this yields
sum_(i<j)(h_i h_j)^(-1/2)>=n(n-1)/(2ell), and hence the exact
average baseline kappa q/(2ell).

This is a lower on the average expectations of the very phases later
updated. The proof does not try to add a gain to a previously maximized
nuclear-norm lower bound belonging to unrelated states.

## 3. Uniform higher-chaos lower with unequal diagonals

Positive Schur multiplication preserves PSD order and fixes scalar
diagonals, giving the common operator cap K_n for every positive
entrywise power of R. In particular tr R^2<=K_n n.

The actual off-diagonal entries of M have magnitude n^(-1/2), so
|tr(M R^{circ k})|<=K_n sqrt(n), uniformly for odd k>=3.
The exact identity |M|=H^(1/2)R H^(1/2)-sM gives the trace sum
sum_(i,j)sqrt(h_i h_j)R_ij^(k+1). Since k+1 is even, this is at
least its diagonal sum n ell. This proves (5.3) without commutation
of H with M or an approximation to either matrix.

For t>=0, the square (|M|-tI/2)^2 is PSD, so
M^2>=t|M|-t^2 I/4. Pairing with the PSD Schur power and choosing
t=2(ell-K_n/sqrt(n))_+ gives precisely the finite lower (5.4).
The same error bound works for every higher odd chaos order.

The Hermite tail has nonnegative coefficient mass 1-kappa and an
operator-convergent covariance series. It can therefore be summed
against the uniform bound. This proves the finite mean-variance lower
in (5.7), not merely a nonuniform per-chaos assertion.

First-chaos orthogonality and Cauchy--Schwarz give c_i^2<=b_i and
w_i=b_i+v_i-c_i^2>=v_i. The actual row-square identity gives
sigma_i^2<=K_n q=2C^2. Taking the trace with M^2<=C^2 I gives
the stronger averaged first-chaos bound avg c_i^2<=kappa C^2.
Neither positivity of the individual c_i nor any common row covariance
is assumed. Since ell<=1 and K_n<=4C^2, the final O_C(n^(-1/2))
interpretation of (5.7) is uniform.

## 4. Joint marginal limit and the correct uniform-integrability use

The actual row coefficients have one zero diagonal entry, maximum
magnitude n^(-1/2), and squared norm q<=1. They meet the generalized
Gaussian-sign lemma's coefficient hypotheses in both phases.

The source correctly restates all required contractions: B^4/n for
partial contractions, B^3/n for full unequal-order contractions, and
B^2/n for a distinguished input coordinate contracted into a higher
kernel. Equal-order full contractions are constants. These bound all
nonconstant terms of the finite-chaos gradient product for any fixed
linear combination of the input Gaussian coordinate and the local field.

The Gaussian integration-by-parts equation and orthogonal sign tail
give the joint two-dimensional limit after taking n first and the
truncation order second. Compact covariance subsequences make the
approximation uniform. No growing-dimensional local-field limit or
independence between the original spins is required.

The new payoff g_C((-sign(u)v)_+) is continuous except possibly at
u=0; the standard Gaussian first marginal assigns that line probability
zero. The payoff is nonnegative and bounded by 2|v|. The uniform
second-moment bound therefore makes this linearly growing payoff
uniformly integrable. It does NOT establish uniform integrability of
F^2, and the proof correctly makes no such inference.

Both the actual pair and the comparison Gaussian pair have compactly
bounded covariance parameters. Applying the same continuity and
integrability argument on a violating subsequence proves (6.2)
uniformly, including zero or degenerate second marginals and either
sign of the covariance c. The Gaussian comparison keeps the ACTUAL
c_i and w_i rather than imposing an aligned representation.

## 5. Exact adaptive Boolean update and normalization

Conditioned on X, all r_i, epsilon_i, Y, and Delta are fixed. The
Bernoulli masks are independent only under that conditioning. Since
the diagonal of sM is zero, the quadratic expansion has exactly the
linear term (sMX)^T z and the quadratic term z^T(sM)z/2,
z_i=epsilon_i Delta_i. There is no missing diagonal variance term.

A nonzero mismatch has F_i Delta_i=2r_i and Delta_i^2=4. Matches
have r_i=epsilon_i=0, and a zero field also has epsilon_i=0 even
if its chosen sign disagrees with X_i. Consequently the spectral lower
is the sum of 2epsilon_i r_i-2C epsilon_i^2 over all coordinates.

The choice min(r_i/(2C),1) is the maximizer over the ADMISSIBLE
interval [0,1]. Its value is r_i^2/(2C) below the cutoff and
2r_i-2C above it, exactly g_C(r_i). Every updated vector is Boolean
on the same source. Multiplication by sqrt(n) converts sM energy
to sA energy, so its normalized improvement is avg g_C(r_i),
as stated in (7.2). No bilinear or polarization factor is introduced.

## 6. Uniform Gaussian radial clipping estimate

An arbitrary comparison pair can be written using two independent
standard Gaussians even if its covariance has rank at most one. Their
radius T is independent of angle with density t exp(-t^2/2).
The mismatch event depends only on angle, and its negative field has
the form bT with 0<=b<=sigma at each such angle.

Expanding (t-a)^2 t and integrating gives

    E(T-a)_+^2 = 2exp(-a^2/2)
                  -2a integral_a^infinity exp(-u^2/2)du.

Division by E T^2=2 proves (8.1). At a=2C/b, its upper is the
unclipped radial second moment times exp(-2C^2/b^2). The actual
variance bound sigma^2<=2C^2 makes this factor at most exp(-1).
Zero b and nonmismatch angles contribute nothing. Averaging therefore
proves the clipping loss bound and (8.3), with the correct 1/(2C)
coefficient. Neither a mismatch-probability denominator nor an
uncontrolled tail removal occurs.

## 7. Negative-field moment and joint convex perspective

The representation r distributed as (-c|Z|+sqrt(w)Z')_+ makes
replacement of negative c by |c| a pointwise decrease. For c>=0,
the two angular mismatch sectors and radial second moment two give
the printed Psi(c^2,w). Its limits at c=0 and w=0 are w/2
and zero respectively. Thus (9.2) holds for all covariance signs.

I differentiated the displayed f directly. Both
f'=[arctan(x^(-1/2))-x^(-1/2)]/pi and
f''=1/[2pi x^(3/2)(1+x)] are correct. The perspective Hessian is
f''(a/w)/w times (1,-a/w)(1,-a/w)^T and is PSD. The continuous
boundary extension is valid; in particular 0<=Psi(a,w)<=w/2
controls the origin. Hence joint convexity extends to the closed quadrant.

Psi decreases in a and has w derivative
arctan(sqrt(w/a))/pi in [0,1/2]. Jensen can therefore be followed
by the upper bound on avg c_i^2 and lower bound on avg w_i in
exactly the directions used in (9.5). Convexity in w alone would not
supply this conclusion. No local parameter distribution need converge.

The common resulting gain lower holds for each phase separately.
Averaging their two valid lower bounds on alpha adds that gain to
the SAME-phase baseline from (4.1), giving (9.6).

The bound ell>=q/C keeps 1/ell uniformly bounded for fixed C.
Replacing q by one costs O_C(1/n). Replacing the second Psi
argument costs O_C(n^(-1/2)), since its w derivative is at most
1/2. Together with the rowwise uniform Gaussianization error this
proves the claimed uniform e_C(n)->0 in (1.2).

## 8. Integral envelope and the entire stated nuclear region

Integrating the w derivative from zero gives the displayed integral
for Psi. The inequality arctan x>=x/(1+x^2) follows by derivative
comparison. Bounding t+c^2 above by w+c^2 inside the integral gives
2c w^(3/2)/[3pi(c^2+w)], with its positive term preserved.

The strict series bound exp(1)>8/3 gives 1-exp(-1)>5/8.
Substituting c=C sqrt(kappa), w=(1-kappa)ell^2, and pi=2/kappa
therefore gives exactly the coefficient 5/48 in J_C. Both gain
functions are positive for fixed C and ell>0.

For C=5/3 and 3/4<=ell<=4/5, the reused kappa interval implies
kappa(1-kappa)>(12/25)^2. The denominator increases in ell,
and at ell=4/5 also increases in kappa. Its stated upper
16/9+144/625=11296/5625<81/40 is correct. The resulting envelope
constant simplifies exactly to 3/1250. The baseline is at least
5kappa/8>35/88.

For ell<=3/4 the baseline alone exceeds 14/33, and
14/33-35/88=7/264>3/1250. Thus the entire smaller-ell region is
covered, not just the endpoint ell=4/5. The final identity
35/88+3/1250=2/5+7/55000 is exact.

Under the two limsup assumptions, any relevant subsequence has a
further ell-convergent subsequence with limit in [3/5,4/5]. The
uniform theorem can first be applied at each fixed C=5/3+epsilon;
only afterward is epsilon sent to zero using continuity. This avoids
an unjustified moving-C error bound and handles alternating nuclear
subregions or approach to 4/5 from above. It proves (1.5) without
an empirical spectral limit or convergence of the original ell sequence.

## Verdict and precise remaining scope

PASS for the ENTIRE frozen 553-line source, with no required correction.
This is an independent check of every new all-law link, with earlier
contribution history disclosed above. It proves the uniform nuclear
gain and excludes the whole stated actual bounded-operator/nuclear-
moment region at normalized original objective tending to 2/5.

The earlier narrower near-flat theorem remains valid with its stronger
gap; it is not rewritten or retracted. This new theorem removes its
spectral-law and diagonal-homogeneity premises, not the ACTUAL operator
cap. No arbitrary optimizer is proved to satisfy that cap or to lie
inside the nuclear region. The complementary source region, paired
all-cell implication, global original MO convergence, and the possible
limit value remain OPEN.

No mathematical program, checker, solver, numerical integration,
optimization, signing construction, or search was run. Tools were used
only for complete reads, hashes and line counts, and this /tmp receipt.
I made no canonical edit and performed no publication or backup operation.
