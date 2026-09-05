# Independent full review: canonical gap and two original quadratic phases

2026-09-05. Reviewer: optimized_profile_docs_gate.
Verdict: PASS; no corrections requested.

## Exact source, dependencies, and independence

The reviewer directly read every line of the 274-line source
`/tmp/original_mo_full_sdp_gap_original_phase_bound.md`, SHA-256
`1d36878bdd157be36b1e935f0e92a0e977cbbabb1bbf23784a645860ac1142c0`.
The reviewer did not contribute to this theorem or its derivation and
made no source edits. This artifact is separate from the two-note gap
compatibility publication that was already being prepared.

The reviewer also read completely the 262-line original-phase spectral
moment prerequisite, SHA-256
`7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`.
The canonical-gap prerequisite, 303 lines and SHA-256
`3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`,
had already received the reviewer's full independent mathematical read.
The current direct-cross normalization prerequisite, SHA-256
`e4919c8e16461c35efdf2963eaf9fdc1b45c07ccfba33ae1549a07e904f7ac8a`,
had also been completely read, including the real/complex interpolation
constant, cube polarization, and the source-scale norm consequence.

## Actual positive phases and original-norm normalization

The matrices K-squared plus or minus K|K| are twice the squares of the
positive and negative spectral parts. Their PSD diagonals imply
|h_i|<=q. The two padding diagonals in (2.1) are nonnegative, and
both phases have EXACT diagonal one after the same coordinatewise
normalization v_i=q+|h_i|. Singular Gaussian correlations are allowed.
Neither phase is replaced by a formal indefinite covariance.

The Gaussian sign identity has the correct sum over unordered edges
and original quadratic half-value. The arcsine Taylor remainder bound
is valid on the entire correlation interval, including both endpoints.
Multiplying its constant by kappa gives rho, so each expected-energy
remainder is at most rho/2 times the off-diagonal Frobenius norm squared.

The lower bound v_i>=q gives the two fourth-moment estimates exactly.
The padding is diagonal and hence absent from both those off-diagonal
estimates and the trace pairing with the zero-diagonal K.
Because the SAME P is used in both phases, their K-squared trace terms
cancel in the difference. The remaining linear term is kappa J and
the sum of the two error bounds is at most 2rho S_4/q-squared.

Each expected ORIGINAL quadratic energy lies between minus Phi(K) and
Phi(K). Their difference is therefore at most 2Phi(K), proving
(2.2) with coefficient kappa/2 and the stated fourth-moment subtraction.
This normalization uses two original source states directly. It is not
obtained by transferring a rectangular Boolean bound through factor four.

## Gap-controlled diagonal mask and weighted loss

The residual R_1=DK-K|K| is one nonnegative summand of the complete
weighted residual theorem, so (3.1) has right side 4qg. Its diagonal
is minus h_i because K has zero diagonal. Weighted Cauchy--Schwarz
then gives both bounds in (3.2), with the exact factor two.

The normalization weights satisfy the stated f_i, a_i, and m_ij
ranges; in particular m_ij<=1/2. Integrating the derivative of
1-(1+x)^(-1/2) gives its upper x/2 bound, hence the displayed A_0
estimate. B_0 sums only off-diagonal entries, exactly q per row, so
the cap B_0<=S/2 has no spurious diagonal contribution.

The dispersion identity and weighted Cauchy--Schwarz give
sum|d_i-dbar|<=S sqrt(delta). Since total signed deviation is zero,
its positive part is exactly half the absolute sum. Combining this
with 0<=a_i<=c_* gives the weighted a_i estimate in the source.

Summing m_ij<=a_i+a_j over OFF-diagonal entries gives (3.5) after
discarding only a nonpositive correction. With
sqrt(Sg/q)=eta N sqrt(gamma) and sqrt(delta)<=2eta sqrt(gamma),
division by S gives exactly
eta(1+N/q+c_*)sqrt(gamma). Retaining the separate half-cap proves
(3.6) with precisely the b_* in (1.1). No maximum diagonal or maximum
phase variance has replaced the actual weighted distribution costs.

## Masked residual, finite inequality, and asymptotic error

The literal sign-square identity is used only for i!=j. Multiplying
the unnormalized canonical objective by f_i f_j=1-m_ij gives the exact
identity (4.1), including the PLUS sign on its masked residual term.
The source explicitly omits the diagonal throughout this identity.

Weighted Cauchy--Schwarz, the residual bound 4qg, and
m_ij-squared<=m_ij/2 give exactly sqrt(2gB_0) in (4.2).
The loss b+sqrt(2gamma b) is increasing for b>=0, which justifies
substitution of the upper b_*. Finally S_4<=L S_3 and
S_3=qS(1-gamma) give precisely (1.2) from the original-phase bound.
The finite inequality remains valid when its right side is negative.

Under a fixed ORIGINAL norm cap, the reviewed inequalities
L-squared<=8Phi(K) and S<=4G Phi(K) bound eta and give
L/q=O_C(N^(-1/4)). The mask loss is O_C(sqrt(gamma)); the additional
gamma and gamma^(3/4) losses are no larger on [0,1). Multiplying by
S yields the exact two error scales N^(3/2)sqrt(gamma) and N^(5/4)
in (1.3). The lower bound S>=N sqrt(q) then permits division by S
to obtain its leading small-gap statement. The rectangular estimate
is used here only to control S under a norm cap, not to create the
original-phase coefficient kappa/2.

## Active cross normalization and the exact large-gap limitation

For a positive actual cross energy c, inversion of the small-gap
source bound gives both inequalities in (5.2) with the correct
directions. The ratio c/Phi(K) is retained on a general shell.
At original internal values zero, c=Phi(K) is a SEPARATE saturation
premise, not an implication of small gap. Only with that premise does
the uniform compatibility estimate give 2c_D/N>=kappa-o_C(1).
The actual weighted cross field and its unevaluated width remain.

For gamma>=1/4, eta>=1 and N/q>1 force the half-cap b_*=1/2.
The main bracket is then 1/2-gamma-sqrt(gamma), at most minus one-quarter.
The nonpositive fourth-moment subtraction cannot improve it. Thus
the displayed gap-only lower bound is genuinely vacuous on that range.
This is a limitation of THIS evaluated formula, not an assertion that
such a gap occurs at every optimizer or that all optimized two-phase
constructions are impossible there.

The theorem establishes a conditional actual small-gap normalization.
It does not infer small gap from original or conditional minimality,
settle the complementary positive-gap branch, evaluate the remaining
weighted Gaussian upper, or prove original MO convergence.

No mathematical computation, certificate execution, phase simulation,
signing construction, census, numerical evaluation, solver, optimization,
or new test was run by this reviewer.
