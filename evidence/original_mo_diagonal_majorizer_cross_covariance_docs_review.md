# Independent full review: diagonal-majorizer cross covariance

2026-09-05. Reviewer: optimized_profile_docs_gate.
Verdict: PASS after the series-domain clarification described below.

## Exact source, complete read, and independence

The reviewer read every line of the initial 381-line source
`/tmp/original_mo_diagonal_majorizer_cross_covariance.md`, SHA-256
`e6af3fc933ed6602fd24a98cf2b1801395f352c54a777f7cb81fb055a20c3961`.
The reviewer then completely checked the revised Section 5, including
both changed locations. The final 384-line source has SHA-256
`0b3921d43d88424457ad2ed777ee158e8ac34c6751c995f0c3b86aee870e95ff`.
This is a full initial read plus a complete final-delta review, not a
claim that the original 381-line hash remains the final version.

The reviewer did not contribute to this construction or its derivation
and made no source edits. The author noticed, and the reviewer agreed
with, the necessary series-domain clarification: epsilon<=1/2 now
precedes the definitions of the separately retained infinite series.
Without that restriction the individual cancellation series need not
converge, even though the original centered covariance is defined.
Section 6 already handles the bounded set of remaining orders directly.
No covariance, constant, theorem scope, or asymptotic rate changed.

The reviewer also completely read these current prerequisites:

- `NOTE_2026-09-05_SHIFTED_SIGN_GAUSSIAN_UNIVERSALITY.md`, 317 lines,
  SHA-256
  `a3ed6d9c3ee73b863c91d069e75baf9973911318a8efe9156ca61e30f55d7e25`.
- `NOTE_2026-09-05_DIRECT_CROSS_COVARIANCE_NORMALIZATION.md`, 265 lines,
  SHA-256
  `e4919c8e16461c35efdf2963eaf9fdc1b45c07ccfba33ae1549a07e904f7ac8a`.

Their actual hypotheses and normalization factors were checked against
the new construction, rather than inferred from similar terminology.

## Attained same-diagonal majorizer and local scales

The bipartite vector SDP dual has the stated factor one-half in its
objective and the stated unscaled off-diagonal copies of K. Symmetry
of K permits swapping the dual diagonal blocks and averaging. The
orthogonal sum/difference decomposition then gives D-K and D+K, with
trace D equal to the dual optimum, not twice or half that optimum.
Strict feasibility and compact nonnegative diagonal trace sublevels
give attainment. A zero diagonal entry would force a nonzero signing
row to vanish in a PSD constraint, so every d_i is strictly positive.

The finite tensor-rounding constant in (1.1) and the cube-polarization
factor four yield exactly S<=4g Phi(K). Conversely the two Loewner
majorizations give |z-transpose K z|<=S on every Boolean z, hence
Phi(K)<=S/2. These bounds do not invoke a numerical SDP calculation.

The normalized symmetric T is a contraction. Its row-square bound
gives d_i>=sum_{j!=i}1/d_j. Cauchy--Schwarz proves both subsequent
lower bounds in (2.1), and summing the first inequality proves (2.2).
Completeness and N>=4 imply d_i d_j>1 for distinct coordinates: the
remaining strictly positive row-square terms exclude equality.
The lower bounds on the two diagonal entries give precisely (2.3).
There is no unstated upper bound d_i=O(sqrt(N)) at each coordinate.

## Symmetric compression, padding, and exact PSD

The symmetric cross-edge basis is orthonormal with its factor sqrt(2).
Compression of T tensor T gives exactly (3.1), including the negative
sign on the second copy of A and the orientation B_il B_kj.
The exchange operator B X-transpose B is self-adjoint on these matrix
coordinates even when B is nonsymmetric.

The tensor contraction makes I-T tensor T positive semidefinite and
bounded by 2I. Its compression therefore gives (3.2), with diagonal
1-q_e-squared. Adding the independent diagonal covariance Q-squared
restores EXACT unit variance and proves (3.3)--(3.5). The Gaussian
law is defined from a PSD matrix at the outset; no indefinite formal
field or positive-part repair is introduced.

For epsilon<=1/2 the diagonal conjugation bounds stated before (3.6)
give an operator difference at most
[2(sqrt(2)+1)+1]epsilon<6epsilon. The literal normalized alternative
is therefore valid, while the proof correctly uses the separable
padded construction for its degreewise estimates.

## Complete weighted Hermite remainder

The mean, first Hermite mass, remaining variance, and orientation
parity in Sections 4--5 are correct for every deterministic real h.
Only the even Hermite coefficients acquire the entrywise sign of B.
The covariance diagonal is already accounted for in v_h I.

For the four-cycle norm, squaring the full ordered-pair matrix gives
the nonnegative entrywise squares claimed in the source. The feature
matrix has V-transpose V=(K-squared) Schur (K-squared). Applying the
positive Schur map to K-squared<=L-squared I gives its operator bound
L-squared times (N-1). Each feature row has squared norm at most N-1,
so the row-sum bound for the squared matrix is correct. Compression
to symmetric off-diagonal coordinates gives TWICE Q_N. This proves
the factor one-half in (5.1) without a conference-scale hypothesis.

The entrywise identity (5.2) is correct on disjoint, adjacent, and
identical cross edges. On disjoint edges the four-cycle sign is the
negative product of the two numerator signs, so magnitude two occurs
exactly for Q_N=1. On adjacent edges the magnitude is one, and the
last identity term cancels the diagonal exactly.

Under the now-explicit epsilon<=1/2 hypothesis, the separate series
in (5.3)--(5.4) converge, and the powers of the local q_e weights are
exact. The retained covariance is a sum of positive rank-one profiles,
not generally a single scalar multiple of b b-transpose. Its diagonal
cancellation with E_even is retained rather than counted twice.

Bounding the three terms of (5.4) separately gives precisely
epsilon-squared times [L(N-1)+4(N-2)+2]. The odd remainder starts at
degree three, and its row count gives at most 8n-squared epsilon-cubed.
Thus (5.6) has the claimed O((L+1)/N) rate and does not require a sign
for its error matrix. Both comparison covariance matrices are PSD.

Each retained coordinate variance is at most 2q_e-to-the-fourth.
Summing Gaussian expected absolute entries bounds the ACTUAL Boolean
norm. Separability of q_e-squared, the two-part product inequality,
and (2.2) give exactly the constant in (5.7), hence O(N). This is a
valid norm-cost estimate on all retained profiles simultaneously;
low rank alone is not used as a norm bound.

## Gaussian comparison, uniform universality, and padding removal

Adding independent delta-I noise gives Gaussian convex order in both
directions for the finite maximum with arbitrary deterministic offsets.
The states have squared norm at most n-squared and cardinality at most
2^(2n+1), including the absolute-value sign. The resulting comparison
cost and (5.7) give exactly O(n sqrt(L+1)+n) in (6.1).

The shifted-sign theorem applies on all n-squared actual cross entries:
the latent covariance has exact diagonal one and operator bound three,
the threshold and orientation are deterministic, the coefficient
states lie in the real unit cube, and fixed internal energies can be
absorbed into the deterministic prior. The complete growing-temperature
estimate, including its h-uniform constants and singular-covariance
allowance, yields the stated n^(16/11) error. No threshold is selected
after seeing the disorder, and no fixed-temperature limit is exchanged.

The real/complex interpolation factor two and cube-polarization factor
four give L-squared<=8Phi(K)<=4S. Therefore n sqrt(L+1) is
O(n^(11/8)), strictly lower order than n^(16/11). The direct bounded
norm treatment of the finitely many smaller orders is uniform in h
and in arbitrary deterministic offsets, and does not use the separated
series outside their justified domain.

The removed padding is independent Gaussian diagonal covariance,
not discarded variance in the actual shifted-sign law. Every Boolean
pairing has variance sum q_e-squared. The finite-maximum estimate gives
the exact bound (6.2), so its expected norm cost is O(n). Lipschitz
comparison proves the second Gaussian model's claim with the SAME
actual sign law. Equation (3.6) also gives O(n) at the Gaussian-maximum
scale, but is not needed for the two asserted models.

## Exact weighted process and conditional-original scope

Contraction of the two covariance tensor terms gives exactly (7.1):
the internal factors are A_L=UAU and A_R=VAV, while the exchanged
cross factors use W_D=UBV in the displayed order. Independent padding
gives exactly (7.2). The differently weighted copies of A cannot be
silently replaced by the unweighted source energies.

The deterministic drift remains s_h B and the actual internal energy
remains Q_A(x)-Q_A(y). No source, cross signing, or conditional objective
was modified by introducing D. The trace cap is not claimed to provide
the still-needed weighted/unweighted energy compatibility.

The independent cross-sign estimate supplies the conditional norm cap
with the correct constants. Every rounded cross signing is admissible
for the SAME original conditional minimization, so optimality gives
the floor direction in (8.1). Gaussian universality transfers that floor
with the already proved error; it does not reverse it into an upper
bound or turn conditional norm optimality into pressure optimality.

No evaluated Gaussian upper or original MO convergence follows from
this source alone. That scope is stated consistently throughout.

No mathematical computation, certificate execution, signing census,
simulation, optimization, solver, or new test was run by this reviewer.
