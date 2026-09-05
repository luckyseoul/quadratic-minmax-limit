# Independent exact review: quartic-penalized balanced profile

Reviewed the complete final 444-line draft
`/tmp/original_mo_quartic_penalized_profile_identity.md`, SHA256
`ad393709abb35ed760986b102e1b86ab4d23c80261efec04f35d03104c821013`.

Result: PASS. All mathematical constants and cancellations checked.

- The physical profile has constant squared row norm
  `d=1-(2-t)/N`; the net-constructed comparison signing gives the
  stated uniform fourth moment budget for actual penalized minimizers.
- The full quartic finite-flip expansion has leading term
  `-16 m (M^3)_ij` and exact remainder `32 d m^2-16 m^4`.
  Penalized optimality concerns `g_e=phi_e+lambda Delta_e V`, not
  the unpenalized gap. The exponential inequality is used with the
  correct direction, even when the quartic flip is negative.
- Block reversal establishes nonnegative internal and cross radial
  derivatives. The resulting signed l1 bound, including both C_0 and
  C_1, follows from the explicit fourth-moment bound without assuming
  any unsigned phase-covariance estimate.
- The new actual row expansion gives (7a) with all negative terms
  retained. It implies the uniform diagonal bound
  `(M^4)_ii<=2+c^2/(8lambda)` and the sharper gap-row bound
  `sum_j g_ij<=2c^2+16lambda`.
- The vector representation of the full entrywise l1 norm of M^3
  has row factors of norms sqrt(N) and sqrt(D_0). The explicitly
  recalled Krivine tensor lift and Boolean polarization give exactly
  `4 K_0 C_Phi sqrt(D_0) N^(3/2)`. Passing to unordered edges gives
  the stated C_* and the upgraded signed sum (7e).
- The improved integrated error is O_c(sqrt(N)), uniformly for
  `0<lambda<=1`: C_Phi is bounded on that interval and
  `lambda sqrt(2+c^2/(8lambda))` is bounded. This does not require
  uniform boundedness of D_0 by itself.
- The leading weighted quartic flip sum is exactly `-V'(t)`.
  Its remainder is exactly `6-2t-(12-6t)/N`, with integral `5-9/N`.
  At N=2 the remainder reduces to t, which provides an additional
  direct consistency check of the normalization.
- The fourth-derivative Taylor remainder bound is valid uniformly in
  r, and its two error terms give exactly the stated inequality (12).
  Active-branch absolute continuity and the resulting identity also
  justify integrability at t=0.
- At t=0 the smaller-order spectral penalties have exactly their own
  normalization. Choosing A and -A proves the half-product endpoint
  equality, not equality with the symmetric-pressure optimum.
- The symmetric-pressure minimizer has the asserted (loose but valid)
  norm cap C_c. Fixed-temperature phase Markov bounds select the SAME
  filler as the operator and Boolean requirements, with positive
  joint probability. The fourth-moment bound is (K+8)^2(N-1), not
  a fourth power of K. Thus K=lambda^(-1/3) yields the claimed
  normalized O_c(lambda^(1/3)) approximation.
- For the fifth moment, the nonnegative diagonal majorizer dominates
  both M and -M. Taking PSD traces against the fourth powers of the
  positive and negative spectral parts gives
  `tr |M|^5<=tr(D M^4)<=D_0 tr D`. The last inequality uses that D
  is diagonal and nonnegative. This proves (16) and the N^(1/5)
  operator bound, with no commutation assumption on D and M.

The final note explicitly defines the order-one zero matrix, so its
N=2 endpoint also has a completely specified domain. All 444 lines
were freshly reread after the substantial rowwise and fifth-moment
extensions and the final domain/constant clarifications.

The mixed reset-gap integral is still unproved. This review supplies
no sign for it, no selected cross-size construction, and no proof of
the original limit. No computation or new experiment was used.
