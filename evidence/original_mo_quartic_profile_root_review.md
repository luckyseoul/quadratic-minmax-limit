# Complete root review: same-order pressure and quartic profile regularization

2026-09-05. Both all-orders proof notes were read completely, including
every displayed inequality, endpoint and scope qualification. No new
mathematical computation was used. The original limit remains OPEN.

## Reviewed artifacts and exact provenance

- `NOTE_2026-09-05_SPECTRAL_REGULARIZATION_PRESSURE_PROFILES.md`, 152 lines,
  SHA-256 `2f9f63f603fcae42a952fbae53a2301eaa6b95bbe7bac2e35bcab8997d28d7d7`.
  The canonical file is byte-identical to the frozen author source
  `/tmp/original_mo_spectral_regularization_pressure_profiles.md`.
- `NOTE_2026-09-05_QUARTIC_PENALIZED_PROFILE_IDENTITY.md`, 444 lines,
  SHA-256 `ad393709abb35ed760986b102e1b86ab4d23c80261efec04f35d03104c821013`.
  The canonical file is byte-identical to the final frozen author source
  `/tmp/original_mo_quartic_penalized_profile_identity.md`.

Separate complete independent reads are recorded in
`original_mo_pressure_profile_exact_review.md` and
`original_mo_quartic_penalized_profile_exact_review.md`. The root's
final complete read includes all substantive extensions and the explicit
order-one, tensor-lift and nonnegative-diagonal conventions.

## Same-filler pressure proof: PASS

The old incident-edge operator and Boolean failures have total probability
less than one half. The two-phase, q-point Markov union costs one quarter;
therefore ONE filler meets all requirements. Its partition-ratio expectation
is exactly the incident-edge cosh product. Principal restriction lowers
each actual phase pressure by conditional Jensen, without bounding the
removed old energy. The simultaneous norm bound supplies both Lipschitz
constants, and floor-grid interpolation has error `C+C'`. This verifies
every temperature in the compact interval, not merely the grid.
The arithmetic and geometric phase means inherit the same bound.

## Quartic profile proof: PASS

1. The physical profile has equal squared row norms
   `d=1-(2-t)/N`. A valid integral comparison signing bounds the actual
   penalized minimum. All discrete edge flips remain admissible.
2. The full finite quartic flip is
   `-16m(M^3)_ij+32dm^2-16m^4`; no edge-dependent Taylor remainder is
   replaced by a common Hessian. It is the PENALIZED gap that is nonnegative.
3. Conditioning on one spin gives the nonnegative actual local-field
   energy. Combining its exact row expansion with Hoeffding retains the
   negative fourth-moment term and proves the displayed bound for EVERY
   diagonal entry of M^4, as well as bounded gap row sums.
4. The entrywise cubic-matrix sum is a vector bilinear pairing. The
   self-contained tensor rounding, Boolean polarization and row-fourth
   bound give `sum|M^3_ij|=O(N^(3/2))`. This yields the stated actual
   signed-Gibbs bound; no unsigned phase-covariance assertion is used.
5. The balanced finite-flip remainder is exactly
   `6-2t-(12-6t)/N`, with integral `5-9/N`. The pressure expansion has
   the displayed fourth-derivative remainder. Its integrated error is
   `O_c(sqrt(N))`, uniformly for `0<lambda<=1`.
6. Finite-branch absolute continuity and the uniform interior error
   justify integration through the zero-cross endpoint. The endpoint
   penalty is exactly block-additive with the smaller-order normalization.
   Phase pairing gives twice the penalized HALF-PRODUCT optimum, which is
   at most twice the symmetric optimum; the two optima are not identified.
7. The same recompletion has quartic trace at most `(K+8)^2(N-1)`.
   Its actual pressure cost is `O_c(N/K)+O(1)`. Taking
   `K=lambda^(-1/3)` proves the normalized approximation claimed.
8. The diagonal SDP majorizer can be traced against the fourth powers
   of the positive and negative spectral parts. It gives
   `tr|M|^5<=tr(D M^4)<=D_0 tr D=O_(c,lambda)(N)`. Commutation of D
   with M is neither true in general nor required.

## Unresolved implication retained

The signed internal/cross gap integral is not controlled by any theorem
in this milestone. Nonnegative gaps, bounded row sums and vanishing largest
entry do not prove its favorable sign. A permutation average cannot replace
the derivative of the optimized envelope. Even a future dyadic comparison
still needs a valid all-orders transfer argument. No original-norm or
half-product optimizer identification, Gaussian universality, convergence
or candidate limiting value has been inferred.

The documentation scanner is only an overclaim/status check. Its receipt
does not certify any mathematical proof. Source and test implementations
are unchanged in this milestone.
