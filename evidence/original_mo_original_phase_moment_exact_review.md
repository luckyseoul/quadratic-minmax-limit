# Independent complete review: original-phase spectral moments

2026-09-05. Reviewer: optimized_profile_exact. Verdict: PASS.

The complete 262-line artifact
`/tmp/original_mo_original_phase_spectral_moment.md` was independently
read and checked, covering every formula (1)--(14).

Approved SHA-256:
`7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`.

## Checked implications

1. Both Gaussian phase matrices are PSD, their separate diagonal
   paddings are nonnegative, and their individual normalizations give
   diagonal-one correlation matrices. Zero diagonal of A removes the
   padding from its trace pairing. The two linear expectations are
   respectively `kappa S3+/v+` and `-kappa S3-/v-`.

2. The positive arcsine Taylor coefficients give the uniform endpoint-
   valid square remainder bound. The sum over unordered pairs gives
   exactly `2 delta S4+/v+^2` and its negative-phase counterpart.
   Multiplication by each phase's own variance gives (2). The bounds
   `2 S2+<=n v+` and `S4+<=a^2 S2+`, with their negative analogues,
   give the final error `delta n(a^2+b^2)` in (3).

3. The minimal admissible phase normalizations are strictly positive,
   and their sum is exactly `2q+osc diag(A|A|)`. The bounded-operator
   interpretation of the error in (9) is correctly restricted and
   does not silently cover arbitrary operator growth.

4. The local spectral measure has mean zero and second moment q.
   Its positive and negative first moments agree and are positive.
   The inequalities `P2<=a t`, `p+>=t/a`, and `p->=t^2/N2` give
   `a^2 N2>=q P2`, hence `N2>=q^2/(a^2+q)` and precisely the printed
   upper bound on h_i. Applying the same argument to -X proves the
   lower bound. The stated two-point equality distribution has the
   required mean and second moment when its support is permitted.

5. The universally admissible variances in (11), their sum in (12),
   the local inequality `ab>=q`, the range of the parenthesized
   factor, and the constants in (13) all check. Neither spectral
   symmetry nor separate assumptions `a,b>=sqrt(q)` are required.

6. The covariances `|A|+A` and `|A|-A` have a common strictly positive
   diagonal. Their normalized correlation difference, multiplied by
   each signing coefficient, is `2/sqrt(d_i d_j)`. The arcsine slope
   inequality therefore yields the first lower bound in (14).
   Pairwise arithmetic-geometric mean followed by Cauchy--Schwarz
   gives exactly `n^2(n-1)/(2 tr|A|)` for the pair sum. The local
   Cauchy--Schwarz bound `d_i<=sqrt(n-1)` recovers the stated CORE
   uniform lower bound.

## Scope

These are exact necessary constraints on the actual ORIGINAL source
signing. They do not assume an unproved optimizer property of a new
posterior, establish an order-transport upper, or close the maximum
over attainable paired shells. No numerical job or finite-order
search is part of this review.
