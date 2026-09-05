# Independent complete review: conditional cross joint-shell upper

2026-09-05. Reviewer: optimized_profile_exact. Verdict: PASS.

## Frozen artifact and review scope

Approved artifact:
`/tmp/original_mo_conditional_cross_joint_shell_upper.md`, 451 lines.

SHA-256:
`64d68bb2feaa59a8049d6bcc42f3ab94c845249c3088fa618916522412d0a68a`.

The complete original 358-line core was independently read and checked
earlier, as was the subsequently added Section 7. This receipt records
those completed reviews. The current complete 451-line artifact has
also been reread and its frozen hash rechecked. Every displayed formula
(1)--(26) is covered.

## Verified implications

1. The actual cross covariance is the orthonormal symmetric-edge
   compression of the whole-source minus law. The signs, exchange
   operator `S_B(X)=B X^T B`, diagonal-one normalization, positivity,
   and operator bound `R<=3I` are correct.

2. On an attainable `(p,q,c)` shell, the stated 2n-dimensional linear
   comparison covariance is PSD. The two-dimensional compression of
   the complete source controls the norm of its nonconstant block.
   This proof specifically uses `D=L^2-1`; it does not silently claim
   positivity under a smaller arbitrary denominator.

3. The full bilinear covariance has exchange term `-d e`. The exact
   increment identity (9)--(10) is correct, including its factor
   `1/(4D)` and nonnegative square `(d-e)^2`. The independent cushion
   contributes `v(n-r_x)(n-r_y)` to the half increment excess.
   Gaussian increment comparison therefore gives a genuine upper.

4. Common marginal variances and the arcsine identity give all three
   reference energies (11), including the opposite cross mean and
   the actual independent-cushion denominator (12).

5. The spectral-midpoint/Hamming argument, subsequent Jensen step,
   and half-normal threshold optimization require no coordinate or
   selection independence. Their constants and the polynomial
   deficit (17) are correct. The optional two-block feasible set
   contains the actual expected flip fractions; its cross constraint
   follows from the stated polarization and joint concavity. Its
   optimized loss is at least the single joint mismatch loss.

6. The count `(2n^2+1)^3` safely bounds the number of nonempty shells.
   Each phase of each shell maximum has Gaussian Lipschitz constant
   at most `n sqrt(||C||)`. The deterministic means and concentration
   remainder in (20) are correct. Actual conditional optimality
   supplies (21), without assuming full order-2n optimality of K.

7. Section 7 applies the previously independently reviewed generic
   shifted-sign theorem to the actual n^2 cross coordinates with
   masked augmented observables and arbitrary deterministic-energy
   prior. It does not substitute a theorem limited to unmasked
   quadratic observables. The covariance norm constant is three.

8. Principal restriction of the complete Hermite remainder gives
   exactly (24). PSD domination plus independent Gaussian noise
   handles the operator error in both directions, at cost
   `O(n^(5/4))`. The retained rank-one term is bounded by its actual
   bilinear norm, at cost `2 n^2/(sqrt(pi) D)=O(n)`.
   These estimates are uniform over the threshold and deterministic
   energy and give (22), with error `C n^(16/11)`.

9. Every actual threshold output is an admissible cross signing
   retaining A and -A. The floor (26) therefore follows in its stated
   direction solely from conditional optimality.

## Scope of PASS

The artifact rigorously combines an actual masked Gaussian-refill
floor with an exchange-preserving Gaussian shell upper. Its remaining
sharp evaluation is explicitly open. This review does not certify
`F_A^* <= 2 sqrt(2) Phi(A)+o(n^(3/2))`, full-order optimality of the
conditional source, or original all-orders convergence.

No numerical experiment or finite-order search is part of this proof
review. Importing the same frozen bytes under a canonical path does
not change any reviewed mathematical statement.
