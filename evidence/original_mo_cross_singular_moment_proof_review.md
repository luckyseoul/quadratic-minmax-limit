# Independent proof review: actual cross singular moments

2026-09-05. The proof agent read the complete 168-line source
`/tmp/original_mo_cross_singular_moment_rounding.md`, SHA-256
`6d5129a1572842c76c8f11a008b0093cb3c340684a40219b7db8828fdeeaf756`.

Result: **PASS**, with no corrections requested.

The general filtered Gram construction is valid: the rows of
`B/sqrt(n)` are unit vectors, the columns of `f(|B|)/sqrt(n)` have
norm at most one by PSD functional calculus, and separate orthogonal
padding preserves all cross inner products. The endpoint-valid arcsine
remainder follows from the positive Taylor coefficients. Both the
signal and Frobenius-error trace identities are exact.

The cubic and clipped corollaries have the correct constants. In
particular clipping bounds the error by `(1-2/pi) K^2 n` without
an operator bound or any change to the source entries. The normalized
constraint (8) has the correct `q/sqrt(n)` finite-order coefficient.

The active-shell restriction is correctly explicit: `beta(B)` can be
replaced by an attained shell magnitude only when that magnitude is
actually `beta(B)`. Conditional optimality alone does not permit that
replacement at every shell. The exact partial-isometry case correctly
has vector SDP value `nd` and asymptotic ratio at most `pi/2` under
the bounded-operator scale; this is not promoted to a universal ratio.

This review is analytic. No computational run was used. The source
does not claim to evaluate the complete conditional Gaussian maximum
or settle original convergence, and neither does this receipt.
