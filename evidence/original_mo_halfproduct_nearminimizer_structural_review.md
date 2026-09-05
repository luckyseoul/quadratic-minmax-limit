# Independent review: leading half-product near-minimizer rigidity

2026-09-05. **PASS**, after a complete read of all 174 lines of
`/tmp/original_mo_halfproduct_nearminimizer_structural_extension.md`.

Input SHA256:
`dccc256d3b7119c666102e54cffe3a2026d31edc1bcd0c4366a15ce92c762f0f`.

The prerequisite field-response and exact-minimizer spectral proofs have
already passed this reviewer's complete independent reads, recorded in
`/tmp/original_mo_field_response_subcritical_spectral_review.md`.
The Riesz--Thorin interpolation argument in Section 4 of the canonical
Gaussian-sign information-scale note was also rechecked directly.
No computation or repository edit was used.

## Verified points

1. Random edge completion adds precisely the allowed optimality gap
   `delta_N N` to the deletion budget. The extreme-state argument then
   contributes `2 delta_N N/c` after multiplying the induced norm by
   `1/sqrt(N)`, as stated in (4).
2. The interpolation inequality `||M||_op^2 <= 16 Phi(M)` applies both
   to the complete host and every nontrivial principal signing. It yields
   the uniform coarse operator bound needed before the truncation step.
3. If normalized eigenvector l1 mass tends to zero, the displayed
   threshold leaves a support of size at most `sqrt(d_N) N` and squared
   Euclidean tail at most `d_N^(3/2)`. The coarse operator bound makes the
   Rayleigh-form error negligible. The approximate hereditary bound makes
   the retained principal operator negligible too, yielding the required
   contradiction. No additive gap is amplified by fractional rounding.
4. The previous sparse-pinning construction uses only the complete-sign
   column norms and the delocalized eigenvector. Its complement norm cap
   is eventually `cC+1`, and its positive extensive response contradicts
   `O(N^(3/4))+delta_N N=o(N)`.
5. In the sparse-energy argument, the opposite-pinning upper bound (7)
   holds uniformly. For the large-field set J, reversing the whole pinned
   block shows that the absolute cross energy is at most
   `Phi(A_(S union J))`; no extra factor two is needed. Choosing `H>b`
   in (8) gives both `|J|=o(N)` and large-field l1 mass `o(N)` uniformly.
6. If the claimed uniform medium-field sparsity failed, a sequence of
   offending supports and pins would have a fixed positive moderate-field
   density and contradict the actual-Gibbs response theorem. Thus that
   step is genuinely uniform. Splitting at a fixed `a>0`, taking limsup,
   and then sending `a` to zero proves the stated uniform l1 and Boolean
   energy conclusions.

## Conclusion and scope

Every leading half-product near-minimizer at fixed `c>0` satisfies both
`||A_N||_op=o(N^(3/4))` and uniform `o(N^(3/2))` Boolean-energy change
under deletion of any `o(N)` vertices. The proof retains the unnormalized
gap throughout. It does not identify an exact-minimizer-only property,
give `O(sqrt(N))` spectral flatness, or establish a positive-fraction
order comparison. The stated compatibility with the existing sparse
module counterfamily and the open original convergence problem is correct.
