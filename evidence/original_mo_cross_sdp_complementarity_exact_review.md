# Independent complete review: actual bipartite SDP complementarity

2026-09-05. Reviewer: optimized_profile_exact. Verdict: PASS.

The complete 293-line artifact
`/tmp/original_mo_cross_sdp_complementarity.md` was independently read
and checked, covering every displayed formula (1)--(19).

Approved SHA-256:
`21deee54db15dfd264106592e866bf57c2e954fadf86f84f87870ce9690ade70`.

## Independently verified points

1. The primal and diagonal dual have the stated factor of two,
   strict feasibility, and attained optima. Positivity of every
   optimal diagonal entry follows from the actual nonzero signing
   entries. Reciprocal block rescaling preserves feasibility and
   proves balance for every optimal diagonal. Schur complements,
   Cauchy--Schwarz, and completeness give exactly (4)--(5).

2. The canonical stacked matrix has unit row norms and is genuinely
   primal-feasible. Its objective is `tr|B|^3/n`, giving the stated
   nonnegative canonical gap and universal lower SDP value.

3. Both weighted diagonal distance estimates in (9), the singular
   first-moment estimate and squared-distance/tail bounds in (10),
   and the rank estimate in (11) follow with the printed constants.
   These statements assume small eta only where explicitly stated.

4. The identity `tr(Z^T Q Z)=2g` and the PSD contraction bound
   `Q D^-1 Q<=2Q` prove the cap-free constant `4g` in (13).
   The separate unweighted estimate uses `||Q||<=d_max+L_B`.
   Its polar intertwining identity remains valid at deficient rank,
   and completeness gives the exact sum-of-diagonal-variances
   identity. The factor in (15) is consequently correct.

5. At zero gap, positivity forces `QZ=0`. Transposing the second
   intertwining equation proves commutation with |B| and avoids
   any invalid inversion on its kernel. Nonzero entries of B then
   force a common scalar diagonal. Conversely, a scaled partial
   isometry gives a feasible common diagonal attaining the primal
   value. The uniqueness claim and the distinction between a
   uniform optimum and zero canonical gap are correctly stated.

6. The tensor lift has unit norm because its squared coefficients
   sum to sinh(c)=1. Its cross inner product is sin(c times the
   original inner product); the range is within the principal
   arcsine interval. Finite Gram realization makes the Gaussian
   rounding finite-dimensional and gives exactly the constant (17).

7. The singular-moment Gaussian rounding in (18) has the correct
   cubic signal and quartic error. Its self-contained bootstrap uses
   `L_B<=n` and `sum sigma^4<=L_B sum sigma^3` to obtain the positive
   coefficient `2 kappa-1`. Then `L_B^3<=sum sigma^3` gives
   `L_B=O_C(n^(5/6))=o(n)`, proving (19) uniformly from the stated
   Boolean norm cap. No additional operator-cap premise is needed.

## Scope

The artifact keeps the Boolean norm, SDP objective, canonical gap,
and original conditional optimum distinct. None of its conclusions
asserts that actual conditional optimality forces small eta or small
canonical gap. The full attainable-shell Gaussian upper and original
all-orders convergence remain unproved. No numerical job or signing
search was used in this independent review.
