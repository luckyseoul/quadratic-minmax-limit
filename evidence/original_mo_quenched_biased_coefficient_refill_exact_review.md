# Independent complete review: quenched biased coefficient refill

2026-09-05. Reviewer: exact-proof agent.

Reviewed the complete source
`/tmp/original_mo_quenched_biased_coefficient_refill.md`, SHA256
`6492829224a667deb420546643f836db838b09f194892c9b49a7c841c8bccd1b`.

**PASS. No corrections required.**

Checks performed analytically, without computational jobs:

1. The full augmented-spin posterior has edge observable
   `tau_e=sigma x_i x_j` in `{-1,1}`. Its third cumulant is exactly
   `-2r(1-r^2)`, so the global derivative bound survives every real
   intermediate Gaussian coefficient.
2. Centered absolute third moments are
   `|m_e|^3(1-rho_e^4)` for the biased sign and
   `2sqrt(2/pi)|m_e|^3(1-rho_e^2)^(3/2)` for its matched Gaussian.
   Conditional telescoping therefore proves the stated quenched
   replacement error, including zero weights and deterministic biases.
3. Every nonzero closed four-walk uses either one edge four times, two
   edges twice each, or four distinct edges once. Only the first type
   contributes a difference, with two orientations per unordered edge.
   The Gaussian fourth moment `(3-2rho_e^4)m_e^4` gives the exact positive
   correction `4 sum_e (1-rho_e^4)m_e^4`.
4. Global optimality applies outcome by outcome only to the admissible
   sign refill. Inserting the replacement and quartic identities then
   yields equation (4), with the spectral correction subtracted in the
   correct direction.
5. The row, incident-set, and all-edge estimates use respectively
   `sum_{e incident i}m_e^2=d`, `sum_{e incident S}m_e^2<=d|S|`, and
   `sum_e m_e^2=Nd/2`, together with `max_e|m_e|<=sqrt(2/N)`.
   All constants in (5)--(6) check.
6. The proof preserves the logarithm and the actual changing posterior.
   It does not evaluate the Gaussian expectation, imply a subextensive
   pressure change, or substitute the correlated canonical boundary law.

This certifies the finite independent-coefficient replacement statement,
not the unresolved order-transport inequality.
