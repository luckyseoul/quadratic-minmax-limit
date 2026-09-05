# Current mathematical status

Updated 2026-09-05. The original MO limit is OPEN.

There is no reviewed proof of convergence or nonconvergence. The value
`1/2` is unproved. These statements are independent of the status of the
optional Paley research program.

## What is established

`CORE.md` gives the route-neutral definitions and proofs. In particular,
`m_n` is nondecreasing, consecutive `alpha_n` differences tend to zero,
the limit-point set is an interval, and

`1/pi <= liminf alpha_n <= limsup alpha_n <= 1/2`.

The upper bound uses conference constructions, not their optimality.
Ratio-dense transfer and two-multiplier Dini amplification are valid
conditional tools; neither supplies its missing hypothesis automatically.

## What changed at the reset

Residual (ii) is parked as an optional route-local open lemma. Its former
status as an unavoidable obstacle was a bookkeeping error. The same applies
to treating an E1/bi-tight conjunction as an if-and-only-if test for the
original limit. The gap-two condition is stronger than the asymptotic
optimality needed to transfer the value `1/2`; neither has been proved
necessary for convergence.

The active global proof registry is `src/original_mo_status.py`.
`src/e1_main_chain_status.py` retains separate optional-route diagnostics.
Historical proof notes, scoped local theorems, and exact certificates remain
available; unproved bridges have not been promoted or declared false.
See `ARTIFACTS.md` for the preserved branches and terminology.

## Fresh mathematical target

The [induced-restriction theorem](evidence/NOTE_2026-09-05_INDUCED_OPTIMIZER_RESTRICTIONS.md)
is proved and independently reviewed. For `n -> infinity` with
`n^2 = o(log N)`, every `n`-vertex signing occurs in every source signing
of norm `O(N^(3/2))`. Nevertheless, a uniform restriction typically has
normalized norm at least `(2/3)*sqrt(2/pi) = 0.531923...`; for globally
optimal sources, `exp(o(n))` uniform-marginal samples cannot preserve the
source's leading constant with probability bounded away from zero.

This proves that smaller-order optimal restrictions exist, not that their
constant matches the source. The latter remains precisely an unproved
cross-order comparison. The theorem excludes only the stated sampling
method and scale, not biased selection or other order comparisons.

Compare different orders using the actual global-minimizer property,
without assuming a Paley model, a limiting value, or a prescribed lift.
Any claimed convergence proof must control the normalized optimum, not just
an arbitrary low-norm signing or a typical induced restriction.

No all-orders comparison with sufficient error control is established.
Finite exclusions and moment identities do not change that status.
