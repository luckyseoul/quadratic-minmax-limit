# Unattempted directions for the three open E(1) units (read at 221e996)

Source: external review of the repo at 221e996 (2026-08-28, "Close collinear
near-line boundaries"). Live gate taken from LONG_HORIZON_GOAL.md +
e1_main_chain_status.py: Lemma D TRUE; spectral floor, residual (ii) k>=4p, and
Type I multi-level all OPEN.

What follows is what appears to be genuinely unattempted, separated from what
the failure graph already bans.

Caveat stated by the reviewer: they read the goal/status/strategy layer, the
three unit modules, the failure graph, and roughly a dozen evidence notes — not
707 source files. If any of the six items below is already killed in a note they
did not open, it belongs in E1_FAILURE_GRAPH.md, which currently stops at F22
and does not cover the boundary-census era.

## The cross-cutting one: the gate is stricter than the theorem

GOAL.md demands E(1) "for every prime p≥5." The limit does not need that.
Denseness says the global limit exists iff it exists along Paley orders; a limit
along a subsequence is unaffected by finitely many terms. E(1) is only needed
for p ≥ p₀.

This is not pedantic — it is where two of the three units are actually stuck:

- R1 is false at p=5 and p=7 and holds at p=11
  (e1_gmin_r1_principal_pge11.py, measured ‖δ‖² = 1536/65 and 19180800/1840091
  vs threshold n(λ̄−6)²/48). Meanwhile λ_min(Φ) ≥ 6 is true at both, certified
  exactly (80/13 at p=5, 7.511 at p=7, floor_lambda_min_recert.json). So the
  sufficient condition fails exactly where the conclusion is verified.
- Type I slack widens with p: μ/L = 0.769, 0.746, 0.307 at p=5,7,11.

Both units are asymptotic statements being forced into uniform form. The right
shape is "proved for p ≥ p₀" ∧ "exactly certified for 5 ≤ p < p₀," with p₀ named
in the unit. That also converts the finite work at p=5,7,11,13 from decoration
into load-bearing, which the current F22 rule correctly forbids under the
current gate. Worth confirming that Prop 6.1's interpolation doesn't secretly
need per-p control before rewiring.

## Spectral floor

Per-stratum QVAR is dead (p=41, k=7, E|Z|²=0). That does not threaten global
mixed-k QVAR: energies are nonnegative, so zero-energy strata are free. What's
missing is the stratum weight law N_k/N at general p — named "the decisive
unknown" in STRATEGY_2026-08-19_E1_remaining.md Step 4, last advanced at p=11
k=6, then dropped when effort moved to the boundary census. Nothing since has
touched it.

Two angles not found anywhere in the repo:

1. **Per-constituent character sums instead of a global L² budget.** The Hecke
   dual SOS dies because D = const − 2Σ dim(c)μ_c² is negative definite
   (coefficient −2, 15.589 D/E). That obstruction is specific to bounding ‖δ‖²
   across all isotypics at once. A per-component bound |λ_c − λ̄| ≤ 2 − ε, with
   λ_c = ⟨T,χ_c⟩/⟨χ_c,χ_c⟩ an explicit sum over PSL(2,q) on the constituents
   already classified in 15.278 F, gives the floor directly and never forms the
   L² aggregate. Weil bounds on principal-series/Weil/Steinberg constituents are
   the standard tool here.

2. **T is positive definite for free.** For fixed y, g ↦ ⟨y,π(g)y⟩ is a diagonal
   matrix coefficient, hence PD; its pointwise square is PD by Schur; T + 2n is
   an average of those. So every λ_c ≥ 0 automatically and Σ dim(c)λ_c is pinned
   by the trace. The floor is then an anti-concentration statement about
   ‖P_c(y⊗y)‖² — no isotypic is energy-deficient. 15.120 uses the Schur square on
   the Gram G⊙G but not the group-side version, which is where the constituent
   structure lives.

Step 3 of the 08-19 plan (T₁ in closed form via Jacobi/Kloosterman sums) also
remains unexecuted.

## Residual (ii)

The current route is a size ladder: size 6 closed for odd p≥5, size 8 for p=7
and finite p=11, uniform range 6 ≤ s ≤ 3(p−1)/4 for p≥17. The first survivor is
the first even s above 3(p−1)/4 — it grows linearly in p, so exact finite models
per prime cannot terminate. 15.669 §5 says as much and lists three exits; the
sessions took the most expensive one.

1. **Higher moments.** 15.669.1 uses only μ and σ² of t = |X∩B| on the middle
   Johnson slice. The parity majorant is a Markov–Krein moment problem; the third
   and fourth slice moments are exact hypergeometric expressions. A quartic
   majorant against four moments extends the feasible range directly and costs
   an afternoon, not a GPU census. This is the cheapest available widening of
   3(p−1)/4.

2. **A size-reversing involution.** Complementation is currently used only to
   fold direction sets to b ≤ (p−1)/2. The nonsquare anti-isometry transfers
   between product signs. Nothing maps large s onto small s. If a boundary
   complement exists (s ↔ p+1−s or similar), the entire open upper range
   collapses onto the closed range and the ladder ends. That's a yes/no question
   worth settling before the next size-ten push.

3. **Weil in the large-s regime.** The k-stratum emptiness criterion p > 4k²
   (15.589 M) is applied on the Max+ side but not on the boundary side. Large
   boundaries have odd-fibre counts governed by character sums over the affine
   plane; Weil deviations O(s√p) should conflict with the exact parity floors
   precisely where the two-moment argument runs out. Parity floors from below,
   Weil from above — that's the pincer that closes a range rather than a point.

## Type I

|μ| ≤ L on |κ|=1, with δ twisted-dead there, so μ = m₄⁺ and this is a single-
ensemble four-point moment — the 4-design defect the literature scan flagged as
unclosed.

1. **The box LP over the K₄-dim kernel.** 15.590 establishes that μ is an affine
   function of a K₄-dimensional δ-vector (K₄ = 1, 2, 4, 6 at p=5,7,11,13,
   computed data-free from orbits alone). The level-4 moment SDP was killed
   (PSD-feasible max|μ| = 31.2 vs true 12 at p=5), but that relaxation compresses
   to Sym²(V±) and is loose in this direction. A plain exact-rational LP
   maximizing |μ| over that low-dimensional affine space under the trivially
   valid |m₄^±(S)| ≤ 1 for every S is a different and much cheaper object,
   calibratable against the known p=5 answer. It does not appear to have been
   run.

2. **Schrijver / Terwilliger SDP.** The Terwilliger-algebra bound is built for
   exactly this — four-point relations in the Hamming scheme — and is strictly
   stronger than the level-4 moment matrix that was killed. It appears nowhere in
   the repo or the literature scan.

3. **Retarget to 2/n.** 15.191 K's |μ| ≤ 2/n implies L for p≥7 and is tight at
   p=7 (109/2863 vs 2/50, ~5%) but slack at p=11 (~30%). Same asymptotic-plus-
   finite split as above; p=5 is already a from-C theorem (15.275 L).

## Ranked

1. Rewire all three units as p ≥ p₀ + finite certification (unblocks R1, which is
   false as currently stated).
2. Quartic parity majorant on the Johnson slice — extends the uniform boundary
   range cheaply.
3. Per-constituent character-sum floor, bypassing the −2 Hecke coefficient.
4. Box LP on the K₄-dimensional kernel for Type I.
5. Resume the stratum-weight law N_k/N; global QVAR cannot be assembled without
   it.
6. Check for a size-reversing boundary involution before another size-ten census.
