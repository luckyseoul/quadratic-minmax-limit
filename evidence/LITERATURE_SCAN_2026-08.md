# Literature scan: MO 413935 / lim α_n (2026-08)

## Target
Does lim α_n exist for min_{A:±1} max_{x:±1} |∑_{i<j} a_ij x_i x_j| / n^{3/2}?
Path C residual (δ² ≤ room_hyp/24 ∀ primes p≥5) or E(1) on ρ=1 Paley family.

## MO status
- https://mathoverflow.net/questions/413935 — still open (0 answers; author comments only). Crawl ~2026-08-01.
- Author (Paata Ivanisvili, UCI): liminf ≥ C>0 (mentions 2^{-5/2}); upper bounds agreed; existence open.

## Author prize / X (2026-07-25)
- @PI010101 prize post for pure-AI solution of MO 413935.
- Explicit: spent ~5 hours without solving; wants existence proof, not numerics ("I am not interested in numerical values").
- Thread: many attempts; no accepted proof in replies reviewed.

## Author recent arXiv (2025–2026) — NONE solves α_n
Sample (author search, newest first):
- 2607.25862 winding number / Fourier summation (Frank–Ivanisvili)
- 2607.19769 Weissler conjecture two-point (Huang–Ivanisvili)
- 2606.31961 Beckmann form of Talagrand discrete cube
- 2606.* Coulomb gases, Riesz transforms Hamming cube, extension norms
- 2605.05193 Grokability in five inequalities (perimeter, Hamming moments, Sidon, Szarek)
- 2605.05192 Almost-orthogonality Lp with Grok
- 2602.20462 Sharp isoperimetric Hamming cube II
- 2510.20013 NICD majority counterexample
- 2507.06115 Optimal Young convolution hypercube
- 2410.19949 Jackson inequality hypercube

Pattern: harmonic analysis / Hamming cube / inequalities / AI-assisted discovery.
**No paper titles or abstracts claim min-max ±1 quadratic form limit or conference Max+ residual.**

## Related but non-closing literature
| Work | Why relevant | Why does NOT close lim α_n |
|------|--------------|----------------------------|
| Kunisky–Yu arXiv:2211.02713 | Deg-4 SOS on Paley clique | Clique ω(G_p), not max x^T C x / n^{3/2} |
| Assumpção et al arXiv:2505.10548 | SDP fractional cut-cover on association schemes | Different parameter; not α_n |
| Krupnik thesis arXiv:2607.17829 | Conference graphs spectral determination | Not boolean quadratic max asymptotic |
| Paley DSD / conference design papers | Conference matrix construction | Statistics designs, not Path C residual |
| Little Grothendieck / maxcut SDP lit | Boolean quadratic maximization | Approximation ratios, not existence of lim min_A max_x |

## Path C residual tools in public lit
- No found closed form for Max+ weight enumerator W_k for general conference order n=p²+1.
- No found Weil/character-sum evaluation of Gegenbauer defect μ_G4 = E[Q₄(⟨z,y⟩/n)] for Max+.
- Spherical design literature treats 2-design property of equiangular lines / conference maximizers as classical; 4-design defect not closed for our residual budget.

## Verdict for implementer
**No published breakthrough closes residual or lim α_n as of this scan.**
Continuing 15.15x reformulation without a new external idea = thrash (F14).
Viable pivots if continuing:
1. E(1) structural argument on ρ=1 family (still open; n=10 matching classification does not lift).
2. Import a specific character-sum / association-scheme theorem if one is identified for conference Max+ — none found here.
3. Honest permanent gap only if disproved — objective is totality close, not gap report.

## Audit
- MO page, author site, arXiv author list (68 papers), X prize thread, targeted web/S2 (S2 rate-limited).
- Date: 2026-08 (session).

## 2026-08-26 post-finding check: Proposition 15.655

After the mod-seven catalog-syndrome exclusion was found, targeted searches
for its exact dimensions/rank and Paley--Johnson description found no
matching theorem. Ghinelli--Key (2011), DOI `10.3934/amc.2011.5.93`, studies
`p`-ary codes from ordinary Paley graph and line-graph incidence matrices;
it is relevant finite-field coding context but does not contain the
`282 x 1225` affine-score matrix, rank 147 over `F_7`, Johnson slack right
sides, or the 2,408-case exclusion. Direct OEIS searches for `1716742440`,
`3939012`, `1372,294,112,294,336`, and `282,1225,147,135` returned no
results. No sequence submission or broad priority claim is proposed.

## 2026-08-26 post-finding check: Proposition 15.656

Searches for `132 x 325`, ranks `67/65`, the combined tuple
`262,325,113,149`, and a Paley full-eigenshell bounded-syndrome theorem found
no matching result. Ghinelli--Key's ordinary Paley incidence codes are
nearby finite-field coding context but use different matrices. Individual
OEIS searches locate `26450`, `15525`, and `10925` in unrelated partition
sequences; no relevant match was found for the structural rank or coverage
tuples. This is a duplication check, not a sequence-submission project.

The requested paper arXiv:2305.03523 was read. It constructs minimal locally
concave Bellman functions on planar domains formed as a difference of
unbounded convex domains. That mechanism may suggest a two-moment envelope
for R1/global QVAR, but the paper has no finite-Paley theorem and does not
contribute to the 15.656 proof.
