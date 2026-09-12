# NOTE 2026-09-12 — Exact Fourier identity for the Paley conference cube-max, and the ρ=1 criterion

**Status of this note:** small verified increment, receipt-backed; no prop number, no
status change. Written after checking `scripts/paley_conference_eval.py` (whose values it
reproduces) and reading the E1/handoff trail.

## 1. Identity (prime q ≡ 1 mod 4)

Let q be an odd prime ≡ 1 (mod 4), C the Paley conference matrix of order n = q+1
(rows indexed by F_q ∪ {∞}, C_{ab} = χ(a−b), C_{a∞} = C_{∞a} = 1, zero diagonal).
For every x ∈ {±1}^{q+1}, with

- S(x) = Σ_{a∈F_q} x_a,
- x̂(ξ) = Σ_{a∈F_q} x_a e^{−2πiξa/q},
- D(x) = Σ_{ξ=1}^{q−1} χ(ξ)·|x̂(ξ)|²,

we have the **exact identity**

    xᵀCx = 2·x_∞·S(x) + D(x)/√q.

*Numerical validation:* max abs deviation 1.4e-14 (q=13) and 2.1e-14 (q=17) over 500 random
x ∈ {±1}^{q+1}; script `/home/nick/scratch/paley_check2.py`. Full enumeration at q=13,17
reproduces the repo values Φ(C14)=21, Φ(C18)=33.

## 2. Budget form and the ρ=1 criterion

Write L(x) = Σ_{χ(ξ)=−1} |x̂(ξ)|²  (nonsquare Fourier mass), s = |S|/√q, and
σ = x_∞·sgn(S) ∈ {±1}. Since D = q² − S² − 2L, the identity becomes exactly

    xᵀCx = √q · [ q + 1 − (s − σ)² ] − 2L/√q .

Consequently:

- max_x xᵀCx ≤ (q+1)√q (the spectral bound, recovered), and
- **ρ(q) := max_x|xᵀCx| / ((q+1)√q) = 1  ⟺  ∃ x∈{±1}^{q+1} with |S(x)| = √q whose nonzero
  Fourier mass lies in a single quadratic class** (all squares, or all nonsquares; the sign
  of D then selects x_∞). CORRECTION to an earlier draft: the squares-only condition is too
  strong — the exact witnesses at q = p² (see §3b) are supported on a *nonsquare* dual line.
  Equivalently: x is a ±1 vector in the top eigenspace: Cx = ±√q·x (the repo's `Max+`,
  up to the global sign).

The second condition needs √q ∈ Z, so *exact* equality on the spectral line can only be
realized at **q = p²**, i.e. at orders n = p²+1 — precisely the E(1) family. The vector x
at equality is equivalent to a ±1 vector in the top eigenspace: Cx = √q·x, i.e. exactly the
repo's **Max+ = {y ∈ {±1}^n : Cy = py}** object. (For prime q, |S| can be taken within 1/2 of
√q, so the balance term is O(1/√q); the entire question is the **leakage term L(x)/q**: how
square-concentrated can a ±1 sequence's Fourier transform be.)

## 3. Measured ρ (full enumerations + new heuristic scan)

Full enumerations (exact) and repo values:

| q | n=q+1 | max xᵀCx | ρ = max/(n√q) | source |
|---|-------|-----------|----------------|--------|
| 13 | 14 | 42 | 0.832 | full enum (verified tonight) |
| 17 | 18 | 66 | 0.889 | full enum (verified tonight) |
| 25 | 26 | 130 | 1.000 | repo `--q 25` (Φ(C26)=65) |
| 29 | 30 | 150 | 0.9285 | repo full enum; ILS reproduces 150 |
| 37 | 38 | 218 | 0.9431 | repo full enum; ILS reproduces 218 |

New local-search (ILS) scan — best found, so these are lower bounds on max xᵀCx
(script `/home/nick/scratch/rho_scan.py`; ~10^5-10^6 evals per q):

- q = 41, 53, 61, 73, 89, 97, 101, 109, 113, 137, 149, 157, 173, 181, 193, 197:
  ρ = 0.960, 0.911, 0.913, 0.940, 0.926, 0.947, 0.915, 0.942, 0.923, 0.923, 0.906,
  0.918, 0.919, 0.932, 0.921, 0.919 (no visible trend toward 1 at fixed search budget).
- q = 293, 397 (longer runs): ρ = 0.921, 0.940.

## 3b. Exact families

- q = p²: the **cylinder** vectors x = f∘L (L an F_p-linear functional on F_{p²}, f: F_p→±1
  with Σf = 1) satisfy Cx = p·x EXACTLY. Ten-line proof: for v ∈ F_{p²},
  (Cx)_v = Σ_{c} f(c)[Σ_{w∈L^{-1}(c), w≠v}χ(v−w)] + 1; cosets not containing v give −1,
  the coset of v gives 0, so (Cx)_v = f(L(v)) if and only if Σf = 1. (Uses only
  Σ_{α∈F_p} χ(α²−c) = −1 for c ≠ 0.) Verified numerically at p = 3, 5, 7, 11, 13.
  Consequence: at every n = p²+1, `Max+` ≠ ∅ — ρ = 1 with an explicit witness, no
  enumeration needed. Verified: ILS at p = 3, 5, 7, 11, 13 returns n·p exactly.
- This is exactly the repo's k=1 stratum / "Cylinders" construction
  (`evidence/share/denseness_path_package.md` §"Cylinders", count N1 = m·C(p,m);
  matches {1:30}, {1:140}, {1:2772} at p = 5, 7, 11). The proof above is an
  independent 10-line route (coset-sum computation), so it can be used as a
  cross-check of `Lemma B`.

## 4. Why this matters

- It converts the 2^{n−1} enumeration into an O(q log q) evaluator and an exact algebraic
  criterion; the equality case is now a clean statement about ±1 sequences with Fourier
  spectrum inside the nonzero quadratic residues (GQR-code/second-order objects).
- It isolates the two loss budgets: **balance** (s−1)² and **leakage** 2L/q. E(1)-type
  statements are statements about leakage; nothing here yet proves leakage → 0 for any
  infinite prime family.
- It changes no global status: the limit question (liminf side) is untouched.

## 5. Non-goals / honesty

This note does not improve the lower bound 1/π, does not prove E(1), and does not touch
existence of lim α_n.
