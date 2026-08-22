# Route of attack: the class-function / character-table reduction for R1

Date: 2026-08-22.  Target: R1 = `‖δ‖² ≤ n/12` ⟺ `Φ_δ ⪰ −(2n+20)/(n−6)·I`
⟺ `λ_min(Φ) ≥ 6`, which closes leftovers 1 AND 3 (15.595).
**Nothing here is proved yet.  This is a plan plus two diagnostics.**

## Why a new route is needed (what the dead ends have in common)

Killed so far: linear 4-point LP; linear 6-point LP; level-4 SoS; Delsarte
2-design + min distance; degree escalation of the contraction kernel;
`(12+ε)n²` majorants; uniform `C/p⁴` scaling; pointwise SOS; equal-density
shortcut; Aut-invariant master + `|m₄|≤1`.  **Every one is a moment
relaxation** — it replaces Max+ by a moment sequence, and the feasible
region always contains fake sequences that meet all constraints without
being sets of ±1 eigenvectors.  Adding moments does not help: degree 6
gained nothing, established twice from independent directions (15.594 and
the 6-point box LP).  A closing argument must use something no moment
relaxation sees.

## The reformulation

For W ∈ Z, `⟨B̃_y, W⟩ = yᵀWy` (15.597), and `π(g)B̃_y = B̃_{g·y}`, so

    A_y(g) := ⟨B̃_y, π(g)B̃_y⟩ = s(y, g·y)² − 2n,     s(y,z) := y·z

is a positive-definite function on G, and the constituent energies are its
Fourier coefficients:  `‖P_c B̃_y‖² = (deg_c/|G|)·Σ_g χ_c(g)·A_y(g)`.

**Diagnostic 1 (done, decisive).**  Is the spread pointwise or only on
average?  `scripts/constituent_energy_diag.py`:

| p | per-vector `‖P_c B̃_y‖²` | verdict |
|---|---|---|
| 5 | CONSTANT across all 260 vectors (160, 288, 176 — exact integers) | pointwise holds |
| 7 | VARIES: min 0, max 1333 on the mult-100 constituent | **pointwise FAILS** |

p=5 is special (Max+ is one G-orbit there; at p=7 it cannot be, since
|G|/|Max+| = 117600/11452 ∉ ℤ).  So a single-vector argument is DEAD —
averaging over the ensemble is essential.  Total energy per y is constant
= n(n−2) at both primes, as proved (15.593 A).

## The route: Γ is a class function, so R1 is ~q numbers, not C(n,4)

Define the **ensemble correlation**

    Γ(g) := E_y[ s(y, g·y)² ] − 2n .

Then `λ_c = (1/|G|)·Σ_g χ_c(g)·Γ(g)` (using mult_c = deg_c, Z
multiplicity-free, 15.589 B).

**Key: Γ is a CLASS FUNCTION.**  Max+ is G-invariant, so for any h,
Γ(hgh⁻¹) = E_y[s(y, hgh⁻¹·y)²] = E_y[s(h⁻¹y, g·h⁻¹y)²] = Γ(g).
PSL(2,q) has ≈ q+4 conjugacy classes, so **Γ is determined by ≈ q numbers
instead of C(n,4) ≈ n⁴/24 four-set moments** — a compression by ~n³.

Then `λ_min(Φ) ≥ 6` becomes: for every irreducible χ_c appearing in Z,

    Σ_K (|K|/|G|)·χ_c(K)·Γ(K)  ≥  6 ,

a finite inequality against the **classical, fully explicit character table
of PSL(2,q)** (trivial, Steinberg, (q−3)/2 principal series, (q−1)/2
discrete series).

## Why this can work where the relaxations could not

* It uses the group action on Max+, which no moment LP encodes.
* The unknowns collapse from ~n⁴ to ~q, and the coefficients (characters)
  are known exactly in closed form for all q.
* Γ(K) plausibly depends only on the fixed-point structure of g on P¹(F_q)
  (elements are classified by trace: split / non-split / unipotent /
  central), which is exactly how Paley/Max+ structure is organised —
  so a closed form for Γ(K) is a realistic target, not a fit.
* `Γ(g) = Σ_{x,x'} d_x(g) d_{x'}(g) m₄(x, g(x), x', g(x'))` expresses Γ
  through m₄ but only via ~q aggregates, so the δ-part enters as ~q
  numbers rather than an unbounded function — the object that made every
  LP feasible-but-wrong.

## Concrete next steps, in order

1. **Compute Γ(K) at p=5, 7** for every conjugacy class, from the exact
   ensembles (cheap; both fit in memory).  Verify it is genuinely a class
   function (invariance check, not assumed).
2. **Reconstruct λ_c from Γ(K) + the PSL(2,q) character table** and check
   it reproduces the measured spectra (6.1538/11.0769/13.5385 at p=5;
   7.511/8.2152/8.9193/9.8582/10.5623 at p=7).  This validates the whole
   pipeline before any general-p claim.
3. **Hunt a closed form for Γ(K)** as a function of the class parameter
   (trace / fixed-point count).  Fit at p=5,7 then PREDICT p=11 and check
   against the stored 37.4M ensemble — out-of-sample, per the
   pre-asymptotic-fit discipline.  A form that survives p=11 is evidence;
   one that fails kills the route cheaply.
4. If a closed form holds, **prove `Σ_K (|K|/|G|)χ_c(K)Γ(K) ≥ 6`** using
   the explicit character values.  This would close R1, hence leftovers
   1 and 3 simultaneously.

## R2 (leftover 2), unchanged and independent

Max− Walsh containment, exact at p=3,5,7,11 (15.596).  Documented trap:
at p=11 `rank(B_U) = 60 < n/2 = 61`, so a general-p proof via "B_U has
full rank" is FALSE — the mechanism is the genuine algebraic containment.
Natural tool: the same signed-orbit/character-sum machinery, applied to
the fixed-edge U/U^c split rather than to four-point moments.
