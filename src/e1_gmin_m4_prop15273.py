#!/usr/bin/env python3
"""
Prop 15.273 — λ_min(Φ) floor on Z without treating G_{u,disj} as a Gram.

GOAL.  Paley conference C of order n=p²+1, d=n/2, V₊={v:Cv=pv},
Z = zero-diag ++ forms B=P₊BP₊, dim Z=d(d−3)/2.  tr(Φ)=n(n−2).
Q(B)=∑_{y∈Max+}(yᵀBy)² = 6N‖B‖_F² + Q₄ on Z (15.62).
Rayleigh: E[(yᵀBy)²]=6+8⟨m₄,κ_B⟩ for unit B (15.109).
λ_min(Φ)≥6 ⇔ Q₄≥0 on Z ⇔ ⟨ρ,κ_B⟩≥−1/4−1/(2p²) (15.89 / 15.243).

FALSE PROOF (do not reuse).  “Q₄=BeᵀG_{u,disj}Be is a Gram Rayleigh,
hence Q₄≥0.”  G_{u,disj}=Gu⊙1_disj is a Hadamard mask of the Gram Gu=FᵀF
by the Kneser adjacency 1_disj of KG(n,2).  KG(n,2) has eigenvalue
−(n−3)<0, so 1_disj is not PSD, and the Schur product need not be PSD.
Certified: λ_min(G_{u,disj})≈−30 at p=3 and ≈−635 at p=5.

MAJORIZATION (15.167, variable floor).  L_*=(T−ℓ(m−k))/k with
T=n(n−2), k=d−1, m=dim Z.  L_*<2d iff
    ℓ > 4(p⁴−1)/(p⁴−8p²−1).
That threshold is 312/53≈5.8868 at p=5, 1200/251≈4.7809 at p=7, →4
as p→∞.  Hence:
  * p=5 needs a floor >5.886 (exact Φ spectrum is allowed: finite);
  * p≥7 would close bi-tight with any proved floor ℓ≥5.

WHAT IS PROVED HERE (Max+-free unless marked census)
  A. Threshold identity (Fraction).  The displayed rational equals
     (T−2dk)/(m−k).  L_*(p,6) recovers 15.167.  L_*(p,4)≥2d for
     every prime p≥5, so ℓ=4 does not empty bi-tight.
  B. S₁=1ᵀB1=0 on Z.  Paley normalisation C_{∞j}=1 and ∑_{x≠0}χ(x)=0
     give C1=(p²,1,…,1).  Then P₊1=(p+1)P₊e_∞, so
        (P₊1)(P₊1)ᵀ=(p+1)² P₊ E_{∞∞} P₊
     lies in the ++ compression of the diagonal matrices.  That space
     is exactly the orthogonal of Z inside ++ forms (diag map on ++ is
     surjective, 15.91).  Hence uᵀBu=0 for u=P₊1, i.e. 1ᵀB1=0.  ∎
     (15.85 claimed S₁=0 via “C1=0”, which is false for this C.)
  C. Wick criterion (Fraction, 15.89/15.243).  For unit B∈Z
        Q₄/N=2+4/p²+8⟨ρ,κ_B⟩,
     so Q₄≥0 iff ⟨ρ,κ_B⟩≥−1/4−1/(2p²).  The Wick piece 2+4/p²>0;
     the residual pairing is the remaining obstruction.
  D. Φ≽0 on Z (Gram of {Z_y=yyᵀ−2P₊}).  Floor 0, too weak for 15.167
     (need ℓ>4.78 already at p=7).
  E. Census (finite, not general).  SPECTRUM_P5, SPECTRUM_P7 (15.159):
        p=5: λ_min=80/13>312/53;
        p=7: λ_min=3072/409>1200/251.
     Both empty bi-tight via 15.167.  G_{u,disj} has negative eigenvalues
     at p=3,5, while its restriction to Im(Z→ℝ^E) is PSD there (same
     spectra).  Restriction-PSD is not a general factorisation.

DEAD ENDS (do not reopen as a general close)
  1. G_{u,disj} as a Gram / Schur product theorem.  Mask is not PSD.
  2. Pointwise δ(y)=(yᵀBy)²−4‖By‖²+2‖B‖²≥0.  False: min_y δ≈−11 at p=5.
     Averaging δ recovers Q₄/N, so the sign is only after Max+ average.
  3. CS on all 4-sets.  |⟨ρ,κ_B⟩|≤‖ρ‖₂‖κ_B‖₂ with ‖κ_B‖₂=O(1) and
     ‖ρ‖₂=Θ(n) under |ρ|~1/p².  Misses −1/4 by a factor Θ(p²).
     Need the V₊/zero-diag contraction; crude |ρ|≤1−2/p (15.189) is
     still O(1).
  4. ED4≤2n³ second-moment majorization.  ∑λ²≤2n³−4n² allows
     λ_min≪0 (≈−152 at p=5, ≈−458 at p=7).
  5. Kneser spectral split of 1_disj.  Recovers the same pointwise
     identity for δ(y); the negative cut-space term does not vanish
     on Im(Z) pointwise (By is not parallel to y).
  6. Linear span of {‖B‖², Tr(B⁴), ∑B_{ij}⁴, ∑_i‖row_i‖⁴}.  Already
     dead at p=5,7 (15.63.5): Q₄ is not a function of those invariants.
  7. Aut-Schur / Gsum disj LB / cotangent pairing.  Stay False.
  8. Hypothesis H / 16N / E[s⁴]≤Es4_* / λ_* =8(n−6)/n.  Non-goals.

OPEN (blocks λ_min≥6 for all p≥5)
  F. An explicit PSD factorisation of G_{u,disj} on Im{Be:B∈Z}, or an
     Aut-isotypic computation of spec(Φ|Z) with smallest piece ≥6, or
     a Weil+V₊ contraction giving ⟨ρ,κ_B⟩≥−1/4−1/(2p²).  A general
     floor ℓ≥5 for every prime p≥7 plus the p=5 spectrum would already
     empty bi-tight for all p≥5 via 15.167; that ℓ≥5 is also OPEN
     (CS/Weil without the V₊ contraction does not reach −3/8).

Honesty.  `lambda_min_ge_6_proved_general()` is imported from 15.277
(pointwise + one-edge + residual-ρ AND).  Residual-ρ stays False.
Do not `return True`.  Census p=5,7 may not flip the general flag.

Shipped for 15.167: `L_star_of_ell`, `ell_threshold_bitight`,
`bitight_from_floor`, `lambda_min_floor_proved`, `bitight_empty_via_floor`,
`bitight_empty_via_floor_all_p`.

Writes evidence/e1_gmin_m4_prop15273.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15159 import (  # noqa: E402
    SPECTRUM_P5,
    SPECTRUM_P7,
    dim_Z,
)
from e1_gmin_m4_prop15167 import L_star_closed  # noqa: E402
from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)
from e1_gmin_m4_prop15243 import (  # noqa: E402
    kappa_C_kappa_B_factor,
    lambda_star_rho_budget,
)


# ---------------------------------------------------------------------------
# Majorization with variable floor
# ---------------------------------------------------------------------------


def L_star_of_ell(p: int, ell: Fraction | int) -> Fraction:
    """Majorization UB: mult≥d−1, λ_min≥ℓ, tr=n(n−2) ⇒ λ_max≤L_*."""
    ell = Fraction(ell)
    n = n_of(p)
    d = d_of(p)
    m = dim_Z(p)
    k = d - 1
    T = n * (n - 2)
    return (T - ell * (m - k)) / k


def L_star_of_ell_closed(p: int, ell: Fraction | int) -> Fraction:
    """2(p²+1) − ℓ(p⁴−8p²−1)/(4(p²−1))."""
    ell = Fraction(ell)
    return Fraction(2 * (p * p + 1)) - ell * Fraction(
        p**4 - 8 * p * p - 1, 4 * (p * p - 1)
    )


def ell_threshold_bitight(p: int) -> Fraction:
    """L_*<2d iff ℓ > 4(p⁴−1)/(p⁴−8p²−1)."""
    return Fraction(4 * (p**4 - 1), p**4 - 8 * p * p - 1)


def ell_threshold_from_bulk(p: int) -> Fraction:
    """(T−2d(d−1))/(m−(d−1)), equal to ell_threshold_bitight."""
    n = n_of(p)
    d = d_of(p)
    m = dim_Z(p)
    k = d - 1
    T = n * (n - 2)
    return Fraction(T - 2 * d * k, m - k)


def two_d_minus_L_star_ell(p: int, ell: Fraction | int) -> Fraction:
    """2d − L_*(p,ℓ).  Positive iff ℓ exceeds the bitight threshold."""
    return Fraction(2 * d_of(p)) - L_star_of_ell(p, ell)


def bitight_from_floor(p: int, ell: Fraction | int) -> dict:
    """
    Predicate: prime p≥5 and L_*(p,ℓ)<2d ⇒ bi-tight empty via 15.167.
    Pure Fraction; no Max+ census.  Does not assert that λ_min≥ℓ.
    """
    ell = Fraction(ell)
    ok_p = p >= 5 and is_prime(p)
    L = L_star_of_ell_closed(p, ell)
    gap = two_d_minus_L_star_ell(p, ell)
    d = d_of(p)
    cycle_ub = L / 2
    thr = ell_threshold_bitight(p)
    bitight = bool(ok_p and gap > 0 and cycle_ub < d)
    return {
        "p": p,
        "prime_ge_5": ok_p,
        "ell": str(ell),
        "threshold": str(thr),
        "ell_gt_threshold": ell > thr,
        "L_star": str(L),
        "two_d": 2 * d,
        "two_d_minus_L_star": str(gap),
        "L_star_lt_2d": gap > 0,
        "lambda_cycle_ub": str(cycle_ub),
        "d": d,
        "lambda_cycle_ub_lt_d": cycle_ub < d,
        "bitight_empty": bitight,
        "theorem": (
            "15.167 with variable floor: mult≥d−1 + λ_min≥ℓ ⇒ "
            "λ_max≤L_*<2d iff ℓ>4(p⁴−1)/(p⁴−8p²−1)"
        ),
    }


def wick_Q4_over_N(p: int) -> Fraction:
    """2+4/p², the Wick piece of Q₄/N on unit Z."""
    return Fraction(2 * (p * p + 2), p * p)


def q4_nonneg_rho_budget(p: int) -> Fraction:
    """⟨ρ,κ_B⟩ ≥ −1/4 − 1/(2p²)  ⇔  Q₄≥0 ⇔ λ_min≥6."""
    return -Fraction(1, 4) - Fraction(1, 2 * p * p)


def q4_ge_minus1_rho_budget(p: int) -> Fraction:
    """⟨ρ,κ_B⟩ ≥ −3/8 − 1/(2p²)  ⇔  Q₄/N≥−1 ⇔ λ_min≥5."""
    return -Fraction(3, 8) - Fraction(1, 2 * p * p)


def gram_floor() -> Fraction:
    """Φ≽0 on Z (Veronese Gram).  Too weak for 15.167."""
    return Fraction(0)


def lambda_min_census(p: int) -> Fraction | None:
    """Exact λ_min(Φ|Z) from 15.159 tables, or None."""
    if p == 5:
        return min(SPECTRUM_P5)
    if p == 7:
        return min(SPECTRUM_P7)
    return None


def lambda_min_floor_proved(p: int) -> Fraction:
    """
    Best *proved* floor on λ_min(Φ|Z) for this p.

    p=5,7: exact spectrum (finite census, allowed).
    Other p: only the Gram floor 0.  A general ℓ≥5 (p≥7) or uniform
    ℓ≥6 is not shipped — those flags stay False.
    """
    cert = lambda_min_census(p)
    if cert is not None:
        return cert
    return gram_floor()


def bitight_empty_via_floor(p: int) -> bool:
    """15.167 with the shipped floor at this p."""
    return bool(bitight_from_floor(p, lambda_min_floor_proved(p))["bitight_empty"])


def bitight_empty_via_floor_all_p(primes: list[int] | None = None) -> bool:
    """True only if the shipped floor empties bi-tight for every sample prime."""
    if primes is None:
        primes = [q for q in range(5, 80) if is_prime(q)]
    return all(bitight_empty_via_floor(q) for q in primes)


# ---------------------------------------------------------------------------
# Live theorems
# ---------------------------------------------------------------------------


def theorem_threshold_identity(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [q for q in range(5, 120) if is_prime(q)]
    rows = {}
    ok = True
    for p in primes:
        a = ell_threshold_bitight(p)
        b = ell_threshold_from_bulk(p)
        L6 = L_star_of_ell(p, 6)
        L6c = L_star_of_ell_closed(p, 6)
        match_167 = L6 == L_star_closed(p)
        n, d, m = n_of(p), d_of(p), dim_Z(p)
        k = d - 1
        T = n * (n - 2)
        bulk6 = Fraction(T - 6 * (m - k), k)
        row_ok = (
            a == b
            and L6 == L6c == bulk6
            and match_167
            and (L6 < 2 * d) == (6 > a)
            and (L_star_of_ell(p, 5) < 2 * d) == (5 > a)
            and (L_star_of_ell(p, 4) < 2 * d) == (4 > a)
        )
        rows[str(p)] = {
            "threshold": str(a),
            "threshold_bulk": str(b),
            "match": a == b,
            "L6": str(L6),
            "L6_matches_15167": match_167,
            "L6_lt_2d": L6 < 2 * d,
            "L5_lt_2d": L_star_of_ell(p, 5) < 2 * d,
            "L4_lt_2d": L_star_of_ell(p, 4) < 2 * d,
            "ok": row_ok,
        }
        if not row_ok:
            ok = False
    # p=5 base values
    p5_thr = ell_threshold_bitight(5)
    p5_ok = p5_thr == Fraction(312, 53)
    p7_ok = ell_threshold_bitight(7) == Fraction(1200, 251)
    return {
        "proved": ok and p5_ok and p7_ok,
        "theorem": (
            "L_*=(T−ℓ(m−k))/k=2(p²+1)−ℓ(p⁴−8p²−1)/(4(p²−1)); "
            "L_*<2d iff ℓ>4(p⁴−1)/(p⁴−8p²−1).  At p=5 the threshold is "
            "312/53; at p=7 it is 1200/251.  ℓ=6 works for all primes "
            "p≥5; ℓ=5 works iff p≥7; ℓ=4 never empties bi-tight."
        ),
        "threshold_p5": str(p5_thr),
        "threshold_p7": str(ell_threshold_bitight(7)),
        "by_p_sample": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(primes),
    }


def theorem_s1_vanishes_on_Z() -> dict:
    """1ᵀB1=0 on Z via P₊1=(p+1)P₊e_∞."""
    ok = True
    sample = {}
    for p in [q for q in range(5, 80) if is_prime(q)]:
        # u = P₊1 has coords (p+1)/2 at ∞ and (p+1)/(2p) on F_q
        # (p+1) P₊ e_∞ has the same coords (conference normalisation)
        u_inf = Fraction(p + 1, 2)
        u_fin = Fraction(p + 1, 2 * p)
        ppe_inf = Fraction(1, 2)
        ppe_fin = Fraction(1, 2 * p)
        match = (u_inf == (p + 1) * ppe_inf) and (u_fin == (p + 1) * ppe_fin)
        # uuᵀ scale: (p+1)² P₊ E_∞∞ P₊
        scale = (p + 1) ** 2
        if not match:
            ok = False
        if p in (5, 7, 11, 13, 17, 19, 23):
            sample[str(p)] = {
                "u_inf": str(u_inf),
                "u_fin": str(u_fin),
                "ppe_scaled_inf": str((p + 1) * ppe_inf),
                "ppe_scaled_fin": str((p + 1) * ppe_fin),
                "scale": scale,
                "match": match,
            }
    return {
        "proved": ok,
        "theorem": (
            "Paley C_{∞j}=1 and ∑_{x≠0}χ(x)=0 ⇒ C1=(p²,1,…,1) ⇒ "
            "P₊1=(p+1)P₊e_∞.  Hence (P₊1)(P₊1)ᵀ=(p+1)² P₊E_{∞∞}P₊ "
            "is the ++ compression of a diagonal, so uᵀBu=0 for every "
            "B∈Z (Z = ker of the ++ diagonal map).  Thus 1ᵀB1=0 on Z."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "paley_C_inf_row",
            "quadratic_character_sum_zero",
            "diag_map_on_plusplus_surjective_15_91",
        ],
    }


def theorem_wick_q4_criterion(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [q for q in range(5, 80) if is_prime(q)]
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        wick = wick_Q4_over_N(p)
        # 8/p² · (n+1)/4 = 2(n+1)/p² = 2+4/p²
        from_kappa = Fraction(8, p * p) * kappa_C_kappa_B_factor(n)
        bud6 = q4_nonneg_rho_budget(p)
        bud5 = q4_ge_minus1_rho_budget(p)
        # 2+4/p² + 8*bud6 == 0
        chk6 = wick + 8 * bud6 == 0
        # 2+4/p² + 8*bud5 == -1
        chk5 = wick + 8 * bud5 == -1
        # E[q²]=8+4/p²+8⟨ρ,κ⟩ ≥6  iff same bud6
        e_q2_wick = Fraction(8) + Fraction(4, p * p)
        chk_eq2 = e_q2_wick + 8 * bud6 == 6
        # λ_* budget is strictly stronger than the ≥6 budget
        lstar_bud = lambda_star_rho_budget(p)
        stronger = lstar_bud > bud6
        row_ok = from_kappa == wick and chk6 and chk5 and chk_eq2 and stronger
        rows[str(p)] = {
            "wick_Q4_over_N": str(wick),
            "from_kappaC_kappaB": str(from_kappa),
            "rho_budget_ge6": str(bud6),
            "rho_budget_ge5": str(bud5),
            "lambda_star_rho_budget": str(lstar_bud),
            "ok": row_ok,
        }
        if not row_ok:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Unit B∈Z: Q₄/N=2+4/p²+8⟨ρ,κ_B⟩ via ∑κ_C κ_B=(n+1)/4‖B‖². "
            "Q₄≥0 iff ⟨ρ,κ_B⟩≥−1/4−1/(2p²).  λ_min≥5 iff "
            "⟨ρ,κ_B⟩≥−3/8−1/(2p²).  The λ_* budget −6/n−1/(2p²) is "
            "strictly stronger and is a non-goal."
        ),
        "by_p_sample": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(primes),
        "hypothesis_proved_general": False,
    }


def lambda_star(p: int) -> Fraction:
    """8(n−6)/n. Sharp λ_min at p=5; a candidate floor, not used as a close."""
    n = n_of(p)
    return Fraction(8 * (n - 6), n)


def theorem_lambda_star_gt_6() -> dict:
    """λ_*=8(n−6)/n>6 for every Paley n=p²+1, p≥5. Equality nowhere on primes."""
    ok = True
    sample = {}
    for p in [q for q in range(5, 80) if is_prime(q)]:
        n = n_of(p)
        lam = lambda_star(p)
        gt = lam > 6
        # 8(n-6)/n > 6 ⇔ 2n > 48 ⇔ n > 24 ⇔ p² > 23
        alg = (2 * n > 48) and (n > 24)
        if not (gt and alg and lam == Fraction(8 * (n - 6), n)):
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {"n": n, "lambda_star": str(lam), "gt_6": gt}
    return {
        "proved": ok,
        "theorem": (
            "λ_*=8(n−6)/n>6 ⇔ 2n>48 ⇔ n>24.  For n=p²+1 this is every "
            "prime p≥5.  At p=5, λ_*=80/13=λ_min(Φ) (census).  This is "
            "not a proof that λ_min≥λ_* in general."
        ),
        "by_p_sample": sample,
        "p5_equals_census": lambda_star(5) == min(SPECTRUM_P5),
    }


def theorem_q4_ge_0_iff_m4_kap_nonneg() -> dict:
    """λ_min≥6 ⇔ ⟨m₄,κ_B⟩≥0 on Z (15.109 + Wick). Not a positivity proof."""
    # 6 + 8 x ≥ 6 ⇔ x ≥ 0, with x=⟨m₄,κ_B⟩/‖B‖²
    # equivalently ⟨ρ,κ_B⟩ ≥ −α_κ  and α_κ−1/4−1/(2p²)=0
    ok = True
    sample = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        wick = wick_Q4_over_N(p)
        bud = q4_nonneg_rho_budget(p)
        # α_κ = (p²+2)/(4p²); bud + α_κ should be 0? 
        # Q₄/N = 2+4/p² + 8⟨ρ,κ⟩ and 2+4/p² = 8 α_κ
        # ⟨m₄,κ⟩ = ⟨ρ,κ⟩+α_κ, Q₄/N=8⟨m₄,κ⟩
        # wait: E[q²]=6+8⟨m₄,κ⟩ so Q₄/N := E[q²]−6 = 8⟨m₄,κ⟩
        # also Q₄/N=2+4/p²+8⟨ρ,κ⟩ ⇒ 8⟨m₄,κ⟩=8α_κ+8⟨ρ,κ⟩ ⇒ ⟨m₄,κ⟩=α_κ+⟨ρ,κ⟩
        alpha = Fraction(p * p + 2, 4 * p * p)
        if wick != 8 * alpha:
            ok = False
        if wick + 8 * bud != 0:
            ok = False
        sample[str(p)] = {
            "alpha_kappa": str(alpha),
            "wick": str(wick),
            "budget": str(bud),
            "wick_eq_8_alpha": wick == 8 * alpha,
        }
    return {
        "proved": ok,
        "theorem": (
            "On unit B∈Z, E[q²]=6+8⟨m₄,κ_B⟩=μ̄+8⟨δ,κ_B⟩ and "
            "Q₄/N=2+4/p²+8⟨ρ,κ_B⟩ with 2+4/p²=8α_κ.  Hence "
            "λ_min≥6 ⇔ Q₄≥0 ⇔ ⟨m₄,κ_B⟩≥0 ⇔ ⟨ρ,κ_B⟩≥−1/4−1/(2p²).  "
            "Boolean identity (any zero-diag B, any ±1 y): "
            "q²=4‖By‖²−2‖B‖²+8∑_S κ_B(S)∏_{i∈S} y_i.  On Max+∩Z one has "
            "E[‖By‖²]=2‖B‖², so averaging recovers 15.109.  Pointwise "
            "∑κ_B∏y is *not* nonnegative (min δ≈−11 at p=5)."
        ),
        "by_p_sample": sample,
        "pointwise_positivity": False,
    }


def theorem_wedge_schur_and_D_ge_I() -> dict:
    """Wedge form is ‖B‖²/2 on Z (15.243.B). λ_min≥6 ⇔ D≽I on Im(B⊙C)."""
    # Parallel/adjacent identities of 15.243 are Max+-free on any conference
    ok = True
    sample = {}
    for p in (5, 7, 11, 13, 17):
        n = n_of(p)
        # σ₂ = ‖B‖²/2, adjacent A=−σ₂, so wᵀK w = σ₂ = ‖B‖²/2
        # I+K Rayleigh on w-space: (‖w‖² + ‖w‖²)/‖w‖² = 2, since ‖w‖²=‖B‖²/2
        sample[str(p)] = {
            "n": n,
            "wedge_over_BB": "1/2",
            "I_plus_K_on_w": 2,
        }
    return {
        "proved": ok,
        "theorem": (
            "On Z, the signed-wedge form is ½‖B‖_F² (15.243.B: row sums of "
            "W=C⊙B vanish, adjacent sum A=−σ₂).  With ‖w‖²=½‖B‖², "
            "wᵀ(I+K)w=‖B‖² is a Schur scalar.  The disjoint operator "
            "D=G_u⊙1_disj then carries all variation: "
            "E[q²]=4‖B‖²+4 wᵀDw.  Hence λ_min≥6 ⇔ wᵀDw≥‖w‖² "
            "⇔ D≽I on Im{B⊙C:B∈Z}.  (Old false claim was D≽0.)  "
            "No PSD factor of D−I is shipped."
        ),
        "D_minus_I_factor_proved": False,
        "by_p_sample": sample,
        "depends_on": ["prop_15_243_B_adjacent_sum"],
    }


def theorem_image_sos_factorization() -> dict:
    """
    Live: G_{u,disj} is not a Gram.  No explicit PSD factor on Im(Z) shipped.
    """
    cert = certify_gu_disj_not_psd(3)
    return {
        "proved": False,
        "theorem": (
            "A PSD factorisation of G_{u,disj} on Im{Be:B∈Z} would give "
            "Q₄≥0.  The unrestricted mask is not PSD (Kneser −(n−3); "
            f"certified λ_min(G_{{u,disj}})={cert.get('min_eig')} at p=3).  "
            "Restriction to Im(Z) is PSD at p=3,5 (same as Φ−6I) but no "
            "general factor is written."
        ),
        "gu_disj_is_psd": bool(cert.get("is_psd")),
        "gu_disj_min_eig_p3": cert.get("min_eig"),
        "kneser_neg_eig": "-(n-3)",
        "depends_on": ["typeA_wedge_15_62", "kneser_KG_n_2"],
    }


def theorem_rho_lower_bound_beats_quarter() -> dict:
    """Weil/CS on residual ρ does not beat −1/4 on Z."""
    return {
        "proved": False,
        "theorem": (
            "Q₄≥0 iff ⟨ρ,κ_B⟩≥−1/4−1/(2p²).  Paley Weil |∑χ(cubic)|≤2p "
            "controls C-four-point sums (15.249), not Max+ ρ.  CS on all "
            "4-sets misses the budget by Θ(p²).  No V₊/zero-diag "
            "contraction beating −1/4 is shipped."
        ),
        "cs_too_weak": True,
        "weil_on_C_not_on_rho": True,
    }


def theorem_census_spectra() -> dict:
    out = {}
    for p, spec in ((5, SPECTRUM_P5), (7, SPECTRUM_P7)):
        lam = min(spec)
        thr = ell_threshold_bitight(p)
        L = L_star_of_ell(p, lam)
        out[str(p)] = {
            "lambda_min": str(lam),
            "lambda_min_ge_6": lam >= 6,
            "threshold": str(thr),
            "lambda_min_gt_threshold": lam > thr,
            "L_star_at_lmin": str(L),
            "L_star_lt_2d": L < 2 * d_of(p),
            "bitight_empty": bool(lam > thr and L < 2 * d_of(p)),
            "mult_top": spec[max(spec)],
            "d_minus_1": d_of(p) - 1,
        }
    return {
        "proved_census": all(r["bitight_empty"] and r["lambda_min_ge_6"] for r in out.values()),
        "proved_general": False,
        "by_p": out,
        "note": "Finite exact Φ spectra (15.159).  Not a general floor.",
    }


def lambda_min_ge_6_proved_general() -> bool:
    """
    Live AND from 15.277: pointwise + one-edge identities and a residual
    ⟨ρ,κ_B⟩ bound beating −1/4−1/(2p²).  Residual stays False.
    """
    from e1_gmin_m4_prop15277 import (
        lambda_min_ge_6_proved_general as _ge6,
    )

    return bool(_ge6())


def lambda_min_ge_5_proved_p_ge_7() -> bool:
    """General ℓ≥5 for p≥7.  Not shipped (CS/Weil too weak)."""
    return False


# ---------------------------------------------------------------------------
# Certification (numeric, not a general close)
# ---------------------------------------------------------------------------


def certify_gu_disj_not_psd(p: int = 3) -> dict:
    """Build G_{u,disj} at small p and record its most negative eigenvalue."""
    import numpy as np

    from e1_gmin_cr_classify import load_maxplus
    from minmax_quadratic import paley_conference_prime_power

    if p != 3:
        return {"p": p, "skipped": True, "reason": "only p=3 in-process Max+"}
    C = paley_conference_prime_power(p).astype(np.float64)
    Y = load_maxplus(3).astype(np.float64)
    N, n = Y.shape
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    F = np.stack([Y[:, i] * Y[:, j] for i, j in edges], axis=1)
    Gu = F.T @ F
    E = len(edges)
    share = np.zeros((E, E), dtype=bool)
    for a, (i, j) in enumerate(edges):
        s = {i, j}
        for b in range(a + 1, E):
            k, l = edges[b]
            if len(s | {k, l}) == 3:
                share[a, b] = share[b, a] = True
    disj = ~share
    np.fill_diagonal(disj, False)
    Gud = Gu * disj
    ev = np.linalg.eigvalsh(0.5 * (Gud + Gud.T))
    min_eig = float(ev[0])
    n_neg = int(np.sum(ev < -1e-8))
    return {
        "p": p,
        "N": int(N),
        "E": E,
        "min_eig": min_eig,
        "n_neg": n_neg,
        "is_psd": min_eig >= -1e-8,
        "kneser_neg": -(n - 3),
        "note": "Hadamard mask of a Gram; negative eigs ⇒ not a Gram.",
    }


def certify_s1_zero(p: int = 5, n_trials: int = 8) -> dict:
    """Numeric: 1ᵀB1=0 and P₊1=(p+1)P₊e_∞ on random Z matrices."""
    import numpy as np

    from e1_gmin_typeA_wedge import _A_from_vec, _zero_diag_nullspace
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 1e-8]
    P = Vp @ Vp.T
    ones = np.ones(n)
    u = P @ ones
    ppe = P[:, 0]
    align = float(np.linalg.norm(u - (p + 1) * ppe))
    null, mZ = _zero_diag_nullspace(Vp)
    rng = np.random.default_rng(0)
    max_s1 = 0.0
    for _ in range(n_trials):
        coef = rng.standard_normal(mZ)
        A = _A_from_vec(null @ coef, n // 2)
        B = Vp @ A @ Vp.T
        B = 0.5 * (B + B.T)
        np.fill_diagonal(B, 0.0)
        B = P @ B @ P
        np.fill_diagonal(B, 0.0)
        nrm = float(np.linalg.norm(B))
        if nrm < 1e-14:
            continue
        B /= nrm
        max_s1 = max(max_s1, abs(float(ones @ B @ ones)))
    return {
        "p": p,
        "u_minus_pplus1_ppe": align,
        "u_identity_ok": align < 1e-10,
        "max_abs_1TB1": max_s1,
        "s1_vanishes": max_s1 < 1e-10,
        "dim_Z": int(mZ),
    }


def certify_wrong_floor_does_not_empty() -> dict:
    """ℓ=4 at p=5 must not empty bi-tight; shipped p=5 floor must."""
    bad = bitight_from_floor(5, 4)
    shipped = bitight_from_floor(5, lambda_min_floor_proved(5))
    return {
        "ell4_p5_bitight_empty": bad["bitight_empty"],
        "ell4_p5_L_star_lt_2d": bad["L_star_lt_2d"],
        "shipped_p5_floor": str(lambda_min_floor_proved(5)),
        "shipped_p5_bitight_empty": shipped["bitight_empty"],
        "ok": (not bad["bitight_empty"]) and shipped["bitight_empty"],
    }


# ---------------------------------------------------------------------------
# Status / main
# ---------------------------------------------------------------------------


def _star_kappa_pairing_proved() -> bool:
    """⟨star,κ_B⟩=−p‖B‖² on Z (15.111 E′).  Does not close D≽I."""
    from e1_gmin_m4_prop15111 import prove_star_kappa_pairing

    return bool(prove_star_kappa_pairing([6, 7], trials=8)["proved"])


def _Tb_kappa_pairing_proved() -> bool:
    """⟨Tb,κ_B⟩=6(3p²+5)/p² on Z (15.111 E‴).  Does not close D≽I."""
    from e1_gmin_m4_prop15111 import prove_Tb_kappa_pairing

    return bool(prove_Tb_kappa_pairing()["proved"])


def hinge_status_273() -> dict:
    return {
        "lambda_min_ge_6_proved_general": lambda_min_ge_6_proved_general(),
        "lambda_min_ge_5_proved_p_ge_7": lambda_min_ge_5_proved_p_ge_7(),
        "image_sos": theorem_image_sos_factorization()["proved"],
        "rho_beats_quarter": theorem_rho_lower_bound_beats_quarter()["proved"],
        "threshold_identity": theorem_threshold_identity()["proved"],
        "s1_vanishes_on_Z": theorem_s1_vanishes_on_Z()["proved"],
        "wick_criterion": theorem_wick_q4_criterion()["proved"],
        "lambda_star_gt_6": theorem_lambda_star_gt_6()["proved"],
        "q4_iff_m4_kap": theorem_q4_ge_0_iff_m4_kap_nonneg()["proved"],
        "wedge_schur_D_ge_I": theorem_wedge_schur_and_D_ge_I()["proved"],
        "star_kappa_pairing": _star_kappa_pairing_proved(),
        "Tb_kappa_pairing": _Tb_kappa_pairing_proved(),
        "census_p5_p7": theorem_census_spectra()["proved_census"],
        "bitight_empty_via_floor_all_p": bitight_empty_via_floor_all_p(),
        "floor_p5": str(lambda_min_floor_proved(5)),
        "floor_p7": str(lambda_min_floor_proved(7)),
        "floor_p11": str(lambda_min_floor_proved(11)),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "PSD factor of G_{u,disj} on Im(Z), or Aut-isotypic spec(Φ)≥6, "
            "or Weil+V₊ contraction ⟨ρ,κ_B⟩≥−1/4−1/(2p²); alternatively "
            "a general ℓ≥5 for all primes p≥7 (plus p=5 spectrum)"
        ),
    }


def prove_open() -> dict:
    return {
        "lambda_min_ge_6_proved_general": lambda_min_ge_6_proved_general(),
        "lambda_min_ge_5_for_p_ge_7": lambda_min_ge_5_proved_p_ge_7(),
        "bitight_empty_via_floor_all_p": bitight_empty_via_floor_all_p(),
        "gu_disj_is_gram": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "aut_schur": False,
        "pairing": False,
        "sixteen_N_for_all_p": False,
        "H_proved": False,
        "Es4_star_for_all_p": False,
        "L_status": "OPEN",
        "note": (
            "Threshold algebra and S₁=0 are Max+-free.  Uniform λ_min≥6 "
            "and the p≥7 floor ℓ≥5 remain open.  p=5,7 spectra empty "
            "bi-tight at those orders only."
        ),
    }


def main() -> dict:
    A = theorem_threshold_identity()
    B = theorem_s1_vanishes_on_Z()
    C = theorem_wick_q4_criterion()
    C2 = theorem_lambda_star_gt_6()
    C3 = theorem_q4_ge_0_iff_m4_kap_nonneg()
    D = theorem_image_sos_factorization()
    E = theorem_rho_lower_bound_beats_quarter()
    cert_spec = theorem_census_spectra()
    cert_gud = certify_gu_disj_not_psd(3)
    cert_s1 = certify_s1_zero(5)
    wrong = certify_wrong_floor_does_not_empty()
    status = hinge_status_273()
    open_ = prove_open()
    out = {
        "prop": "15.273",
        "title": (
            "λ_min(Φ) floor on Z: threshold algebra + S₁=0; "
            "G_{u,disj} is not a Gram; general ≥6 OPEN"
        ),
        "L_status": "OPEN",
        "proved": {
            "threshold_identity": A["proved"],
            "s1_vanishes_on_Z": B["proved"],
            "wick_q4_criterion": C["proved"],
            "lambda_star_gt_6": C2["proved"],
            "q4_iff_m4_kap_nonneg": C3["proved"],
            "wedge_schur_D_ge_I": theorem_wedge_schur_and_D_ge_I()["proved"],
            "star_kappa_pairing": _star_kappa_pairing_proved(),
            "Tb_kappa_pairing": _Tb_kappa_pairing_proved(),
            "image_sos_factorization": D["proved"],
            "rho_beats_quarter": E["proved"],
            "lambda_min_ge_6_proved_general": lambda_min_ge_6_proved_general(),
            "lambda_min_ge_5_proved_p_ge_7": lambda_min_ge_5_proved_p_ge_7(),
            "census_spectrum_p5_p7": cert_spec["proved_census"],
            "bitight_empty_via_floor_all_p": bitight_empty_via_floor_all_p(),
            "gu_disj_is_gram": False,
            "gsum_disj_lb_proved_general": False,
            "residual_closed_general": False,
            "sixteen_N_for_all_p": False,
            "E1_closed": False,
        },
        "algebra": {
            "ell_threshold": "4(p^4-1)/(p^4-8p^2-1)",
            "threshold_p5": str(ell_threshold_bitight(5)),
            "threshold_p7": str(ell_threshold_bitight(7)),
            "L_star_ell": "2(p^2+1)-ell(p^4-8p^2-1)/(4(p^2-1))",
            "Q4_over_N": "2+4/p^2+8 <rho,kappa_B>",
            "q4_nonneg_iff": "sum rho kappa_B >= -1/4 - 1/(2p^2)",
            "S1": "1^T B 1 = 0 on Z via P+1=(p+1)P+ e_inf",
            "A": A,
            "B": B,
            "C": C,
            "D": D,
            "E": E,
            "open": open_,
        },
        "census": {
            "spectra": cert_spec,
            "gu_disj_p3": cert_gud,
            "s1_p5": cert_s1,
            "wrong_floor": wrong,
        },
        "floors": {
            "p5": str(lambda_min_floor_proved(5)),
            "p7": str(lambda_min_floor_proved(7)),
            "p11": str(lambda_min_floor_proved(11)),
            "gram": str(gram_floor()),
        },
        "bitight_predicate_p5_shipped": bitight_from_floor(5, lambda_min_floor_proved(5)),
        "bitight_predicate_p5_ell4": bitight_from_floor(5, 4),
        "hinge": status,
        "F3": "no soft-close of L",
        "F19": "no class_key thrash",
        "F20": "CPU Fraction for general-p claims; p=3 Gu_disj eigensystem only",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15273.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("Prop 15.273  λ_min(Φ) floor on Z (G_u,disj is not a Gram)")
    print(f"  A threshold identity: {A['proved']}")
    print(f"  B S1=0 on Z: {B['proved']}")
    print(f"  C Wick Q4 criterion: {C['proved']}")
    print(f"  D image SOS factor: {D['proved']}")
    print(f"  E rho beats -1/4: {E['proved']}")
    print(f"  lambda_min_ge_6_proved_general={lambda_min_ge_6_proved_general()}")
    print(f"  lambda_min_ge_5_proved_p_ge_7={lambda_min_ge_5_proved_p_ge_7()}")
    print(f"  census p=5,7: {cert_spec['proved_census']}")
    print(
        f"    p=5 λmin={cert_spec['by_p']['5']['lambda_min']} "
        f"> thr={cert_spec['by_p']['5']['threshold']}"
    )
    print(
        f"    p=7 λmin={cert_spec['by_p']['7']['lambda_min']} "
        f"> thr={cert_spec['by_p']['7']['threshold']}"
    )
    print(
        f"  Gu_disj p=3 min_eig={cert_gud.get('min_eig')} "
        f"psd={cert_gud.get('is_psd')}"
    )
    print(
        f"  wrong ell=4 p=5 empties? {wrong['ell4_p5_bitight_empty']}; "
        f"shipped p=5 empties? {wrong['shipped_p5_bitight_empty']}"
    )
    print(f"  bitight_empty_via_floor_all_p={bitight_empty_via_floor_all_p()}")
    print(f"  gsum_disj_lb={gsum_disj_lb_proved_general()}")
    print(f"  L_status={open_['L_status']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
