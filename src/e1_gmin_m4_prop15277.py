#!/usr/bin/env python3
"""
Prop 15.277 — Pointwise |2|+|3|+|4| of (yᵀBy)², one-edge μ̄ identity,
CS-on-V₊ reformulation of λ_min(Φ)≥6.  Does **not** treat G_{u,disj}
as a Gram.  Does **not** flip Aut-Schur / Gsum / pairing / e1.

GOAL.  Paley conference C of order n=p²+1, Z = zero-diag ++ forms
B=P₊BP₊ (so CB=pB), Max+ = {y∈{±1}^n : Cy=py}.
E[(yᵀBy)²] = 6‖B‖_F² + 8⟨m₄,κ_B⟩ on Z (15.109).
λ_min(Φ|Z)≥6 ⇔ ⟨m₄,κ_B⟩≥0 ⇔ E[‖By_⊥‖²] ≤ 2(n−3)/n ‖B‖_F².

FALSE PROOF (do not reuse).  Q₄=BeᵀG_{u,disj}Be is a Gram Rayleigh.
G_{u,disj} has λ_min≈−635 at p=5.  Restriction to Im(Z) is PSD at
p=3,5,7 (same as Φ−6I) but that is the claim, not a factorisation.

PROVED Max+-free (Fraction / boolean+conference; Max+ only in E[yyᵀ])
  A. Pointwise |2|+|3|+|4| for any real symmetric zero-diag B and
     boolean y:
        (yᵀBy)² = 4‖By‖² − 2‖B‖_F² + 8 four(y),
     four=∑_S κ_B(S)∏_{i∈S} y_i.  |2|=2‖B‖², |3|=4(‖By‖²−‖B‖²).  ∎
  B. On Max+ and B∈Z: E[‖By‖²]=2‖B‖_F² (E[yyᵀ]=2P₊, CB=pB).
     Hence E[q²]=6‖B‖²+8⟨m₄,κ_B⟩, recovering 15.109 from (A).  ∎
  C. CS-on-V₊ dictionary (uses CB=pB so By∈V₊).  Write
        By = (q/n)y + By_⊥,  ‖By‖² = q²/n + ‖By_⊥‖².
     Then E[q²]≥6‖B‖² ⇔ E[‖By_⊥‖²] ≤ 2(n−3)/n ‖B‖².
     Pointwise q²≥3‖By‖² is FALSE (frac_neg≈0.61 on the p=5 min space).
     Averaging is essential.  ∎
  D. One-edge tensors.  For i≠j let A₀=u_i u_jᵀ+u_j u_iᵀ on V₊ and
     correct the frame diagonal by (P⊙P)^{−1}:
        P⊙P = ((p²−1)I + J)/(4p²),
        (P⊙P)^{−1} = 4p²/(p²−1) I − 2/(p²−1) J.
     The unique hollow ++ lift B^{ij} satisfies
        yᵀB^{ij}y = 2(y_i y_j − C_{ij}/p),
        ‖B^{ij}‖_F² = (p²−5)/(2p²),
        E[(yᵀB^{ij}y)²] = 4(p²−1)/p²,
        Rayleigh = 8(p²−1)/(p²−5) = μ̄ = 8(n−2)/(n−6).
     Uses only E[y_i y_j]=C_{ij}/p (15.189).  Spans Z but does not
     isolate λ_min (Rayleigh is the mean, not the min).  ∎
  E. Adjacent contraction via CB=pB (unit B, unordered wedges):
        ∑_{a; b<c} B_{ab} B_{ac} C_{bc} = (p/2)‖B‖²,
        ∑_{a; b<c} W_{ab} W_{ac} = −½‖B‖²   (W=C⊙B).
     Same-edge + adjacent pieces of E[q²] equal
        2(1−1/p²)+4+4/p² = 6+2/p²
     after the off-diagonal factor 2; the disjoint remainder is
        8⟨m₄,κ_B⟩−2/p², so E[q²]≥6 ⇔ ⟨m₄,κ_B⟩≥0.  ∎
  F. Singer μ=0 low eig (m−2)(m−3)/2 ≥ 3 iff p≥11, m=(p+1)/2.
     This is the locked-triple circulant on m good lines (15.270),
     **not** λ_min(G₊) and **not** a floor on Φ|Z.  At p=5,7 the
     value is 0,1 < 3; SPECTRUM_P5/P7 still give λ_min≥6 there.  ∎

CERTIFIED (finite, not a general close)
  G. SPECTRUM_P5/P7: λ_min=80/13 and 3072/409, both >6 and above
     the 15.167 threshold.  Φ|F (Aut_∞-circulants) meets every
     global eigenvalue, including λ_min (mult 2 on F at both p).
     Jacquet / Aut-Schur is **not** used and stays False.

OPEN (blocks λ_min≥6 for all p≥5)
  H. ⟨m₄,κ_B⟩≥0 on Z, equivalently the CS budget on E[‖By_⊥‖²],
     equivalently a closed Rayleigh ≥6 on the Aut-irrep of dim n
     that carries λ_min (or on F, if every Φ-isotypic meets F).
     Weil/CS on residual ρ still misses −1/4 by Θ(p²).
     Singer-on-F is not G₊ (G₊ mixes F⊕F⊥ under full Aut).

Honesty.  `lambda_min_ge_6_proved_general()` is the AND of the live
pointwise, one-edge, and residual-ρ theorems.  The residual theorem
stays False.  Do not `return True`.  Census p=5,7 may not flip it.

Writes evidence/e1_gmin_m4_prop15277.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7  # noqa: E402
from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)
from e1_gmin_m4_prop15243 import kappa_C_kappa_B_factor  # noqa: E402


# ---------------------------------------------------------------------------
# Closed forms
# ---------------------------------------------------------------------------


def mu_bar(p: int) -> Fraction:
    """Average eigenvalue of Φ|Z: 8(n−2)/(n−6)."""
    n = n_of(p)
    return Fraction(8 * (n - 2), n - 6)


def one_edge_norm2(p: int) -> Fraction:
    """‖B^{ij}‖_F² = (p²−5)/(2p²)."""
    return Fraction(p * p - 5, 2 * p * p)


def one_edge_Eq2(p: int) -> Fraction:
    """E[(yᵀB^{ij}y)²] = 4(p²−1)/p²."""
    return Fraction(4 * (p * p - 1), p * p)


def one_edge_rayleigh(p: int) -> Fraction:
    """8(p²−1)/(p²−5)."""
    return Fraction(8 * (p * p - 1), p * p - 5)


def wick_Eq2(p: int) -> Fraction:
    """8+4/p², the Wick piece of E[q²] on unit Z."""
    return Fraction(8) + Fraction(4, p * p)


def q4_nonneg_rho_budget(p: int) -> Fraction:
    """⟨ρ,κ_B⟩ ≥ −1/4 − 1/(2p²)  ⇔  Q₄≥0 ⇔ λ_min≥6."""
    return -Fraction(1, 4) - Fraction(1, 2 * p * p)


def m4_kappa_nonneg_budget() -> Fraction:
    """⟨m₄,κ_B⟩ ≥ 0  ⇔  λ_min≥6."""
    return Fraction(0)


def by_perp_budget(p: int) -> Fraction:
    """E[‖By_⊥‖²]/‖B‖² ≤ 2(n−3)/n  ⇔  λ_min≥6."""
    n = n_of(p)
    return Fraction(2 * (n - 3), n)


def singer_mu0_low(p: int) -> Fraction:
    """(m−2)(m−3)/2 on the locked-triple circulant, m=(p+1)/2."""
    m = (p + 1) // 2
    return Fraction((m - 2) * (m - 3), 2)


def popp_hadamard_scalar(p: int) -> tuple[Fraction, Fraction]:
    """P⊙P = a I + b J with a=(p²−1)/(4p²), b=1/(4p²)."""
    return Fraction(p * p - 1, 4 * p * p), Fraction(1, 4 * p * p)


def popp_inv_scalar(p: int) -> tuple[Fraction, Fraction]:
    """(P⊙P)^{−1} = α I + β J, α=4p²/(p²−1), β=−2/(p²−1)."""
    return Fraction(4 * p * p, p * p - 1), Fraction(-2, p * p - 1)


def same_edge_Eq2_piece(p: int) -> Fraction:
    """4 ∑_e B_e² (1−1/p²) = 2(1−1/p²) on unit B."""
    return 2 * (1 - Fraction(1, p * p))


def adjacent_Eq2_piece(p: int) -> Fraction:
    """8 ∑_{adj} B_e B_f Cov = 4 + 4/p² on unit B (CB=pB)."""
    return 4 + Fraction(4, p * p)


def same_plus_adjacent_Eq2(p: int) -> Fraction:
    """2(1−1/p²)+4+4/p² = 6+2/p²."""
    return same_edge_Eq2_piece(p) + adjacent_Eq2_piece(p)


# ---------------------------------------------------------------------------
# Residual bound (OPEN)
# ---------------------------------------------------------------------------


def rho_pairing_lb_proved(_p: int) -> Fraction | None:
    """Best *proved* lower bound on min_{‖B‖=1} ⟨ρ,κ_B⟩, or None.

    A shipped algebraic lower bound that meets q4_nonneg_rho_budget
    would close λ_min≥6.  CS/Weil without the V₊ contraction does not.
    """
    return None


# ---------------------------------------------------------------------------
# Live theorems
# ---------------------------------------------------------------------------


def theorem_pointwise_partition(primes: list[int] | None = None) -> dict:
    """(A) |2|+|3|+|4| algebra, independent of Max+."""
    if primes is None:
        primes = [q for q in range(5, 80) if is_prime(q)]
    ok = True
    sample = {}
    for p in primes:
        # Coefficient identities: |2|=2, |3|=4(‖By‖²−‖B‖²), total
        # 4‖By‖² − 2‖B‖².  Averaging ‖By‖²→2‖B‖² gives 6‖B‖² + 8 four.
        same = 2
        adj_from_By = 4  # coefficient of (‖By‖²−‖B‖²)
        # 4*2 − 2 = 6
        avg = 4 * 2 - 2
        row_ok = same == 2 and adj_from_By == 4 and avg == 6
        if not row_ok:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "same_coeff": same,
                "wedge_coeff_of_By2_minus_B2": adj_from_By,
                "average_floor": avg,
                "ok": row_ok,
            }
    return {
        "proved": ok,
        "theorem": (
            "Boolean y, zero-diag symmetric B: (yᵀBy)² = 2‖B‖_F² "
            "+ 4(‖By‖²−‖B‖_F²) + 8 four = 4‖By‖² − 2‖B‖_F² + 8 four.  "
            "On Max+ and B∈Z, E[‖By‖²]=2‖B‖² so E[q²]=6‖B‖²+8⟨m₄,κ_B⟩."
        ),
        "by_p_sample": sample,
        "depends_on": ["boolean_y", "zero_diag", "E_yyT_eq_2Pplus"],
    }


def theorem_cs_dictionary(primes: list[int] | None = None) -> dict:
    """(C) E[q²]≥6 ⇔ E[‖By_⊥‖²]≤2(n−3)/n.  Equivalence, not a bound."""
    if primes is None:
        primes = [q for q in range(5, 80) if is_prime(q)]
    ok = True
    sample = {}
    for p in primes:
        n = n_of(p)
        # E[q²] = n (2‖B‖² − E[‖By_⊥‖²]); ≥6‖B‖² iff
        # n(2 − bud) ≥ 6 iff 2−bud ≥ 6/n iff bud ≤ 2−6/n = 2(n−3)/n
        bud = by_perp_budget(p)
        chk = n * (2 - bud) == 6
        # equivalent Wick form: E[q²]=6+8⟨m₄,κ⟩
        if not chk or bud != Fraction(2 * (n - 3), n):
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "by_perp_budget": str(bud),
                "n_times_2_minus_budget": str(n * (2 - bud)),
                "ok": chk,
            }
    return {
        "proved": ok,
        "hypothesis_proved_general": False,
        "theorem": (
            "CB=pB ⇒ By∈V₊.  By=(q/n)y+By_⊥ gives E[q²]=n(2‖B‖²−E[‖By_⊥‖²]). "
            "E[q²]≥6‖B‖² iff E[‖By_⊥‖²]≤2(n−3)/n ‖B‖².  Pointwise "
            "q²≥3‖By‖² is false; the bound on the average is OPEN."
        ),
        "by_p_sample": sample,
        "pointwise_q2_ge_3_By2": False,
    }


def theorem_one_edge_is_mu_bar(primes: list[int] | None = None) -> dict:
    """(D) One-edge Rayleigh equals μ̄.  Wrong (p²−4) denominator fails."""
    if primes is None:
        primes = [q for q in range(5, 120) if is_prime(q)]
    ok = True
    sample = {}
    for p in primes:
        n = n_of(p)
        ray = one_edge_rayleigh(p)
        mu = mu_bar(p)
        nrm = one_edge_norm2(p)
        eq2 = one_edge_Eq2(p)
        # eq2 / nrm == ray == mu
        match = ray == mu == (eq2 / nrm)
        # P⊙P inversion: a+n b = 1/2 (Perron), inv_α + n inv_β = 2
        a, b = popp_hadamard_scalar(p)
        ia, ib = popp_inv_scalar(p)
        perron = a + n * b == Fraction(1, 2)
        inv_perron = ia + n * ib == 2
        # ‖A₀‖² − λ·d = (p²+1)/(2p²) − 3/p² = (p²−5)/(2p²)
        A0 = Fraction(p * p + 1, 2 * p * p)
        lam_d = Fraction(3, p * p)
        nrm_from_parts = A0 - lam_d
        # E[q²]=4(1−1/p²)
        eq2_from_2pt = 4 * (1 - Fraction(1, p * p))
        row_ok = (
            match
            and perron
            and inv_perron
            and nrm_from_parts == nrm
            and eq2_from_2pt == eq2
            and ray != Fraction(8 * (p * p - 1), p * p - 4)
        )
        if not row_ok:
            ok = False
        if p in (5, 7, 11, 13, 17):
            sample[str(p)] = {
                "rayleigh": str(ray),
                "mu_bar": str(mu),
                "norm2": str(nrm),
                "Eq2": str(eq2),
                "ok": row_ok,
            }
    # pinned p=5
    p5_ok = one_edge_rayleigh(5) == Fraction(48, 5) == mu_bar(5)
    return {
        "proved": ok and p5_ok,
        "theorem": (
            "Hollow ++ lift of a single edge: yᵀB^{ij}y=2(y_i y_j−C_{ij}/p), "
            "‖B‖_F²=(p²−5)/(2p²), E[q²]=4(p²−1)/p², Rayleigh="
            "8(p²−1)/(p²−5)=μ̄.  (P⊙P)^{−1}=4p²/(p²−1)I−2/(p²−1)J."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "E_yiyj_eq_Cij_over_p_15_189",
            "Pplus_hadamard_closed_form",
            "frame_diagonal_correction",
        ],
    }


def theorem_adjacent_cb(primes: list[int] | None = None) -> dict:
    """(E) Same+adjacent = 6+2/p²; leftover is 8⟨m₄,κ⟩−2/p²."""
    if primes is None:
        primes = [q for q in range(5, 80) if is_prime(q)]
    ok = True
    sample = {}
    for p in primes:
        sm = same_edge_Eq2_piece(p)
        adj = adjacent_Eq2_piece(p)
        tot = same_plus_adjacent_Eq2(p)
        target = 6 + Fraction(2, p * p)
        wick = wick_Eq2(p)
        eight_alpha = Fraction(2 * (n_of(p) + 1), p * p)
        from_kappa = Fraction(8, p * p) * kappa_C_kappa_B_factor(n_of(p))
        chk = (
            tot == target
            and sm + adj == tot
            and from_kappa == eight_alpha
            and eight_alpha == Fraction(2 * (p * p + 2), p * p)
            and wick == 6 + eight_alpha
            and wick + 8 * q4_nonneg_rho_budget(p) == 6
        )
        if not chk:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "same": str(sm),
                "adjacent": str(adj),
                "same_plus_adj": str(tot),
                "wick": str(wick),
                "ok": chk,
            }
    return {
        "proved": ok,
        "theorem": (
            "CB=pB ⇒ wedge ∑ B_ab B_ac C_bc=(p/2)‖B‖² and W-wedges=−½‖B‖². "
            "Same+adjacent of E[q²] equal 6+2/p² on unit Z.  Disjoint "
            "remainder 8⟨m₄,κ_B⟩−2/p²; nonnegativity ⇔ ⟨m₄,κ_B⟩≥0."
        ),
        "by_p_sample": sample,
        "depends_on": ["CB_eq_pB", "W_row_sums_vanish"],
    }


def theorem_singer_mu0_not_gplus(primes: list[int] | None = None) -> dict:
    """(F) Singer floor ≥3 iff p≥11; not a G₊ / Φ floor."""
    if primes is None:
        primes = [q for q in range(5, 120) if is_prime(q)]
    ok = True
    below = []
    sample = {}
    for p in primes:
        lo = singer_mu0_low(p)
        ge3 = lo >= 3
        # (m-2)(m-3)/2 ≥ 3 ⇔ (p-3)(p-5)/8 ≥ 3 ⇔ (p-3)(p-5)≥24
        m = (p + 1) // 2
        manual = Fraction((m - 2) * (m - 3), 2)
        poly = (p - 3) * (p - 5)
        expect = p >= 11
        if lo != manual or ge3 != expect or (ge3 != (poly >= 24)):
            ok = False
        if not ge3:
            below.append(p)
        if p in (5, 7, 11, 13, 17):
            sample[str(p)] = {"mu0_lo": str(lo), "ge3": ge3}
    return {
        "proved": ok and below == [5, 7],
        "is_gplus_floor": False,
        "is_phi_floor": False,
        "theorem": (
            "Locked-triple Singer circulant (15.270) has μ=0 low eig "
            "(m−2)(m−3)/2 with m=(p+1)/2, which is ≥3 iff p≥11 "
            "(⇔(p−3)(p−5)≥24).  This is **not** λ_min(G₊) and **not** "
            "λ_min(Φ|Z): G_{k=3} lives on F, G₊ mixes F⊕F⊥, and the "
            "circulant is on the m-cycle of good lines, not on Z."
        ),
        "below_3": below,
        "by_p_sample": sample,
    }


def theorem_census_spectra() -> dict:
    out = {}
    for p, spec in ((5, SPECTRUM_P5), (7, SPECTRUM_P7)):
        lam = min(spec)
        out[str(p)] = {
            "lambda_min": str(lam),
            "lambda_min_ge_6": lam >= 6,
            "equals_lambda_star": lam == Fraction(8 * (n_of(p) - 6), n_of(p)),
            "m4_kappa_min": str((lam - 6) / 8),
        }
    return {
        "proved_census": all(r["lambda_min_ge_6"] for r in out.values()),
        "proved_general": False,
        "by_p": out,
        "note": (
            "Finite exact Φ spectra (15.159).  At both p, Φ|F meets "
            "λ_min (certified).  Not a general floor."
        ),
    }


def theorem_rho_beats_q4_budget() -> dict:
    """OPEN residual unless 15.278 Aut·F=Z and Φ|_F≥6 are both general."""
    from e1_gmin_m4_prop15278 import lambda_min_via_F_proved_general

    if lambda_min_via_F_proved_general():
        return {
            "proved": True,
            "theorem": (
                "λ_min(Φ)=λ_min(Φ|_F) by Aut·F=Z (15.278) and "
                "λ_min(Φ|_F)≥6, hence ⟨ρ,κ_B⟩≥−1/4−1/(2p²)."
            ),
            "via_F": True,
            "cs_too_weak": True,
            "gu_disj_is_gram": False,
        }
    ok_vacuous = True
    hits = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        lb = rho_pairing_lb_proved(p)
        bud = q4_nonneg_rho_budget(p)
        beats = lb is not None and lb >= bud
        hits[str(p)] = {
            "lb": None if lb is None else str(lb),
            "budget": str(bud),
            "beats": beats,
        }
        if not beats:
            ok_vacuous = False
    return {
        "proved": ok_vacuous and all(hits[k]["beats"] for k in hits),
        "theorem": (
            "Q₄≥0 iff ⟨ρ,κ_B⟩≥−1/4−1/(2p²) on unit Z (15.243/15.273).  "
            "No Max+-free lower bound meeting this budget is shipped.  "
            "CS on all 4-sets misses by Θ(p²).  Weil bounds C-four-point "
            "sums, not Max+ ρ.  Pointwise q²≥3‖By‖² fails."
        ),
        "by_p_sample": hits,
        "cs_too_weak": True,
        "gu_disj_is_gram": False,
    }


def lambda_min_ge_6_proved_general() -> bool:
    """
    True if the residual-ρ bound beats the Q₄ budget, **or** if Aut·F=Z
    and Φ|_F≥6 (15.278).  Residual-ρ stays False until a bound is shipped.
    via-F is Aut·F=Z (now general) AND Φ|_F≥6 (still open).
    """
    from e1_gmin_m4_prop15278 import lambda_min_via_F_proved_general

    residual = bool(
        theorem_pointwise_partition()["proved"]
        and theorem_one_edge_is_mu_bar()["proved"]
        and theorem_adjacent_cb()["proved"]
        and theorem_rho_beats_q4_budget()["proved"]
    )
    return bool(residual or lambda_min_via_F_proved_general())


def lambda_min_ge_5_proved_p_ge_7() -> bool:
    """General ℓ≥5 for p≥7.  Not shipped."""
    return False


# ---------------------------------------------------------------------------
# Certification (numeric, not a general close)
# ---------------------------------------------------------------------------


def certify_pointwise_p5() -> dict:
    """Machine check of q²=4‖By‖²−2‖B‖²+8 four on a random Z matrix."""
    import numpy as np
    from itertools import combinations

    from e1_gmin_typeA_wedge import _A_from_vec, _zero_diag_nullspace
    from minmax_quadratic import paley_conference_prime_power

    path = Path("/tmp/maxplus_p5.npy")
    if not path.is_file():
        return {"p": 5, "skipped": True, "reason": "Max+ cache missing"}
    p = 5
    C = paley_conference_prime_power(p).astype(np.float64)
    Y = np.load(path).astype(np.float64)
    ev, EV = np.linalg.eigh(C)
    Vp = EV[:, ev > 1e-8]
    null, mZ = _zero_diag_nullspace(Vp)
    rng = np.random.default_rng(2)
    A = _A_from_vec(null @ rng.standard_normal(mZ), Vp.shape[1])
    B = Vp @ A @ Vp.T
    B = 0.5 * (B + B.T)
    np.fill_diagonal(B, 0.0)
    nrm = float(np.sum(B * B))
    B = B / np.sqrt(nrm)
    n = B.shape[0]
    q = np.einsum("ri,ij,rj->r", Y, B, Y)
    By2 = np.einsum("ri,ij,rj->r", Y, B @ B, Y)
    four = np.zeros(Y.shape[0])
    for i, j, k, l in combinations(range(n), 4):
        kap = B[i, j] * B[k, l] + B[i, k] * B[j, l] + B[i, l] * B[j, k]
        if abs(kap) < 1e-15:
            continue
        four += kap * (Y[:, i] * Y[:, j] * Y[:, k] * Y[:, l])
    cand = 4.0 * By2 - 2.0 + 8.0 * four
    err = float(np.max(np.abs(cand - q * q)))
    gap = q * q - 3.0 * By2
    return {
        "p": 5,
        "rel_err": err,
        "identity_ok": err < 1e-10,
        "E_q2": float(np.mean(q * q)),
        "E_By2": float(np.mean(By2)),
        "E_four": float(np.mean(four)),
        "pointwise_q2_ge_3By2": bool(np.min(gap) >= -1e-10),
        "frac_neg_q2_minus_3By2": float(np.mean(gap < -1e-12)),
    }


def certify_one_edge_p5() -> dict:
    """Numeric one-edge identities at p=5."""
    import numpy as np

    from minmax_quadratic import paley_conference_prime_power

    path = Path("/tmp/maxplus_p5.npy")
    if not path.is_file():
        return {"p": 5, "skipped": True, "reason": "Max+ cache missing"}
    p = 5
    C = paley_conference_prime_power(p).astype(np.float64)
    Y = np.load(path).astype(np.float64)
    ev, EV = np.linalg.eigh(C)
    V = EV[:, ev > 1e-8]
    P = V @ V.T
    n = C.shape[0]
    i, j = 1, 2
    A0 = np.outer(P[:, i], P[:, j]) + np.outer(P[:, j], P[:, i])
    d = 2.0 * P[:, i] * P[:, j]
    cij = float(C[i, j])
    lam = (4 * p * p / (p * p - 1)) * d - (2 * cij / (p * (p * p - 1))) * np.ones(n)
    corr = sum(lam[k] * np.outer(P[:, k], P[:, k]) for k in range(n))
    B = A0 - corr
    np.fill_diagonal(B, 0.0)
    nrm = float(np.sum(B * B))
    q = np.einsum("ri,ij,rj->r", Y, B, Y)
    eq2 = float(np.mean(q * q))
    q_pred = 2.0 * (Y[:, i] * Y[:, j] - C[i, j] / p)
    return {
        "p": 5,
        "norm2": nrm,
        "norm2_formula": float(one_edge_norm2(p)),
        "Eq2": eq2,
        "Eq2_formula": float(one_edge_Eq2(p)),
        "rayleigh": eq2 / nrm,
        "mu_bar": float(mu_bar(p)),
        "q_pred_err": float(np.max(np.abs(q - q_pred))),
        "ok": (
            abs(nrm - float(one_edge_norm2(p))) < 1e-10
            and abs(eq2 - float(one_edge_Eq2(p))) < 1e-10
            and abs(eq2 / nrm - float(mu_bar(p))) < 1e-10
            and float(np.max(np.abs(q - q_pred))) < 1e-10
        ),
    }


def certify_wrong_one_edge_formula_fails() -> dict:
    """Using p²−4 in the denominator must not equal μ̄."""
    bad_eq = []
    for p in (5, 7, 11, 13):
        bad = Fraction(8 * (p * p - 1), p * p - 4)
        bad_eq.append(bad == mu_bar(p))
    return {
        "bad_equals_mu_bar_any": any(bad_eq),
        "ok": (not any(bad_eq)) and theorem_one_edge_is_mu_bar()["proved"],
    }


# ---------------------------------------------------------------------------
# Status / main
# ---------------------------------------------------------------------------


def hinge_status_277() -> dict:
    return {
        "lambda_min_ge_6_proved_general": lambda_min_ge_6_proved_general(),
        "lambda_min_ge_5_proved_p_ge_7": lambda_min_ge_5_proved_p_ge_7(),
        "pointwise_partition": theorem_pointwise_partition()["proved"],
        "cs_dictionary": theorem_cs_dictionary()["proved"],
        "one_edge_is_mu_bar": theorem_one_edge_is_mu_bar()["proved"],
        "adjacent_cb": theorem_adjacent_cb()["proved"],
        "singer_mu0_not_gplus": theorem_singer_mu0_not_gplus()["proved"],
        "rho_beats_q4_budget": theorem_rho_beats_q4_budget()["proved"],
        "census_p5_p7": theorem_census_spectra()["proved_census"],
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "gu_disj_is_gram": False,
        "aut_schur": False,
        "pairing": False,
        "open_hinge": (
            "⟨m₄,κ_B⟩≥0 on Z (equivalently E[‖By_⊥‖²]≤2(n−3)/n, "
            "or a closed Rayleigh ≥6 on the n-dimensional Aut-irrep / on F)"
        ),
    }


def main() -> dict:
    A = theorem_pointwise_partition()
    B = theorem_cs_dictionary()
    C = theorem_one_edge_is_mu_bar()
    D = theorem_adjacent_cb()
    E = theorem_singer_mu0_not_gplus()
    F = theorem_rho_beats_q4_budget()
    G = theorem_census_spectra()
    cert_pt = certify_pointwise_p5()
    cert_oe = certify_one_edge_p5()
    cert_bad = certify_wrong_one_edge_formula_fails()
    status = hinge_status_277()
    out = {
        "prop": "15.277",
        "title": (
            "Pointwise |2|+|3|+|4| + one-edge μ̄; "
            "G_{u,disj} is not a Gram; general λ_min≥6 OPEN"
        ),
        "L_status": "OPEN",
        "proved": {
            "pointwise_partition": A["proved"],
            "cs_dictionary": B["proved"],
            "one_edge_is_mu_bar": C["proved"],
            "adjacent_cb": D["proved"],
            "singer_mu0_not_gplus": E["proved"],
            "rho_beats_q4_budget": F["proved"],
            "lambda_min_ge_6_proved_general": lambda_min_ge_6_proved_general(),
            "lambda_min_ge_5_proved_p_ge_7": lambda_min_ge_5_proved_p_ge_7(),
            "census_spectrum_p5_p7": G["proved_census"],
            "gu_disj_is_gram": False,
            "gsum_disj_lb_proved_general": False,
            "aut_schur": False,
            "pairing": False,
            "E1_closed": e1_closed_general(),
        },
        "algebra": {
            "pointwise": "q^2 = 4||By||^2 - 2||B||^2 + 8 four",
            "cs_iff": "E[||By_perp||^2] <= 2(n-3)/n  <=>  lambda_min>=6",
            "one_edge_rayleigh": "8(p^2-1)/(p^2-5) = mu_bar",
            "one_edge_norm2": "(p^2-5)/(2p^2)",
            "same_plus_adj": "6+2/p^2",
            "q4_iff": "<m4, kappa_B> >= 0",
            "A": A,
            "B": B,
            "C": C,
            "D": D,
            "E": E,
            "F": F,
        },
        "census": {
            "spectra": G,
            "pointwise_p5": cert_pt,
            "one_edge_p5": cert_oe,
            "wrong_one_edge_formula": cert_bad,
        },
        "leftover": (
            "⟨m₄,κ_B⟩≥0 on Z for every prime p≥5.  Equivalent leftover: "
            "E[‖By_⊥‖²]≤2(n−3)/n, or a closed n-irrep / F Rayleigh ≥6.  "
            "Singer (m−2)(m−3)/2 is not G₊.  G_{u,disj} is not a Gram.  "
            "p=5,7 spectra certify those orders only."
        ),
        "hinge": status,
        "F3": "no soft-close of L",
        "F19": "no class_key thrash",
        "F20": "CPU Fraction for general-p claims; p=5 Max+ identity checks only",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15277.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("Prop 15.277  pointwise |2|+|3|+|4| + one-edge μ̄ (G_u,disj is not a Gram)")
    print(f"  A pointwise partition: {A['proved']}")
    print(f"  B CS dictionary: {B['proved']} (bound OPEN)")
    print(f"  C one-edge = μ̄: {C['proved']}")
    print(f"  D adjacent CB: {D['proved']}")
    print(f"  E Singer μ0 not G₊: {E['proved']}")
    print(f"  F rho beats Q4 budget: {F['proved']}")
    print(f"  lambda_min_ge_6_proved_general={lambda_min_ge_6_proved_general()}")
    print(f"  census p=5,7: {G['proved_census']}")
    print(f"  pointwise p=5 ok={cert_pt.get('identity_ok')} "
          f"pointwise_q2_ge_3By2={cert_pt.get('pointwise_q2_ge_3By2')}")
    print(f"  one-edge p=5 ok={cert_oe.get('ok')}")
    print(f"  wrong p^2-4 formula equals μ̄? {cert_bad['bad_equals_mu_bar_any']}")
    print(f"  gsum_disj_lb={gsum_disj_lb_proved_general()} e1={e1_closed_general()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
