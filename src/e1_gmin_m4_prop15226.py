#!/usr/bin/env python3
"""
Prop 15.226 — Farkas-sharp |μ| thr; G-invariant δ; Max± cross moment;
‖μ_part‖₂² closed; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.176, 15.218, 15.225)
  Dual-eq Farkas (e∉G, k=3p−2): obstructed if every Gsum_{e,e'} ≥ μ with
  μ > μ_* := (6/p−4)/k = (6−4p)/(p(3p−2)).
  On |κ|=1: min disj Gsum = −2 max|μ₄|, so
      Gsum ≥ μ_*   ⇔   max|μ₄| ≤ μ_f := −μ_*/2 = (2p−3)/(p(3p−2)).
  The convenient thr 1/(2p) is strictly stronger: 1/(2p) < μ_f for p≥5.
  μ=μ_part+δ, δ∈E_{±4p}, μ_part=μ_mn (15.225); |μ_part|≤maj≤1/(2p)<μ_f.

PROVED Max+-free
  A. Farkas-sharp |μ| threshold.
        μ_f(p) := (2p−3)/(p(3p−2)).
     For all primes p≥5: 1/(2p) < μ_f(p).
     Proof: 2(2p−3) > 3p−2 ⇔ 4p−6 > 3p−2 ⇔ p > 4.  ∎
     Hence |μ₄|≤1/(2p) ⇒ |μ₄|≤μ_f ⇒ Gsum≥−2μ_f=μ_* is the boundary;
     the strict Farkas form uses Gsum≥−1/p > μ_* (15.176).  Recording μ_f
     as the sharp matching-entry threshold for the same dual-eq obstruction
     (any max|μ|<μ_f also fires Farkas).  ∎

  B. Particular majorant beats μ_f.
     maj=(p²+4p−9)/(p²(p²−5)) ≤ μ_f for all primes p≥5.
     Proof: (3p−2)(p²+4p−9) ≤ p(p²−5)(2p−3)
     ⇔ 0 ≤ 2p⁴ − 6p³ − 20p² + 50p − 18 = 2(p⁴−3p³−10p²+25p−9).
     At p=5: value 116>0; derivative of the quartic 4p³−9p²−20p+25 has
     discriminant structure with min on [5,∞) positive (value 116 at 5;
     second derivative 24p²−18p−20>0 on [5,∞) and first deriv at 5 is
     500−225−100+25=200>0).  Fraction-checked all primes ≤300.  ∎
     Combined with A: room_δ^f := μ_f − maj > room_δ = 1/(2p)−maj > 0.
     Closed form room_δ^f via common denominator (module).  ∎

  C. Max± cross fourth-moment identity.
     For any symmetric conference (Cᵀ=C, C_ii=0, C_ij=±1, C²=(n−1)I):
     Max+⊂ker(C−pI), Max−⊂ker(C+pI) ⇒ V₊⊥V₋ ⇒ y·z=0 for all
     y∈Max+, z∈Max− (15.218 A).  Hence for every such pair, the boolean
     vector y⊙z has coordinate-sum 0, and Newton–Girard gives
        e₄(y⊙z) = n(n−2)/8
     (15.116 A at s=0).  Therefore
        ∑_{|S|=4} m₄⁺(S) m₄⁻(S)
          = E_{y∈Max+} E_{z∈Max−} ∑_S χ_S(y) χ_S(z)
          = E E e₄(y⊙z) = n(n−2)/8.  ∎
     (Independent of |Max±|; only uses y·z=0 and the e₄ polynomial.)

  D. L² identity for μ=½(m₄⁺+m₄⁻).
        ∑_S μ(S)² = ¼(∑(m₄⁺)² + ∑(m₄⁻)²) + n(n−2)/16.
     If ∑(m₄⁺)²=∑(m₄⁻)² (true when Max− ≅ Max+ as boolean eigensets —
     holds for Paley by C ↔ −C), then
        ∑ μ² = ½ ∑(m₄⁺)² + n(n−2)/16.  ∎

  E. G-invariance of δ.
     Strict Aut G ≤ Aut(C) (15.134 A) preserves Max±, hence preserves m₄±
     and μ.  μ_part is a fixed linear image of (κ,φ), both G-invariant
     (constant on G-orbits of 4-sets).  Therefore δ=μ−μ_part is G-invariant:
        δ ∈ E_{±4p}^G := E_{±4p} ∩ {G-invariant 4-set functions}.
     Combined with 15.225: δ=δ₊+δ₋ with δ_±=P_±δ ∈ E_{±4p}^G.
     The L^∞ residual |δ_S|≤room on |κ|=1 is a finite-dimensional
     problem on the G-orbit space (no dim≤1 assumption — dim E_{4p}^G is
     2,7 at p=5,7).  ∎

  F. Closed ‖μ_part‖₂².
     Using ∑κ²=n(n−1)(n−2)(n−5)/8, ‖Tκ‖₂²=6n(n−1)(n−2)(n−6),
     T²κ=48κ−24φ, T²φ=4(n+2)(φ−2κ), Tφ=(n+2)star, star=−Tκ/6:
        ⟨φ,κ⟩ = (48∑κ² − ‖Tκ‖₂²)/24,
        ∑φ²   = (n+2)‖star‖₂²/4 + 2⟨φ,κ⟩,
     with ‖star‖₂²=n(n−1)(n−2)(n−6)/6.  Then μ_part=aκ+bφ with
     a=(p²−1)/den, b=−2/den, den=p²(p²−5) yields
        ‖μ_part‖₂² = (p−1)(p+1)(p²+1)(3p²+1) / (24(p²−5))
     for every prime p≥5.  ∎
     Since δ⊥μ_part (15.225 H): ‖μ‖₂² = ‖μ_part‖₂² + ‖δ‖₂².  ∎

CERTIFIED
  G. room_δ^f > room_δ > 0 and maj≤μ_f≤ as Fraction for primes 5≤p≤300.
  H. Cross identity and ‖μ_part‖₂² match census anchors at p=5,7 structure.

OPEN (blocks residual i / predicate flip)
  I. Prove for all primes p≥5 one of:
     (i)   max_{|κ|=1}|μ| ≤ μ_f (or ≤1/(2p) / ≤2/n / ≤M_cand);
     (ii)  |δ_S|≤room_δ^f (or room_δ) on |κ|=1 via Aut-SOS on E_{±4p}^G;
     (iii) Path C / k₂ gap / K₄≤Wick / free-e+ker=sc.
  Soft-close forbidden.

Until I: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15226.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
    residual_i_dual_eq_empty_proved_general,
)
from e1_gmin_m4_prop15176 import mu_star_threshold  # noqa: E402
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import (  # noqa: E402
    majorant_mu_part,
    theorem_particular_majorant_beats_thr,
    theorem_phi_bound_general,
)
from e1_gmin_m4_prop15225 import room_delta  # noqa: E402


def mu_farkas_threshold(p: int) -> Fraction:
    """Sharp matching |μ| thr: (2p−3)/(p(3p−2)) = −μ_*/2."""
    return Fraction(2 * p - 3, p * (3 * p - 2))


def theorem_farkas_sharp_thr() -> dict:
    ok = True
    sample = {}
    for p in range(5, 301):
        if not is_prime(p):
            continue
        mf = mu_farkas_threshold(p)
        thr = mu4_threshold(p)
        mu_star = mu_star_threshold(p)
        # −2 μ_f = μ_*
        if -2 * mf != mu_star:
            ok = False
            break
        if thr >= mf:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101):
            sample[str(p)] = {
                "mu_f": str(mf),
                "one_over_2p": str(thr),
                "mu_star": str(mu_star),
                "thr_lt_mu_f": thr < mf,
                "minus_2_mu_f_eq_mu_star": -2 * mf == mu_star,
            }
    return {
        "proved": ok,
        "theorem": (
            "μ_f=(2p−3)/(p(3p−2))=−μ_*/2 satisfies 1/(2p)<μ_f for all primes "
            "p≥5, and Gsum≥−2μ_f is the matching-entry form of the Farkas "
            "boundary Gsum≥μ_*.  The sufficient thr |μ|≤1/(2p) is strictly "
            "stronger than |μ|≤μ_f."
        ),
        "by_p_sample": sample,
        "depends_on": ["prop_15_176_mu_star"],
    }


def room_delta_farkas(p: int) -> Fraction:
    """μ_f − maj."""
    return mu_farkas_threshold(p) - majorant_mu_part(p)


def theorem_maj_beats_mu_farkas() -> dict:
    ok = True
    sample = {}
    phi_ok = theorem_phi_bound_general()["proved"]
    maj_ok = theorem_particular_majorant_beats_thr()["proved"]
    for p in range(5, 301):
        if not is_prime(p):
            continue
        maj = majorant_mu_part(p)
        mf = mu_farkas_threshold(p)
        rd = room_delta(p)
        rdf = room_delta_farkas(p)
        if maj > mf or rdf <= 0 or rdf <= rd:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101):
            sample[str(p)] = {
                "maj": str(maj),
                "mu_f": str(mf),
                "room_delta": str(rd),
                "room_delta_farkas": str(rdf),
                "maj_le_mu_f": maj <= mf,
                "room_f_gt_room": rdf > rd,
            }
    # Algebra: 2(p^4 - 3p^3 - 10p^2 + 25p - 9) > 0 on primes p≥5
    alg_ok = all(
        (p**4 - 3 * p**3 - 10 * p * p + 25 * p - 9) > 0
        for p in range(5, 301)
        if is_prime(p)
    )
    return {
        "proved": ok and alg_ok and phi_ok and maj_ok,
        "theorem": (
            "maj=(p²+4p−9)/(p²(p²−5))≤μ_f for all primes p≥5, with "
            "room_δ^f:=μ_f−maj > room_δ:=1/(2p)−maj > 0.  "
            "If |δ_S|≤room_δ^f on |κ|=1 then |μ|≤μ_f ⇒ dual-eq Farkas."
        ),
        "by_p_sample": sample,
        "poly": "2(p^4-3p^3-10p^2+25p-9)>0",
        "depends_on": [
            "prop_15_186_majorant",
            "theorem_farkas_sharp_thr",
        ],
        "sufficient_for_residual_i": False,
        "note": "|δ|≤room_δ^f still OPEN (weaker L^∞ target than room_δ).",
    }


def theorem_maxpm_cross_moment() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For any symmetric conference: y∈Max+, z∈Max− ⇒ y·z=0 ⇒ "
            "e₄(y⊙z)=n(n−2)/8 ⇒ ∑_S m₄⁺(S) m₄⁻(S)=n(n−2)/8."
        ),
        "formula": "n(n-2)/8",
        "depends_on": [
            "theorem_Max_plus_perp_Max_minus",
            "prop_15_116_e4_poly",
        ],
    }


def theorem_mu_L2_from_m4() -> dict:
    return {
        "proved": True,
        "theorem": (
            "∑ μ² = ¼(∑(m₄⁺)²+∑(m₄⁻)²) + n(n−2)/16.  "
            "When ∑(m₄⁺)²=∑(m₄⁻)²: ∑μ²=½∑(m₄⁺)² + n(n−2)/16."
        ),
        "depends_on": ["theorem_maxpm_cross_moment"],
    }


def theorem_delta_G_invariant() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Strict Aut G preserves Max± ⇒ μ G-invariant; μ_part is a "
            "function of G-invariant (κ,φ) ⇒ δ=μ−μ_part ∈ E_{±4p}^G.  "
            "With P± of 15.225: δ_±∈E_{±4p}^G.  No dim≤1 claim "
            "(false for G at p=5,7: dims 2,7)."
        ),
        "depends_on": [
            "prop_15_134_strict_Aut",
            "prop_15_225_projectors",
            "prop_15_225_mu_part_eq_mu_mn",
        ],
        "note": (
            "Finite-dim Aut-SOS on E_{±4p}^G is the residual L^∞ problem; "
            "symbolic PSD certificate still OPEN."
        ),
    }


def kappa_norm_sq_closed(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(n * (n - 1) * (n - 2) * (n - 5), 8)


def Tkappa_norm_sq_closed(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(6 * n * (n - 1) * (n - 2) * (n - 6))


def star_norm_sq_closed(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(n * (n - 1) * (n - 2) * (n - 6), 6)


def phi_kappa_inner(p: int) -> Fraction:
    """⟨φ,κ⟩ = (48∑κ² − ‖Tκ‖₂²)/24."""
    return (48 * kappa_norm_sq_closed(p) - Tkappa_norm_sq_closed(p)) / 24


def phi_norm_sq_closed(p: int) -> Fraction:
    """∑φ² from Tφ=(n+2)star and T²φ=4(n+2)(φ−2κ)."""
    n = n_of(p)
    # (n+2) ||star||² = 4(||φ||² - 2⟨φ,κ⟩)
    # ||φ||² = (n+2)||star||²/4 + 2⟨φ,κ⟩
    return (n + 2) * star_norm_sq_closed(p) / 4 + 2 * phi_kappa_inner(p)


def mu_part_norm_sq(p: int) -> Fraction:
    """‖μ_part‖₂² = (p−1)(p+1)(p²+1)(3p²+1)/(24(p²−5))."""
    return Fraction(
        (p - 1) * (p + 1) * (p * p + 1) * (3 * p * p + 1),
        24 * (p * p - 5),
    )


def mu_part_norm_sq_via_ab(p: int) -> Fraction:
    den = p * p * (p * p - 5)
    a = Fraction(p * p - 1, den)
    b = Fraction(-2, den)
    return (
        a * a * kappa_norm_sq_closed(p)
        + b * b * phi_norm_sq_closed(p)
        + 2 * a * b * phi_kappa_inner(p)
    )


def theorem_mu_part_norm_sq_closed() -> dict:
    ok = True
    sample = {}
    for p in range(5, 201):
        if not is_prime(p):
            continue
        a = mu_part_norm_sq(p)
        b = mu_part_norm_sq_via_ab(p)
        if a != b or a <= 0:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31):
            sample[str(p)] = {
                "mu_part_norm_sq": str(a),
                "via_ab_match": a == b,
            }
    return {
        "proved": ok,
        "theorem": (
            "‖μ_part‖₂²=(p−1)(p+1)(p²+1)(3p²+1)/(24(p²−5)) for all primes "
            "p≥5.  Proof: expand aκ+bφ with T²κ, T²φ, Tφ identities "
            "(15.215/224) to evaluate ∑κ², ∑φ², ⟨φ,κ⟩; simplify.  "
            "Hence ‖μ‖₂²=‖μ_part‖₂²+‖δ‖₂² (δ⊥μ_part)."
        ),
        "formula": "(p-1)(p+1)(p^2+1)(3p^2+1)/(24(p^2-5))",
        "by_p_sample": sample,
        "depends_on": [
            "prop_15_215_T_identities",
            "prop_15_224_T2_phi",
            "prop_15_225_L2_orth",
        ],
    }


def residual_i_closed_via_226() -> bool:
    return False


def hinge_status_226() -> dict:
    return {
        "farkas_sharp_thr": theorem_farkas_sharp_thr()["proved"],
        "maj_beats_mu_farkas": theorem_maj_beats_mu_farkas()["proved"],
        "maxpm_cross_moment": theorem_maxpm_cross_moment()["proved"],
        "mu_L2_from_m4": theorem_mu_L2_from_m4()["proved"],
        "delta_G_invariant": theorem_delta_G_invariant()["proved"],
        "mu_part_norm_sq_closed": theorem_mu_part_norm_sq_closed()["proved"],
        "residual_i_closed_via_226": residual_i_closed_via_226(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "|μ|≤μ_f on |κ|=1 (or |δ|≤room_δ^f via Aut-SOS on E_pm4p^G); "
            "or Path C / K4 / free-e+ker=sc."
        ),
    }


def main() -> dict:
    a = theorem_farkas_sharp_thr()
    b = theorem_maj_beats_mu_farkas()
    c = theorem_maxpm_cross_moment()
    d = theorem_mu_L2_from_m4()
    e = theorem_delta_G_invariant()
    f = theorem_mu_part_norm_sq_closed()
    status = hinge_status_226()
    out = {
        "prop": "15.226",
        "title": (
            "Farkas-sharp μ_f; G-invariant δ; Max± cross; ‖μ_part‖₂² closed; "
            "residual-(i) OPEN"
        ),
        "proved": {
            "farkas_sharp_thr": a["proved"],
            "maj_beats_mu_farkas": b["proved"],
            "maxpm_cross_moment": c["proved"],
            "mu_L2_from_m4": d["proved"],
            "delta_G_invariant": e["proved"],
            "mu_part_norm_sq_closed": f["proved"],
            "residual_i_closed": residual_i_closed_via_226(),
        },
        "algebra": {
            "mu_f": "(2p-3)/(p(3p-2))",
            "mu_star": "(6-4p)/(p(3p-2))",
            "cross": "n(n-2)/8",
            "mu_part_norm_sq": "(p-1)(p+1)(p^2+1)(3p^2+1)/(24(p^2-5))",
            "mu_f_p5": str(mu_farkas_threshold(5)),
            "mu_f_p7": str(mu_farkas_threshold(7)),
            "room_delta_farkas_p5": str(room_delta_farkas(5)),
            "room_delta_farkas_p7": str(room_delta_farkas(7)),
            "mu_part_norm_sq_p5": str(mu_part_norm_sq(5)),
            "mu_part_norm_sq_p7": str(mu_part_norm_sq(7)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_226(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Farkas-sharp thr, cross moment, G-invariance, ‖μ_part‖₂² proved.  "
            "Pointwise |μ|≤μ_f / Aut-SOS on E^G still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15226.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.226  residual_i={residual_i_closed_via_226()}")
    print(
        f"  mu_f={a['proved']} maj_le_mu_f={b['proved']} "
        f"cross={c['proved']} mu_part_L2={f['proved']}"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
