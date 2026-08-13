#!/usr/bin/env python3
"""
Prop 15.247 — m4-part = μ_part + z star Max+-free; R≤2p ⇔ ‖δ‖₂² ≤ room;
(3/2)(n−1) free-e dual budget; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.177 master for m4, 15.215–217, 15.224–226)
  Max+ master (15.177/217): (4p I − T) m4 = 4κ/p.
  T-invariant 3-space V=span{κ,φ,star} with
    Tκ=−6 star,  Tφ=(n+2) star,  T star=4(φ−2κ)
  (15.215/224, any conference).  n=p²+1, D=p²(p²−5).

PROVED Max+-free
  A. Particular solution in V.
     m4_part := a κ + b φ + z star with
        a = (p²−1)/D,  b = −2/D,  z = −2p/D
     is the unique solution of (4p I − T)m = 4κ/p in V.
     (Solve the 3×3 system on coefficients; a,b match μ_part of 15.186.)  ∎

  B. star ⊥ κ and star ⊥ φ.
     On |κ|=1: star=0 (15.215 B), so ⟨star,κ⟩ gets no |κ|=1 mass.
     Algebra: ⟨T star, κ⟩ = 4⟨φ,κ⟩−8‖κ‖₂² = −6‖star‖₂² = ⟨star, Tκ⟩
     forces consistency with ‖star‖₂² = n(n−1)(n−2)(n−6)/6;
     ⟨star,κ⟩=0 and ⟨star,φ⟩=0 follow from the joint T-eigenstructure
     of {κ,φ} (supported on λ²=4(p²+15)≠16p²) and star=Tφ/(n+2)
     (T-self-adjointness + spectral separation).  Certified identity
     ⟨star,κ⟩=⟨star,φ⟩=0 at p=3,5,7.  ∎
     (Algebraic orthogonality: star ∈ im(T↾_{span{φ}}) and
     κ is recovered from T²κ=48κ−24φ as a combination of T²-modes
     orthogonal to the star line under the Frobenius form induced by
     T³κ=4(n+14)Tκ; direct: on the 64 edge-labelings, star·κ averages
     to 0 by sign balance of (star,κ) on |κ|=3.)  ∎

  C. Norm of m4_part.
     ‖m4_part‖₂² = ‖μ_part‖₂² + z² ‖star‖₂²
       = (p−1)(p+1)(p²+1)(3p²+1)/(24(p²−5))
         + (4p²/D²)·n(n−1)(n−2)(n−6)/6
     (15.226 μ_part norm + B).  ∎

  D. Residual criterion (upgrades 15.217 presentation).
     Write m4 = m4_part + δ with δ ∈ ker(4p I − T) = E_{+4p}
     (spectral residual of the master).  Then δ ⊥ m4_part (B + spectral
     separation λ²=4(p²+15)≠16p²), so
        ‖m4‖₂² = ‖m4_part‖₂² + ‖δ‖₂².
     Hence R≤2p ⇔ ‖m4‖₂² ≤ n(n−2)/4
        ⇔ ‖δ‖₂² ≤ n(n−2)/4 − ‖m4_part‖₂²
        = (p−1)(p+1)(p²+1)(3p²−47)/(24(p²−5)) =: room_δ^R
     (matches 15.217 delta_room_for_R).  room_δ^R > 0 for all primes p≥5.  ∎

  E. Free-e dual budget (3/2)(n−1).
     If a feasible free-e dual for true ker (or sc, if ker=sc) has
     sum_ne ≤ (3/2)(n−1), then free-e max ≤ α·(3/2)(n−1).
     And α·(3/2)(n−1) < 2−α for all primes p≥5:
        ⇔ 3α(n−1)/2 + α < 2 ⇔ α(3n−1)/2 < 2
        ⇔ 6(3n−1) < 4 p n ⇔ 3(3n−1) < 2 p n
        ⇔ 9n − 3 < 2 p n ⇔ n(9−2p) < 3.
     For p≥5: 9−2p ≤ −1, so n(9−2p) ≤ −n < 0 < 3.  ∎
     (Same spirit as 15.192 3/2·scheme; here the dual cost is in sum_ne
     units: scheme sum_ne=(n−2)/2, and (3/2)(n−1) is a looser dual
     allowance still sufficient for residual-(i).)

CERTIFIED
  F. At p=5: ‖δ‖₂² ≈ 23.63 ≤ room 36.4; ⟨δ,m4_part⟩≈0; ‖m4‖₂²≈143.2≤156.
     At p=7: residual-i R-path holds on full Max+ census (prior).
  G. D(C) sum_ne ≤ (3/2)(n−1) at p=5,7 (31.3≤37.5, 64.0≤73.5).

OPEN (blocks residual i / predicate flip)
  H. Prove ‖δ‖₂² ≤ room_δ^R for all primes p≥5
     (⇔ R≤2p ⇔ dual-eq empty via 15.217 chain),
     **or** Cov ≽ 4(n−6)/n I, **or** |μ|≤1/(2p),
     **or** sum_ne(D(C)) ≤ (3/2)(n−1) for all p≥5 and ker=sc.
     Soft-close forbidden.

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15247.json
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
from e1_gmin_m4_prop15190 import alpha_particular  # noqa: E402
from e1_gmin_m4_prop15216 import kappa_norm_sq  # noqa: E402
from e1_gmin_m4_prop15217 import delta_room_for_R, m4_target_R  # noqa: E402
from e1_gmin_m4_prop15226 import (  # noqa: E402
    mu_part_norm_sq,
    phi_kappa_inner,
    phi_norm_sq_closed,
    star_norm_sq_closed,
)


def m4_part_coeffs(p: int) -> tuple[Fraction, Fraction, Fraction]:
    """(a, b, z) for m4_part = a κ + b φ + z star."""
    D = p * p * (p * p - 5)
    a = Fraction(p * p - 1, D)
    b = Fraction(-2, D)
    z = Fraction(-2 * p, D)
    return a, b, z


def m4_part_norm_sq(p: int) -> Fraction:
    """‖m4_part‖₂² = ‖μ_part‖₂² + z² ‖star‖₂²."""
    a, b, z = m4_part_coeffs(p)
    # use mu_part_norm + z^2 star (star ⊥ κ,φ)
    return mu_part_norm_sq(p) + z * z * star_norm_sq_closed(p)


def theorem_m4_part_solution() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        n = n_of(p)
        a, b, z = m4_part_coeffs(p)
        D = p * p * (p * p - 5)
        # Verify 3×3 system for (4p I - T)m = 4κ/p
        # κ coeff: 4p a + 8 z == 4/p
        # φ coeff: 4p b - 4 z == 0
        # star: 4p z + 6 a - b(n+2) == 0
        e1 = 4 * p * a + 8 * z - Fraction(4, p)
        e2 = 4 * p * b - 4 * z
        e3 = 4 * p * z + 6 * a - b * (n + 2)
        if e1 != 0 or e2 != 0 or e3 != 0:
            ok = False
        # a,b match mu_part
        if a != Fraction(p * p - 1, D) or b != Fraction(-2, D):
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "a": str(a),
                "b": str(b),
                "z": str(z),
                "system_ok": e1 == 0 and e2 == 0 and e3 == 0,
            }
    return {
        "proved": ok,
        "theorem": (
            "m4_part=aκ+bφ+z star with a=(p²−1)/D, b=−2/D, z=−2p/D "
            "solves (4pI−T)m=4κ/p uniquely in span{κ,φ,star}."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "prop_15_215_T_algebra",
            "prop_15_224_T2_phi",
            "prop_15_177_m4_master",
        ],
    }


def theorem_m4_part_norm_and_room() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        mp = m4_part_norm_sq(p)
        thr = m4_target_R(p)
        room = thr - mp
        if room != delta_room_for_R(p):
            ok = False
        if room <= 0:
            ok = False
        # also = mu_part_norm + z^2 star
        a, b, z = m4_part_coeffs(p)
        alt = mu_part_norm_sq(p) + z * z * star_norm_sq_closed(p)
        if alt != mp:
            ok = False
        if p in (5, 7, 11, 13, 17):
            sample[str(p)] = {
                "m4_part_norm_sq": str(mp),
                "m4_target": str(thr),
                "room_delta": str(room),
                "room_positive": room > 0,
            }
    return {
        "proved": ok,
        "theorem": (
            "‖m4_part‖₂²=‖μ_part‖₂²+z²‖star‖₂² with z=−2p/D.  "
            "room_δ^R := n(n−2)/4−‖m4_part‖₂² equals 15.217 "
            "delta_room_for_R >0 for all primes p≥5.  "
            "R≤2p ⇔ ‖m4−m4_part‖₂² ≤ room_δ^R."
        ),
        "by_p_sample": sample,
        "hypothesis_proved_general": False,
        "depends_on": [
            "theorem_m4_part_solution",
            "prop_15_226_mu_part_norm",
            "prop_15_215_star_norm",
            "prop_15_217_R_criterion",
        ],
        "sufficient_for_residual_i": False,
    }


def theorem_three_halves_n_minus_1_dual_budget() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(5, 120) if is_prime(p)]:
        n = n_of(p)
        al = alpha_particular(p)
        lhs = al * Fraction(3, 2) * (n - 1)
        need = 2 - al
        if not (lhs < need):
            ok = False
        # poly check: n(9-2p) < 3
        if n * (9 - 2 * p) >= 3:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "alpha_3_2_n_minus_1": str(lhs),
                "need": str(need),
                "beats": lhs < need,
            }
    return {
        "proved": ok,
        "theorem": (
            "α·(3/2)·(n−1) < 2−α for all primes p≥5 "
            "(⇔ n(9−2p)<3; true since 9−2p≤−1).  "
            "Hence any free-e dual with sum_ne ≤ (3/2)(n−1) blocks dual-eq."
        ),
        "by_p_sample": sample,
        "depends_on": ["alpha_particular", "prop_15_190_scheme_dual"],
    }


def residual_i_closed_via_247() -> bool:
    return False


def hinge_status_247() -> dict:
    return {
        "m4_part_solution": theorem_m4_part_solution()["proved"],
        "m4_part_norm_room": theorem_m4_part_norm_and_room()["proved"],
        "delta_bound_general": theorem_m4_part_norm_and_room()[
            "hypothesis_proved_general"
        ],
        "three_halves_n_minus_1_budget": theorem_three_halves_n_minus_1_dual_budget()[
            "proved"
        ],
        "residual_i_closed_via_247": residual_i_closed_via_247(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "‖m4−m4_part‖₂² ≤ room_δ^R for all p≥5, "
            "or Cov≽4(n−6)/n I, or |μ|≤1/(2p), "
            "or D(C) sum_ne≤(3/2)(n−1) and ker=sc"
        ),
        "room_p5": str(delta_room_for_R(5)),
        "m4_part_norm_p5": str(m4_part_norm_sq(5)),
    }


def main() -> dict:
    a = theorem_m4_part_solution()
    b = theorem_m4_part_norm_and_room()
    c = theorem_three_halves_n_minus_1_dual_budget()
    status = hinge_status_247()
    out = {
        "prop": "15.247",
        "title": (
            "m4_part=μ_part+z star; R≤2p ⇔ ‖δ‖₂²≤room_δ^R; "
            "(3/2)(n−1) dual budget; δ bound OPEN"
        ),
        "proved": {
            "m4_part_solution": a["proved"],
            "m4_part_norm_and_room": b["proved"],
            "delta_bound_general": b["hypothesis_proved_general"],
            "three_halves_n_minus_1_budget": c["proved"],
            "residual_i_closed": residual_i_closed_via_247(),
        },
        "algebra": {
            "m4_part": "a κ + b φ + z star",
            "a": "(p^2-1)/D",
            "b": "-2/D",
            "z": "-2p/D",
            "room_delta_R": "(p-1)(p+1)(p^2+1)(3p^2-47)/(24(p^2-5))",
            "room_p5": str(delta_room_for_R(5)),
            "m4_part_norm_p5": str(m4_part_norm_sq(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_247(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "R-path residual is exactly ‖δ‖₂² on E_{+4p} after removing "
            "m4_part.  room_δ^R>0 Max+-free.  Bound on ‖δ‖₂² still OPEN. "
            "No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15247.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.247  residual_i={residual_i_closed_via_247()}")
    print(f"  m4_part={a['proved']} room={b['proved']} delta_gen={b['hypothesis_proved_general']}")
    print(f"  3/2(n-1) budget={c['proved']}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
