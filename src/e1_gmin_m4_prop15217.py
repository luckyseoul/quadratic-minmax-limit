#!/usr/bin/env python3
"""
Prop 15.217 — Φ identity on Max+; residual-(i) reduces to m4₂ ≤ n(n−2)/4.

Repairs the 15.216 Rayleigh gap (convex combo of R_ke and 4p) by computing
⟨m₄,κ⟩ exactly on Max+={y∈{±1}ⁿ : Cy=py}.

PROVED Max+-free (any conference C with Cᵀ=C, C_ii=0, C_ij=±1, C²=(n−1)I,
n=p²+1, and Max+ = boolean +p-eigenvectors):
  A. Φ(y) := ∑_{|S|=4} κ(S) ∏_{i∈S} y_i equals n(n−1)(n−2)/8 for every y∈Max+.
  B. Hence ⟨m₄,κ⟩ = Φ = n(n−1)(n−2)/8 (m₄ = E_Max+[∏ y]).
  C. Master (4pI−T)m₄=4κ/p ⇒ R := ⟨m₄,Tm₄⟩/‖m₄‖₂² = 4p − 4Φ/(p ‖m₄‖₂²).
     R≤2p ⇔ ‖m₄‖₂² ≤ n(n−2)/4.
  D. With Φ exact: ‖η‖₂² = ‖m₄‖₂² − n(n−2)/4 + ‖κ‖₂²/p⁴.
     So R≤2p ⇔ ‖η‖₂² ≤ ‖κ‖₂²/p⁴ (crude) ⇔ (with crude≤thr-η-budget) K₄≤thr
     ⇔ Q_e≤10 ⇔ dual-eq empty (15.216 reductions B,C,D that remain valid).
  E. Equivalent δ-form: m₄=ke+δ, δ∈E_{4p}, ⟨ke,δ⟩=0 ⇒
        R≤2p ⇔ ‖δ‖₂² ≤ n(n−2)/4 − ‖ke‖₂²
     with closed room (p−1)(p+1)(p²+1)(3p²−47)/(24(p²−5)) >0 for p≥5.

OPEN (blocks residual i / predicate flip):
  F. Prove ‖m₄‖₂² ≤ n(n−2)/4 for all primes p≥5 (equivalently K₄ ≤ n(15n−22),
     or ‖δ‖₂² ≤ room above).  Census: holds at p=5 (m4₂≈143.2 ≤ 156).
     Path-C residual δ²≤room_hyp/24 is sufficient (room_hyp/24 ≤ R-room) but open.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

Writes evidence/e1_gmin_m4_prop15217.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15197 import K4_threshold_for_Q10  # noqa: E402
from e1_gmin_m4_prop15198 import C_diag  # noqa: E402
from e1_gmin_m4_prop15215 import eta_part_norm_sq  # noqa: E402
from e1_gmin_m4_prop15216 import (  # noqa: E402
    R_ke,
    crude_eta_ub,
    kappa_norm_sq,
    residual_i_dual_eq_empty_proved_general,
    thr_eta_budget,
    theorem_Q10_blocks_dual_eq_all_p,
    theorem_crude_le_thr_eta_budget,
)


# ---------------------------------------------------------------------------
# A. Φ identity
# ---------------------------------------------------------------------------


def Phi_closed(p: int) -> Fraction:
    """n(n−1)(n−2)/8 with n=p²+1."""
    n = n_of(p)
    return Fraction(n * (n - 1) * (n - 2), 8)


def theorem_Phi_identity() -> dict:
    """
    For conference C (Cᵀ=C, zero diagonal, C_ij=±1, C²=(n−1)I) and
    y∈{±1}ⁿ with Cy=λy, λ²=n−1:

      Q := yᵀ C y = λ n.
      Define M_ij = C_ij y_i y_j (i≠j), M_ii=0.  Then M is symmetric,
      M_ij=±1 off-diagonal, M1=λ1 (row sums λ), and Q=1ᵀM1=λn.

      Expand Q² = ∑_{i≠j} ∑_{k≠l} M_ij M_kl by support size r=|{i,j,k,l}|:

      r=2: only (i,j)∈{(k,l),(l,k)}.  Contribution 2n(n−1)
            (n(n−1) ordered edges, each pairs with itself and reverse).

      r=3: all vanish when λ²=n−1:
        • V-shape same source: ∑_i ∑_{j≠k≠i} M_ij M_ik
          = ∑_i[(∑_{j≠i} M_ij)² − ∑_{j≠i} M_ij²] = ∑_i(λ² − (n−1)) = 0.
        • V-shape same target: identical by M=Mᵀ.
        • Paths of length 2: ∑_{i,j,k dist} M_ij M_jk
          = ∑_{i≠k} y_i y_k ∑_{j≠i,k} C_ij C_jk = ∑_{i≠k} y_i y_k (C²)_{ik}=0.

      r=4: for each 4-set S and each of the three matchings, 8 ordered
            directed edge pairs contribute M_e M_f = κ_matching(S)·∏_{i∈S} y_i.
            Total 8 Φ(y) with Φ=∑_S κ(S)∏ y_i.

      Therefore Q² = 2n(n−1) + 8 Φ, so
        Φ = (λ² n² − 2n(n−1))/8 = n(n−1)(n−2)/8
      when λ²=n−1.  On Max+ (λ=p) this is constant, independent of y.  ∎
    """
    # Algebra check: Φ formula matches (p² n² - 2n(n-1))/8
    ok = True
    sample = {}
    for p in [p for p in range(3, 60) if is_prime(p)]:
        n = n_of(p)
        Phi = Phi_closed(p)
        alt = Fraction(p * p * n * n - 2 * n * (n - 1), 8)
        if Phi != alt or Phi != Fraction(n * (n - 1) * (n - 2), 8):
            ok = False
        if p in (3, 5, 7, 11):
            sample[str(p)] = {"Phi": str(Phi), "match_alt": Phi == alt}
    # Census: Φ constant on Max+ at p=3,5 (from prior session / identity)
    return {
        "proved": ok,
        "theorem": (
            "For every y∈Max+={boolean Cy=py}: Φ(y)=∑_S κ(S)∏_{i∈S} y_i "
            "equals n(n−1)(n−2)/8.  Proof: Q=yᵀCy=pn; Q²=2n(n−1)+8Φ via "
            "M_ij=C_ij y_i y_j (row sum p, C² off-diag 0 kills r=3).  "
            "Hence ⟨m₄,κ⟩=Φ constant."
        ),
        "formula": "n(n-1)(n-2)/8",
        "by_p_sample": sample,
        "depends_on": ["C_squared_eq_n_minus_1_I", "Max_plus_boolean_evecs"],
    }


# ---------------------------------------------------------------------------
# B–D. R≤2p ⇔ m4₂ ≤ n(n−2)/4 ⇔ η ≤ crude
# ---------------------------------------------------------------------------


def m4_target_R(p: int) -> Fraction:
    """n(n−2)/4 = 2Φ/p²."""
    n = n_of(p)
    return Fraction(n * (n - 2), 4)


def delta_room_for_R(p: int) -> Fraction:
    """‖δ‖₂² upper for R≤2p: n(n−2)/4 − ‖ke‖₂²."""
    # Closed form (p−1)(p+1)(p²+1)(3p²−47)/(24(p²−5))
    return Fraction(
        (p - 1) * (p + 1) * (p * p + 1) * (3 * p * p - 47),
        24 * (p * p - 5),
    )


def ke2_of(p: int) -> Fraction:
    n = n_of(p)
    Tk2 = 6 * n * (n - 1) * (n - 2) * (n - 6)
    k2 = kappa_norm_sq(p)
    ep2 = eta_part_norm_sq(p)
    kap_eta = Fraction(Tk2, p * p * 12 * (p * p - 5))
    return Fraction(k2, p**4) + ep2 + 2 * kap_eta / Fraction(p * p)


def theorem_R_le_2p_iff_m4_bound() -> dict:
    """
    Master (4pI−T)m₄=4κ/p ⇒ ⟨m₄,(4pI−T)m₄⟩=4⟨m₄,κ⟩/p.
    With ⟨m₄,κ⟩=Φ: R=4p−4Φ/(p‖m₄‖₂²).
    R≤2p ⇔ Φ ≥ (p²/2)‖m₄‖₂² ⇔ ‖m₄‖₂² ≤ 2Φ/p² = n(n−2)/4
    (using p²=n−1).  ∎
    """
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        Phi = Phi_closed(p)
        tgt = m4_target_R(p)
        # 2Φ/p² == n(n-2)/4
        if 2 * Phi / Fraction(p * p) != tgt:
            ok = False
        room = delta_room_for_R(p)
        # room == tgt - ke2
        if room != tgt - ke2_of(p):
            ok = False
        # room > 0 for p≥5
        if room <= 0:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "Phi": str(Phi),
                "m4_target": str(tgt),
                "ke2": str(ke2_of(p)),
                "delta_room": str(room),
                "R_ke": str(R_ke(p)),
            }
    return {
        "proved": ok,
        "theorem": (
            "For all primes p≥5: R(m₄)≤2p ⇔ ‖m₄‖₂² ≤ n(n−2)/4 "
            "⇔ ‖δ‖₂² ≤ (p−1)(p+1)(p²+1)(3p²−47)/(24(p²−5)), "
            "where m₄=ke+δ, δ∈E_{4p}.  (Uses Φ identity A.)"
        ),
        "by_p_sample": sample,
        "depends_on": ["theorem_Phi_identity", "prop_15_177_master"],
    }


def theorem_eta_formula_with_Phi() -> dict:
    """‖η‖₂² = ‖m₄‖₂² − n(n−2)/4 + ‖κ‖₂²/p⁴."""
    return {
        "proved": True,
        "theorem": (
            "η=m₄−κ/p² ⇒ ‖η‖₂²=‖m₄‖₂²−2⟨m₄,κ⟩/p²+‖κ‖₂²/p⁴.  "
            "With ⟨m₄,κ⟩=Φ=n(n−1)(n−2)/8 and 2Φ/p²=n(n−2)/4: "
            "‖η‖₂²=‖m₄‖₂²−n(n−2)/4+‖κ‖₂²/p⁴.  "
            "Hence R≤2p (‖m₄‖₂²≤n(n−2)/4) ⇒ ‖η‖₂²≤‖κ‖₂²/p⁴."
        ),
        "depends_on": ["theorem_Phi_identity", "theorem_R_le_2p_iff_m4_bound"],
    }


def theorem_R_implies_K4_thr() -> dict:
    """Conditional: m4₂≤n(n−2)/4 ⇒ K₄≤16n²−14n ⇒ dual-eq empty."""
    crude_ok = theorem_crude_le_thr_eta_budget()
    q10_ok = theorem_Q10_blocks_dual_eq_all_p()
    return {
        "proved": bool(crude_ok["proved"] and q10_ok["proved"]),
        "conditional_on": "‖m₄‖₂² ≤ n(n−2)/4 (OPEN in general)",
        "theorem": (
            "IF ‖m₄‖₂² ≤ n(n−2)/4 for all primes p≥5, THEN "
            "‖η‖₂²≤‖κ‖₂²/p⁴≤thr-η-budget ⇒ K₄≤16n²−14n ⇒ Q_e≤10 "
            "⇒ dual-eq empty (15.216 D) ⇒ residual (i) Type I closed."
        ),
        "depends_on": [
            "theorem_eta_formula_with_Phi",
            "theorem_crude_le_thr_eta_budget",
            "theorem_Q10_blocks_dual_eq_all_p",
        ],
    }


# ---------------------------------------------------------------------------
# F. OPEN hinge: m4₂ bound
# ---------------------------------------------------------------------------


def theorem_m4_bound_open() -> dict:
    """The remaining Max+-free step."""
    # Algebra: K4 target equivalent
    # m4₂ ≤ n(n−2)/4 ⇔ K₄ ≤ n(15n−22) (via e4 identity + E[s²]=2n)
    ok_equiv = True
    sample = {}
    for p in [5, 7, 11, 13, 17]:
        n = n_of(p)
        # K4_need = n(15n - 22)
        k4_need = n * (15 * n - 22)
        thr = K4_threshold_for_Q10(p)
        # need < thr: thr - need = (p²+1)(p²+9)
        diff = thr - k4_need
        expected_diff = (p * p + 1) * (p * p + 9)
        if diff != expected_diff:
            ok_equiv = False
        sample[str(p)] = {
            "K4_need_for_R": k4_need,
            "K4_thr_Q10": thr,
            "thr_minus_need": diff,
            "formula_diff": expected_diff,
        }
    # room vs Path-C
    pathc_suff = True
    for p in [p for p in range(5, 60) if is_prime(p)]:
        # room_hyp/24 = 4(p²-9)(p²-1)² / (3(p²-5)(p²+1))
        rh = Fraction(
            4 * (p * p - 9) * (p * p - 1) ** 2,
            3 * (p * p - 5) * (p * p + 1),
        )
        if rh > delta_room_for_R(p):
            pathc_suff = False
    return {
        "proved": False,
        "theorem": (
            "OPEN: ‖m₄‖₂² ≤ n(n−2)/4 for all primes p≥5.  Equivalent forms: "
            "K₄ ≤ n(15n−22); ‖δ‖₂² ≤ (p−1)(p+1)(p²+1)(3p²−47)/(24(p²−5)).  "
            "Path-C residual δ²≤room_hyp/24 is sufficient (room_hyp/24 ≤ R-room "
            f"for all primes p≥5: {pathc_suff}) but itself OPEN.  "
            "Census: m4₂ ≤ target at p=5 (boolean evec enumeration)."
        ),
        "K4_equivalence_algebra_ok": ok_equiv,
        "pathc_residual_sufficient": pathc_suff,
        "by_p_sample": sample,
        "depends_on": ["theorem_R_le_2p_iff_m4_bound", "E_s2_eq_2n_2design"],
    }


def residual_i_closed_via_217() -> bool:
    """False until m4₂ bound proved."""
    return False


def hinge_status_217() -> dict:
    return {
        "Phi_identity": theorem_Phi_identity()["proved"],
        "R_iff_m4_bound": theorem_R_le_2p_iff_m4_bound()["proved"],
        "eta_formula_Phi": theorem_eta_formula_with_Phi()["proved"],
        "R_implies_K4_thr_conditional": theorem_R_implies_K4_thr()["proved"],
        "m4_bound": theorem_m4_bound_open()["proved"],
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "residual_i_closed_via_217": residual_i_closed_via_217(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "note": (
            "15.217 repairs Rayleigh gap: R≤2p is equivalent to m4₂≤n(n−2)/4 "
            "(via exact Φ).  That m4₂ bound remains OPEN.  Residual (i) OPEN."
        ),
    }


def main() -> dict:
    a = theorem_Phi_identity()
    b = theorem_R_le_2p_iff_m4_bound()
    c = theorem_eta_formula_with_Phi()
    d = theorem_R_implies_K4_thr()
    e = theorem_m4_bound_open()
    status = hinge_status_217()
    out = {
        "prop": "15.217",
        "title": "Φ identity; residual-(i) reduces to m4₂ ≤ n(n−2)/4 (OPEN)",
        "proved": {
            "Phi_identity": a["proved"],
            "R_iff_m4_bound": b["proved"],
            "eta_formula_with_Phi": c["proved"],
            "R_implies_dual_eq_conditional": d["proved"],
            "m4_bound": e["proved"],
            "residual_i_closed": residual_i_closed_via_217(),
        },
        "open_hinge": e["theorem"],
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_217(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "soft_close": False,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15217.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {path}")
    for k, v in out["proved"].items():
        print(f"  {k}={v}")
    print(f"  residual_i={residual_i_closed_via_217()}")
    return out


if __name__ == "__main__":
    main()
