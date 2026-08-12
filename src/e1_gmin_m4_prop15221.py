#!/usr/bin/env python3
"""
Prop 15.221 — (4pI−T)χ=D on Max+; ‖D‖₂² closed form; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP
  For y∈Max+={boolean Cy=py}, χ∈R^{binom(n,4)} with χ_S=∏_{i∈S} y_i.
  T = C-weighted 4-set replacement operator (15.177).
  q_S(y)=∑_{i≠j∈S} C_{ij} y_i y_j,   D_S=χ_S q_S.
  m₄=E_Max+[χ]; master (4pI−T)m₄=4κ/p; m₄=m₄^{mn}+δ, δ∈E_{4p}.

PROVED Max+-free (any conference C: Cᵀ=C, C_ii=0, C_ij=±1, C²=(n−1)I=p²I;
y∈{±1}ⁿ with Cy=py)
  A. Identity (4p I − T) χ = D.
     Expand (Tχ)_S=∑_{v∈S,r∉S} C_{vr} χ_{S_{v→r}}.
     ∑_{r∉S} C_{vr} y_r = p y_v − ∑_{w∈S∖{v}} C_{vw} y_w.
     Hence (Tχ)_S = ∑_v (∏_{u∈S∖v} y_u)(p y_v − ∑_{w∈S∖v} C_{vw} y_w)
       = 4p χ_S − χ_S ∑_{v≠w∈S} C_{vw} y_v y_w
       = 4p χ_S − D_S.  ∎

  B. Wedge vanishing for edge features f_{ij}=C_{ij} y_i y_j.
     ∑_{i} ∑_{j<k, j,k≠i} f_{ij} f_{ik}
       = ½ ∑_i[(∑_{j≠i} f_{ij})² − ∑_{j≠i} f_{ij}²].
     ∑_{j≠i} f_{ij}=y_i(Cy)_i=p,  ∑_{j≠i} f_{ij}²=n−1=p².
     Hence each vertex contributes ½(p²−p²)=0; sum of f_e f_{e'} over
     all wedge (share-vertex) edge pairs is 0.  ∎

  C. Closed form ‖D‖₂² = ∑_S q_S².
     Let F=∑_{i<j} f_{ij}=pn/2,  N_e=binom(n,2).
     From B and (∑f_e)²: sum of f_e f_{e'} over disjoint edge pairs
       = (F² − N_e)/2.
     Double-count ∑_S(∑_{e⊂S} f_e)²
       = N_e·binom(n−2,2) + 2·0·(n−3) + 2·(F²−N_e)/2
       = N_e·binom(n−2,2) + F² − N_e.
     With q_S=2∑_{e⊂S} f_e:
       ∑_S q_S² = 4 N_e binom(n−2,2) + 4F² − 4 N_e
         = n(n−1)(n−2)(n−3) + p² n² − 2n(n−1).  ∎
     Constant on all of Max+ (independent of y).

  D. Consequences for the master residual δ.
     (A) ⇒ for each y, χ(y)=χ^{part}(y)+k(y) with k∈E_{4p}(T) when
     solving (4pI−T)x=D(y) with min-norm part.
     Averaging: m₄=m₄^{mn}+δ, δ=E[k], δ∈E_{4p}.
     ‖D‖₂² exact (C) is the starting point for spectral control of ‖k‖₂
     and hence ‖δ‖₂; the pointwise / L² bound on δ that would force
     ‖m₄‖₂²≤n(n−2)/4 (15.217) is still OPEN.

CERTIFIED
  E. Formula C matches full enum at p=3,5 (‖D‖₂²=5760, 374400).
  F. At p=5: ‖δ‖₂²≈23.63 ≤ room 36.4 (15.217); R(m₄)≈9.11<2p=10.

OPEN (blocks residual i / predicate flip)
  G. Prove Max+-free ‖δ‖₂² ≤ (p−1)(p+1)(p²+1)(3p²−47)/(24(p²−5))
     (the 15.217 room), or any other hinge |μ|≤1/(2p) / K₄ thr / free-e.
     Spectral path: use (C)+(A) and σ(T) to bound avg‖P_{E_{4p}} χ‖₂².

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15221.json
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
    residual_i_dual_eq_empty_proved_general,
)
from e1_gmin_m4_prop15217 import delta_room_for_R, m4_target_R  # noqa: E402


def D_energy_closed(p: int) -> Fraction:
    """∑_S q_S² = n(n−1)(n−2)(n−3) + p² n² − 2n(n−1)."""
    n = n_of(p)
    return Fraction(
        n * (n - 1) * (n - 2) * (n - 3) + p * p * n * n - 2 * n * (n - 1)
    )


def theorem_Tchi_eq_D() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For y∈Max+, χ_S=∏_{i∈S} y_i, q_S=∑_{i≠j∈S} C_{ij} y_i y_j, "
            "D_S=χ_S q_S: (4p I − T)χ = D.  "
            "Proof: ∑_{r∉S} C_{vr} y_r = p y_v − ∑_{w∈S∖v} C_{vw} y_w; "
            "expand Tχ and use ∏_{S∖v} y = χ y_v."
        ),
        "depends_on": ["Cy_eq_py", "T_definition_replacement"],
    }


def theorem_wedge_ff_vanishes() -> dict:
    return {
        "proved": True,
        "theorem": (
            "f_{ij}=C_{ij} y_i y_j: for each i, ∑_{j≠i} f_{ij}=p and "
            "∑_{j≠i} f_{ij}²=n−1=p², so ∑_{j<k≠i} f_{ij} f_{ik}=0.  "
            "Hence ∑_{wedge edge pairs} f_e f_{e'}=0."
        ),
        "depends_on": ["Cy_eq_py", "C_squared_eq_n_minus_1_I", "y_i_squared_1"],
    }


def theorem_D_energy_closed_form() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(3, 60) if is_prime(p)]:
        n = n_of(p)
        D2 = D_energy_closed(p)
        # rebuild from F, Ne
        F2 = Fraction(p * p * n * n, 4)
        Ne = Fraction(n * (n - 1), 2)
        b22 = Fraction((n - 2) * (n - 3), 2)
        alt = 4 * (Ne * b22 + F2 - Ne)
        if D2 != alt or D2 <= 0:
            ok = False
        if p in (3, 5, 7, 11, 13):
            sample[str(p)] = {
                "D2": str(D2),
                "binom_n_4": comb(n, 4),
                "D2_over_binom": float(D2 / comb(n, 4)),
            }
    return {
        "proved": ok,
        "theorem": (
            "For every y∈Max+ and every conference of order n=p²+1: "
            "∑_{|S|=4} q_S(y)² = n(n−1)(n−2)(n−3)+p²n²−2n(n−1).  "
            "Proof: wedge products of f vanish (B); disjoint products from "
            "(∑f)²; double-count 4-sets containing edge pairs."
        ),
        "formula": "n(n-1)(n-2)(n-3)+p^2 n^2-2n(n-1)",
        "by_p_sample": sample,
        "depends_on": ["theorem_wedge_ff_vanishes", "yTCy_eq_pn"],
    }


def theorem_delta_path_open() -> dict:
    """Record that D-energy enables but does not finish the δ bound."""
    room_ok = all(
        delta_room_for_R(p) > 0 for p in range(5, 40) if is_prime(p)
    )
    return {
        "proved": False,
        "room_positive_p_ge_5": room_ok,
        "theorem": (
            "OPEN: from (4pI−T)χ=D and exact ‖D‖₂², prove "
            "‖δ‖₂² ≤ room of 15.217 for m₄=E[χ], so ‖m₄‖₂²≤n(n−2)/4."
        ),
        "conditional_chain": (
            "‖δ‖₂²≤room ⇒ ‖m₄‖₂²≤n(n−2)/4 ⇒ R≤2p ⇒ K₄ thr ⇒ dual-eq empty "
            "(15.217/216) ⇒ residual_i."
        ),
        "depends_on": [
            "theorem_Tchi_eq_D",
            "theorem_D_energy_closed_form",
            "prop_15_217_room",
        ],
    }


def residual_i_closed_via_221() -> bool:
    return False


def hinge_status_221() -> dict:
    return {
        "Tchi_eq_D": theorem_Tchi_eq_D()["proved"],
        "wedge_ff_vanishes": theorem_wedge_ff_vanishes()["proved"],
        "D_energy_closed": theorem_D_energy_closed_form()["proved"],
        "delta_bound": theorem_delta_path_open()["proved"],
        "residual_i_closed_via_221": residual_i_closed_via_221(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "‖δ‖₂²≤15.217-room (via D-energy+σ(T)), or |μ|≤1/(2p), or K₄ thr, "
            "or free-e on true ker."
        ),
    }


def main() -> dict:
    a = theorem_Tchi_eq_D()
    b = theorem_wedge_ff_vanishes()
    c = theorem_D_energy_closed_form()
    d = theorem_delta_path_open()
    status = hinge_status_221()
    out = {
        "prop": "15.221",
        "title": (
            "(4pI−T)χ=D; ‖D‖₂²=n(n−1)(n−2)(n−3)+p²n²−2n(n−1); "
            "δ-bound OPEN"
        ),
        "proved": {
            "Tchi_eq_D": a["proved"],
            "wedge_ff_vanishes": b["proved"],
            "D_energy_closed": c["proved"],
            "delta_bound": d["proved"],
            "residual_i_closed": residual_i_closed_via_221(),
        },
        "algebra": {
            "D_energy_formula": str(D_energy_closed(5)),
            "D_energy_by_p": {
                str(p): str(D_energy_closed(p))
                for p in (3, 5, 7, 11, 13)
            },
            "delta_path": d,
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_221(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "D-energy and Tχ=D proved Max+-free.  "
            "δ²≤room still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15221.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.221  residual_i={residual_i_closed_via_221()}")
    print(f"  Tχ=D={a['proved']} wedge={b['proved']} D2={c['proved']}")
    print(f"  D2(p=5)={D_energy_closed(5)}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
