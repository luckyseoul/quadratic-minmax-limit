#!/usr/bin/env python3
"""
Prop 15.222 — ⟨D,TD⟩ closed form; CS bound ‖δ‖₂² ≤ k₂/|Max+|; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (from 15.221)
  y∈Max+={boolean Cy=py}, χ_S=∏_{i∈S} y_i, q_S=∑_{i≠j∈S} F_{ij} with
  F_{ij}=C_{ij} y_i y_j (so F1=p1, F²=p²I, F_ii=0, F_ij=±1).
  D_S=χ_S q_S, (4pI−T)χ=D, ‖D‖₂² known (15.221).
  m₄=E[χ]=m₄^{mn}+δ with δ∈E_{4p}; χ=χ_part+k, k∈E_{4p},
  m0=‖χ_part‖₂², k₂=binom(n,4)−m0 (constant on Max+ when m0 is).

PROVED Max+-free (any conference: Cᵀ=C, C_ii=0, C_ij=±1, C²=p²I; y∈{±1}ⁿ, Cy=py)
  A. Expansion of ⟨D,TD⟩ over common 3-sets.
     For S={v}∪A, S'={r}∪A with |A|=3, v≠r outside A:
       q_S=q_A+2s(v,A), q_{S'}=q_A+2s(r,A), and
       ⟨D,TD⟩=∑_{v≠r}∑_A F_{vr}(q_A+2s(v,A))(q_A+2s(r,A))
             =T1+T2+T3+T4 with
       T1=∑ F_{vr} q_A²,  T2=2∑ F_{vr} q_A s(v,A),
       T3=2∑ F_{vr} q_A s(r,A),  T4=4∑ F_{vr} s(v,A)s(r,A).

  B. Triangle sums on 3-sets A (σ_A=∑_{edges of A} F_e, q_A=2σ_A, π_A=∏_{e⊂A} F_e):
       ∑_A q_A = n p (n−2),
       ∑_A q_A² = 12 binom(n,3)   (wedge products of F vanish: ∑_i(p²−(n−1))=0),
       ∑_A π_A = 0                (Tr(F³)=0 ⇒ 6∑π=0),
       ∑_A σ_A³ = 7 ∑_A σ_A       (σ³=7σ+6π on each triangle),
       ∑_A q_A³ = 28 n p (n−2).

  C. T1. For fixed A: ∑_{v≠r∉A} F_{vr}=p(n−6)+q_A.
       T1=p(n−6)∑q_A²+∑q_A³=2np(n−2)(n²−7n+20).

  D. T2=T3. ∑_{r∉A,r≠v} F_{vr}=p−s(v,A), so
       T2=2∑_A q_A[p(3p−q_A)−∑_{v∉A}s(v,A)²].
       ∑_{v∉A}s²=3(n−3)−2P_A with ∑_A P_A=0 and P_A=(σ²−3)/2;
       q_A P_A=σ³−3σ ⇒ ∑ q_A P_A=2np(n−2).
       Using n−1=p²: T2=T3=4np(n−2)(6−n).

  E. T4. z_v=s(v,A) on Aᶜ gives zᵀFz=σ³−(2n−1)σ on each A;
       T4=4∑(σ³−(2n−1)σ)=4np(n−2)(4−n).

  F. Closed form.
       ⟨D,TD⟩=T1+2T2+T4=2np(n−2)(n²−13n+52)
       =2p(p−1)(p+1)(p²+1)(p⁴−11p²+40)   (n=p²+1).

  G. Cauchy–Schwarz averaging for δ.
       δ=E_Max+[k], k=χ−χ_part∈E_{4p}.
       ‖∑_y k(y)‖₂² ≤ |Max+| ∑_y ‖k(y)‖₂², hence
       ‖δ‖₂² ≤ E[‖k‖₂²]/|Max+|.
       If m0 (hence k₂=binom−m0) is constant on Max+:
       ‖δ‖₂² ≤ k₂/|Max+|.
       Therefore k₂ ≤ room·|Max+| ⇒ ‖δ‖₂²≤room ⇒ R≤2p (15.217)
       ⇒ residual (i) dual-eq empty.

CERTIFIED
  H. Formula F matches full T-matvec at p=3,5 (and single-y T-build at p=7).
  I. At p=5: m0=8038, k₂=6912, |Max+|=260, k₂/|Max+|≈26.58 ≤ room=182/5=36.4;
     actual ‖δ‖₂²≈23.63 ≤ room (strict CS slack).
  J. At p=7: m0≈36595, k₂≈193705, |Max+|=11452, k₂/|Max+|≈16.91 ≤ room=2500/11.

OPEN (blocks residual i / predicate flip)
  K. Max+-free lower bound m0 ≥ binom(n,4) − room(p)·|Max+|
     (equivalently upper bound k₂ ≤ room·|Max+|), including closed forms
     or sufficient bounds for m0 and |Max+| for all primes p≥5.
     (CS step G is proved; the scalar comparison is not.)

Until K: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15222.json
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
from e1_gmin_m4_prop15221 import D_energy_closed  # noqa: E402


def dtd_closed(p: int) -> Fraction:
    """⟨D,TD⟩ = 2p(p−1)(p+1)(p²+1)(p⁴−11p²+40)."""
    return Fraction(
        2 * p * (p - 1) * (p + 1) * (p * p + 1) * (p**4 - 11 * p * p + 40)
    )


def dtd_closed_via_n(p: int) -> Fraction:
    """⟨D,TD⟩ = 2 n p (n−2)(n²−13n+52) with n=p²+1."""
    n = n_of(p)
    return Fraction(2 * n * p * (n - 2) * (n * n - 13 * n + 52))


def T1_closed(p: int) -> Fraction:
    """T1 = 2 n p (n−2)(n²−7n+20)."""
    n = n_of(p)
    return Fraction(2 * n * p * (n - 2) * (n * n - 7 * n + 20))


def T2_closed(p: int) -> Fraction:
    """T2 = T3 = 4 n p (n−2)(6−n)."""
    n = n_of(p)
    return Fraction(4 * n * p * (n - 2) * (6 - n))


def T4_closed(p: int) -> Fraction:
    """T4 = 4 n p (n−2)(4−n)."""
    n = n_of(p)
    return Fraction(4 * n * p * (n - 2) * (4 - n))


def mean_T_on_D(p: int) -> Fraction:
    """⟨D,TD⟩/‖D‖₂²."""
    return dtd_closed(p) / D_energy_closed(p)


def ADA_energy(p: int) -> Fraction:
    """⟨D,(4pI−T)D⟩ = 4p‖D‖₂² − ⟨D,TD⟩."""
    return 4 * p * D_energy_closed(p) - dtd_closed(p)


def jensen_m0_lower(p: int) -> Fraction:
    """
    Convex f(λ)=1/(4p−λ)² on (−∞,4p): m0 ≥ ‖D‖₂² / (4p − λ̄)²
    with λ̄=⟨D,TD⟩/‖D‖₂² (Jensen).  Too weak alone for room at p=5.
    """
    D2 = D_energy_closed(p)
    lam = mean_T_on_D(p)
    gap = 4 * p - lam
    return D2 / (gap * gap)


def theorem_dtd_closed_form() -> dict:
    """Max+-free closed form for ⟨D,TD⟩ via T1–T4."""
    ok = True
    sample = {}
    for p in [pp for pp in range(3, 60) if is_prime(pp)]:
        a = dtd_closed(p)
        b = dtd_closed_via_n(p)
        tsum = T1_closed(p) + 2 * T2_closed(p) + T4_closed(p)
        if a != b or a != tsum or a <= 0:
            ok = False
        if p in (3, 5, 7, 11, 13):
            sample[str(p)] = {
                "dtd": str(a),
                "T1": str(T1_closed(p)),
                "T2": str(T2_closed(p)),
                "T4": str(T4_closed(p)),
                "dtd_over_D2": str(mean_T_on_D(p)),
                "ADA": str(ADA_energy(p)),
            }
    return {
        "proved": ok,
        "theorem": (
            "For every y∈Max+ and conference of order n=p²+1: "
            "⟨D,TD⟩=2p(p−1)(p+1)(p²+1)(p⁴−11p²+40).  "
            "Proof: expand over common 3-sets A as T1+T2+T3+T4; "
            "triangle F-algebra (wedge sum 0, Tr(F³)=0, σ³=7σ+6π) "
            "closes each Ti; assemble and substitute n=p²+1."
        ),
        "formula": "2p(p-1)(p+1)(p^2+1)(p^4-11p^2+40)",
        "formula_n": "2np(n-2)(n^2-13n+52)",
        "by_p_sample": sample,
        "depends_on": [
            "prop_15_221_D_energy",
            "F1_eq_p1",
            "F_squared_eq_p2_I",
            "T_definition_replacement",
        ],
    }


def theorem_CS_delta_le_k2_over_N() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Let k(y)=χ(y)−χ_part(y)∈E_{4p} and δ=E_{Max+}[k].  "
            "By Cauchy–Schwarz on ∑_y k(y) in the Hilbert space R^{binom(n,4)}: "
            "‖δ‖₂² ≤ E[‖k‖₂²]/|Max+|.  "
            "If ‖k(y)‖₂²=k₂ is constant on Max+, then ‖δ‖₂² ≤ k₂/|Max+|."
        ),
        "conditional_close": (
            "k₂ ≤ room(p)·|Max+| with room=δ-room of 15.217 "
            "⇒ ‖δ‖₂²≤room ⇒ ‖m₄‖₂²≤n(n−2)/4 ⇒ residual (i)."
        ),
        "depends_on": ["m4_equals_E_chi", "delta_in_E_4p"],
    }


def theorem_m0_k2_bound_open() -> dict:
    room_ok = all(
        delta_room_for_R(p) > 0 for p in range(5, 40) if is_prime(p)
    )
    # Jensen gap at p=5: need m0 ≥ binom - room*|Max+| with |Max+|=260
    # = 14950 - 36.4*260 = 5486; Jensen only ~2054
    n5 = n_of(5)
    need5 = Fraction(comb(n5, 4)) - delta_room_for_R(5) * 260
    jen5 = jensen_m0_lower(5)
    return {
        "proved": False,
        "room_positive_p_ge_5": room_ok,
        "theorem": (
            "OPEN: Max+-free m0=‖χ_part‖₂² (or lower bound) and |Max+| "
            "(or lower bound) so that binom(n,4)−m0 ≤ room(p)·|Max+| "
            "for all primes p≥5.  CS (G) is proved; comparison open."
        ),
        "jensen_insufficient_p5": {
            "jensen_m0_lb": str(jen5),
            "need_m0_if_Nmax_260": str(need5),
            "jensen_lt_need": jen5 < need5,
        },
        "certified_small_p": {
            "p5": {
                "m0": 8038,
                "k2": 6912,
                "Nmax": 260,
                "k2_over_N": str(Fraction(6912, 260)),
                "room": str(delta_room_for_R(5)),
                "CS_closes_numerically": Fraction(6912, 260) <= delta_room_for_R(5),
            },
            "p7": {
                "m0": 36595,
                "k2": 193705,
                "Nmax": 11452,
                "k2_over_N": str(Fraction(193705, 11452)),
                "room": str(delta_room_for_R(7)),
                "CS_closes_numerically": (
                    Fraction(193705, 11452) <= delta_room_for_R(7)
                ),
            },
        },
        "depends_on": [
            "theorem_dtd_closed_form",
            "theorem_CS_delta_le_k2_over_N",
            "prop_15_217_room",
        ],
    }


def residual_i_closed_via_222() -> bool:
    return False


def hinge_status_222() -> dict:
    return {
        "dtd_closed": theorem_dtd_closed_form()["proved"],
        "CS_delta": theorem_CS_delta_le_k2_over_N()["proved"],
        "m0_k2_bound": theorem_m0_k2_bound_open()["proved"],
        "residual_i_closed_via_222": residual_i_closed_via_222(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "m0 ≥ binom − room·|Max+| Max+-free (with |Max+| LB), "
            "or |μ|≤1/(2p), or K₄ thr, or free-e on true ker."
        ),
    }


def main() -> dict:
    a = theorem_dtd_closed_form()
    b = theorem_CS_delta_le_k2_over_N()
    c = theorem_m0_k2_bound_open()
    status = hinge_status_222()
    out = {
        "prop": "15.222",
        "title": (
            "⟨D,TD⟩=2p(p−1)(p+1)(p²+1)(p⁴−11p²+40); "
            "‖δ‖₂²≤k₂/|Max+| by CS; m0-bound OPEN"
        ),
        "proved": {
            "dtd_closed": a["proved"],
            "CS_delta": b["proved"],
            "m0_k2_bound": c["proved"],
            "residual_i_closed": residual_i_closed_via_222(),
        },
        "algebra": {
            "dtd_formula": str(dtd_closed(5)),
            "dtd_by_p": {
                str(p): str(dtd_closed(p)) for p in (3, 5, 7, 11, 13)
            },
            "T_pieces_p5": {
                "T1": str(T1_closed(5)),
                "T2": str(T2_closed(5)),
                "T4": str(T4_closed(5)),
            },
            "m0_path": c,
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_222(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "dtd and CS averaging proved Max+-free.  "
            "m0/|Max+| comparison still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15222.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.222  residual_i={residual_i_closed_via_222()}")
    print(f"  dtd={a['proved']} CS={b['proved']} m0_bound={c['proved']}")
    print(f"  dtd(p=5)={dtd_closed(5)}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
