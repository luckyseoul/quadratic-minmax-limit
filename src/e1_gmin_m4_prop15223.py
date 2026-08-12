#!/usr/bin/env python3
"""
Prop 15.223 — Path-C residual ⇒ Wick/K₄ residual-(i); CS δ vs Wick/R rooms;
hinge still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP
  From 15.215: ∑_S η_S² = ‖η_part‖₂² + ‖δ‖₂² (δ∈E_{4p}, η_part ⊥ δ),
  with closed
      ‖η_part‖₂² = 5(p²−1)(p²+1)(p²+3)/(6 p² (p²−5)),
  and K₄≤Wick_hi ⇔ ∑η² ≤ B_wick := n(13n−10)/(6(n−1)).
  Path-C residual (15.117): δ² ≤ room_hyp/24 with
      room_hyp/24 = 4(p²−9)(p²−1)²/(3(p²−5)(p²+1)).
  15.217 R-room: δ² ≤ (p−1)(p+1)(p²+1)(3p²−47)/(24(p²−5)).
  15.222: ‖δ‖₂² ≤ k₂/|Max+| by CS (k₂=‖χ−χ_part‖₂² constant on Max+).

PROVED Max+-free
  A. Wick gap closed form.
        gap_wick(p) := B_wick − ‖η_part‖₂²
                     = 4(p⁴ − 8p² − 9)/(3(p²−5)).
     Proof.  n=p²+1, n−1=p²:
        B_wick = n(13n−10)/(6(n−1)) = (p²+1)(13p²+3)/(6p²).
        Subtract η_part (15.215 H) over common denominator 6p²(p²−5);
        numerator simplifies to 4p²(p⁴−8p²−9) / [6p²(p²−5)] wait
        direct Fraction algebra (module checks all primes ≤200):
        B − η = 4(p⁴−8p²−9)/(3(p²−5)).  ∎

  B. Path-C budget strictly inside Wick gap for all primes p≥5.
        gap_wick − room_hyp/24
          = 16 p² (p²−9) / [3(p²−5)(p²+1)]
          = 16 p² (p−3)(p+3) / [3(p²−5)(p²+1)].
     For p≥5: numerator >0 and denominator >0, so
        room_hyp/24 < gap_wick.
     Equality of the two budgets fails only at p=3 (both sides 0 for hyp;
     gap_wick−hyp=0).  ∎

  C. Conditional residual-(i) close via Path C.
     Path C (δ² ≤ room_hyp/24) + (B) ⇒ δ² < gap_wick
     ⇒ ∑η² = η_part + δ² < B_wick ⇒ K₄ ≤ Wick_hi (15.198)
     ⇒ Q_e < 10 ⇒ dual-eq empty for all primes p≥5 (15.198 A–B).
     Same Path C also ⇒ δ² ≤ room_hyp/24 ≤ R-room (15.217 Path-C
     sufficiency, already recorded) ⇒ ‖m₄‖₂² ≤ n(n−2)/4.  ∎

  D. CS comparison (structure).
     By 15.222 G: ‖δ‖₂² ≤ k₂/|Max+|.
     Hence either of
        k₂/|Max+| ≤ gap_wick    or    k₂/|Max+| ≤ R-room
     closes residual (i) (Wick path or 15.217 path).  ∎

CERTIFIED
  E. At p=5: k₂=6912, |Max+|=260, k₂/N=6912/260=1728/65≈26.585
     ≤ gap_wick=416/15≈27.733 and ≤ R-room=182/5=36.4;
     actual δ²=room_hyp/24=1536/65≈23.631.
  F. At p=7: k₂≈193705, |Max+|=11452, k₂/N≈16.91
     ≤ gap_wick=2000/33≈60.61 and ≤ R-room=2500/11≈227.27.
  G. Identity A,B verified as exact Fraction for all primes 5≤p≤300.

OPEN (blocks residual i / predicate flip)
  H. Prove Max+-free one of:
     (i)   δ² ≤ room_hyp/24 for all p≥5 (Path C primary — then C fires);
     (ii)  k₂/|Max+| ≤ gap_wick (or ≤ R-room), with closed forms / bounds
           for k₂=binom(n,4)−‖χ_part‖₂² and |Max+|;
     (iii) |μ|≤1/(2p), K₄≤Wick by other means, free-e+ker=sc, etc.
  Soft-close forbidden.  Census does not flip general predicates.

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15223.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
    residual_i_dual_eq_empty_proved_general,
)
from e1_gmin_m4_prop15198 import eta_sumsq_budget_for_wick_hi  # noqa: E402
from e1_gmin_m4_prop15215 import eta_part_norm_sq  # noqa: E402
from e1_gmin_m4_prop15217 import delta_room_for_R  # noqa: E402


def gap_wick_closed(p: int) -> Fraction:
    """B_wick − ‖η_part‖₂² = 4(p⁴ − 8p² − 9)/(3(p² − 5))."""
    return Fraction(4 * (p**4 - 8 * p * p - 9), 3 * (p * p - 5))


def hyp_inside_wick_gap(p: int) -> Fraction:
    """gap_wick − room_hyp/24 = 16 p² (p²−9) / [3(p²−5)(p²+1)]."""
    return Fraction(
        16 * p * p * (p * p - 9),
        3 * (p * p - 5) * (p * p + 1),
    )


def theorem_gap_wick_closed_form() -> dict:
    ok = True
    sample = {}
    for p in range(5, 201):
        if not is_prime(p):
            continue
        g = gap_wick_closed(p)
        alt = eta_sumsq_budget_for_wick_hi(p) - eta_part_norm_sq(p)
        if g != alt or g <= 0:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101, 199):
            sample[str(p)] = {
                "gap_wick": str(g),
                "match_B_minus_eta": g == alt,
            }
    return {
        "proved": ok,
        "theorem": (
            "For every prime p≥5: B_wick − ‖η_part‖₂² "
            "= 4(p⁴−8p²−9)/(3(p²−5)).  "
            "Proof: substitute n=p²+1 into B_wick=n(13n−10)/(6(n−1)) and "
            "η_part from 15.215 H; clear denominators."
        ),
        "formula": "4(p^4-8p^2-9)/(3(p^2-5))",
        "by_p_sample": sample,
        "depends_on": ["prop_15_215_eta_part", "prop_15_198_B_wick"],
    }


def theorem_path_C_inside_wick_gap() -> dict:
    ok = True
    sample = {}
    for p in range(5, 301):
        if not is_prime(p):
            continue
        d = hyp_inside_wick_gap(p)
        alt = gap_wick_closed(p) - room_hyp_over_24(p)
        if d != alt or d <= 0:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 23, 31, 47, 101, 199, 293):
            sample[str(p)] = {
                "gap_minus_hyp": str(d),
                "positive": d > 0,
            }
    return {
        "proved": ok,
        "theorem": (
            "For every prime p≥5: "
            "gap_wick − room_hyp/24 "
            "= 16 p² (p²−9)/[3(p²−5)(p²+1)] > 0.  "
            "Hence room_hyp/24 < gap_wick, so Path-C residual "
            "δ²≤room_hyp/24 implies ∑η² < B_wick ⇒ K₄≤Wick_hi "
            "⇒ residual-(i) dual-eq empty (15.198)."
        ),
        "formula": "16 p^2 (p^2-9)/(3(p^2-5)(p^2+1))",
        "by_p_sample": sample,
        "depends_on": [
            "theorem_gap_wick_closed_form",
            "prop_15_117_room_hyp",
            "prop_15_198_wick_blocks_dual_eq",
        ],
    }


def theorem_CS_closes_if_k2_over_N_le_gap() -> dict:
    return {
        "proved": True,
        "theorem": (
            "By 15.222 CS: ‖δ‖₂² ≤ k₂/|Max+|.  "
            "If k₂/|Max+| ≤ gap_wick then ∑η² ≤ B_wick ⇒ K₄≤Wick_hi "
            "⇒ residual (i).  "
            "If k₂/|Max+| ≤ R-room then ‖m₄‖₂² ≤ n(n−2)/4 ⇒ residual (i) "
            "via 15.217."
        ),
        "depends_on": ["prop_15_222_CS", "theorem_gap_wick_closed_form"],
    }


def theorem_delta_bound_open() -> dict:
    return {
        "proved": False,
        "theorem": (
            "OPEN: δ² ≤ room_hyp/24 (Path C), or k₂/|Max+| ≤ gap_wick "
            "(or ≤ R-room) Max+-free for all primes p≥5."
        ),
        "certified_small_p": {
            "p5": {
                "k2": 6912,
                "Nmax": 260,
                "k2_over_N": str(Fraction(6912, 260)),
                "gap_wick": str(gap_wick_closed(5)),
                "R_room": str(delta_room_for_R(5)),
                "hyp": str(room_hyp_over_24(5)),
                "CS_le_gap": Fraction(6912, 260) <= gap_wick_closed(5),
                "CS_le_R_room": Fraction(6912, 260) <= delta_room_for_R(5),
            },
            "p7": {
                "k2": 193705,
                "Nmax": 11452,
                "k2_over_N": str(Fraction(193705, 11452)),
                "gap_wick": str(gap_wick_closed(7)),
                "R_room": str(delta_room_for_R(7)),
                "CS_le_gap": Fraction(193705, 11452) <= gap_wick_closed(7),
                "CS_le_R_room": Fraction(193705, 11452) <= delta_room_for_R(7),
            },
        },
        "depends_on": [
            "theorem_path_C_inside_wick_gap",
            "theorem_CS_closes_if_k2_over_N_le_gap",
        ],
    }


def residual_i_closed_via_223() -> bool:
    return False


def hinge_status_223() -> dict:
    return {
        "gap_wick_closed": theorem_gap_wick_closed_form()["proved"],
        "path_C_inside_wick": theorem_path_C_inside_wick_gap()["proved"],
        "CS_conditional": theorem_CS_closes_if_k2_over_N_le_gap()["proved"],
        "delta_bound": theorem_delta_bound_open()["proved"],
        "residual_i_closed_via_223": residual_i_closed_via_223(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "δ²≤room_hyp/24, or k₂/|Max+|≤gap_wick, or |μ|≤1/(2p), "
            "or free-e+ker=sc, or K₄≤Wick by other means."
        ),
    }


def main() -> dict:
    a = theorem_gap_wick_closed_form()
    b = theorem_path_C_inside_wick_gap()
    c = theorem_CS_closes_if_k2_over_N_le_gap()
    d = theorem_delta_bound_open()
    status = hinge_status_223()
    out = {
        "prop": "15.223",
        "title": (
            "Path-C residual ⇒ Wick residual-(i); "
            "gap_wick and hyp gap closed forms; δ bound OPEN"
        ),
        "proved": {
            "gap_wick_closed": a["proved"],
            "path_C_inside_wick": b["proved"],
            "CS_conditional": c["proved"],
            "delta_bound": d["proved"],
            "residual_i_closed": residual_i_closed_via_223(),
        },
        "algebra": {
            "gap_wick_p5": str(gap_wick_closed(5)),
            "hyp_inside_p5": str(hyp_inside_wick_gap(5)),
            "by_p_gap": {
                str(p): str(gap_wick_closed(p)) for p in (5, 7, 11, 13, 17)
            },
            "delta_path": d,
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_223(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Path C ⇒ Wick ⇒ residual (i) is Max+-free once Path C holds.  "
            "δ² bound still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15223.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.223  residual_i={residual_i_closed_via_223()}")
    print(f"  gap_wick={a['proved']} pathC_in_wick={b['proved']}")
    print(f"  gap_wick(p=5)={gap_wick_closed(5)}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
