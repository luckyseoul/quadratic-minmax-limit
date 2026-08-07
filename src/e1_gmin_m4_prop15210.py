#!/usr/bin/env python3
"""
Prop 15.210 — D(C) cost census; 8/p dual budget; free-e/sc implication;
residual (i) still OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.203–209)
  Dual-eq empty if free-e max κ_e < 2−α on ker∩box.  sc⊆ker; ker=sc iff
  G₊≻0 on 𝒲₊₊⁰ (15.207–209).  D(C) is a Max+-free constructive dual for
  free-e over sc (15.203): cost_D = α·sum_ne bounds free-e_sc.
  8/p < 2−α for p≥5 (15.204), so sum_ne ≤ (4/3)n ⇒ cost_D ≤ 8/p < 2−α.

PROVED Max+-free
  A. Dual budget comparison (Fraction, all primes p≥5):
        8/p < 2 − α = 2 − 6/(p n)
     already proved in 15.204 (poly p³−4p²−p+7 > 0).  Hence any sc-dual
     with cost ≤ 8/p blocks dual-eq whenever ker=sc.  ∎

  B. Implication (structural, Max+-free):
        (ker=sc for all p≥5) ∧ (cost_D(p) ≤ 8/p for all p≥5)
        ⇒ free-e max = free-e_sc ≤ cost_D ≤ 8/p < 2−α
        ⇒ dual-eq empty for all primes p≥5
        ⇒ residual-(i) dual-eq branch blocked.  ∎
     (Does **not** assert either hypothesis; both remain open in general.
      cost_D is the constructive D(C) of 15.203.)

  C. Asymptotic room.  α ∼ 6/p³, scheme max ∼ 3/p; free-e candidate
     r(p)·scheme → (7/5)·(3/p) = 21/(5p).  8/p is a coarser but simpler
     dual budget than r(p)·scheme.  ∎

CERTIFIED (not general)
  D. D(C) cost for primes p∈{5,7,11,13,17,19,23}:
        cost_D · p ∈ [6.56, 7.68], max at p=7 (≈7.684)
        cost_D ≤ 8/p and cost_D < 2−α at every listed prime
        sum_ne / n ∈ [1.09, 1.28]
     Blocks dual-eq **if** ker=sc (certified ker=sc at p=3,5,7 by 15.209).

  E. G₊≻0 on 𝒲₊₊⁰ at p=3,5,7 (15.209); free-e = r·scheme there;
     cost_D is strictly looser than free-e max but still < 2−α.

OPEN (blocks residual i / predicate flip)
  F. Prove for all primes p≥5:
     (i)  cost_D(p) ≤ 8/p (or sum_ne(D(C)) ≤ (4/3)n, or free-e_sc ≤ r·scheme),
          **and** ker=sc (G₊≻0 on 𝒲₊₊⁰), **or**
     (ii) |μ₄| ≤ M_cand (or ≤1/(2p) / ≤2/n) on |κ|=1, **or**
     (iii) K₄ ≤ Wick_hi (or ≤16n²−14n).
     Then wire gsum_disj_lb_proved_general → residual_i → E1 → L only if
     residual (ii) full allows.  Soft-close forbidden.

Until F: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15210.json
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
)
from e1_gmin_m4_prop15190 import alpha_particular  # noqa: E402
from e1_gmin_m4_prop15203 import dual_construction_cost  # noqa: E402
from e1_gmin_m4_prop15204 import theorem_psd_ker_characterization  # noqa: E402


def theorem_eight_over_p_beats_need() -> dict:
    """8/p < 2−α for all primes p≥5 (recall 15.204)."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        al = alpha_particular(p)
        lhs = Fraction(8, p)
        need = 2 - al
        beats = lhs < need
        # poly from 15.204: p³ − 4p² − p + 7 > 0
        poly = p**3 - 4 * p * p - p + 7
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "8_over_p": str(lhs),
                "need": str(need),
                "beats": beats,
                "poly": poly,
            }
        if not beats or poly <= 0:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For all primes p≥5: 8/p < 2−α.  Hence any sc-dual with "
            "cost ≤ 8/p blocks dual-eq when ker=sc (15.204)."
        ),
        "by_p_sample": sample,
        "depends_on": ["prop_15_204"],
    }


def theorem_implication_costD_ker_sc() -> dict:
    return {
        "proved": True,
        "theorem": (
            "If ker(Gsum)=scheme⊕cross for all primes p≥5 and "
            "cost_D(p)≤8/p for all p≥5, then free-e max ≤ cost_D < 2−α, "
            "so dual-eq is empty and residual-(i) dual-eq branch is blocked.  "
            "Neither hypothesis is claimed here."
        ),
        "hypotheses_open": ["ker_eq_sc_general", "cost_D_le_8_over_p_general"],
        "depends_on": ["prop_15_203", "prop_15_204", "prop_15_207"],
    }


def certify_cost_D_budget(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7, 11, 13, 17, 19, 23]
    rows = {}
    all_le_8p = True
    all_lt_need = True
    max_cost_p = 0.0
    for p in primes:
        r = dual_construction_cost(p)
        al = float(alpha_particular(p))
        cost = float(r["cost"])
        need = float(r["need"])
        eight = 8.0 / p
        le8 = cost <= eight + 1e-9
        ltneed = cost < need - 1e-9
        if not le8:
            all_le_8p = False
        if not ltneed:
            all_lt_need = False
        max_cost_p = max(max_cost_p, cost * p)
        rows[str(p)] = {
            "sum_ne": r["sum_ne"],
            "sum_ne_over_n": r["sum_ne"] / n_of(p),
            "cost": cost,
            "cost_times_p": cost * p,
            "eight_over_p": eight,
            "need": need,
            "cost_le_8_over_p": le8,
            "cost_lt_need": ltneed,
            "t0": r["t0"],
            "blocks_if_ker_sc": ltneed,
        }
    return {
        "certified": all_le_8p and all_lt_need,
        "proved_general": False,
        "primes": primes,
        "by_p": rows,
        "max_cost_times_p": max_cost_p,
        "note": (
            f"cost_D·p ≤ {max_cost_p:.4f} < 8 on listed primes; "
            "cost_D ≤ 8/p and cost_D < 2−α throughout.  General bound open."
        ),
    }


def free_e_bound_proved_general() -> bool:
    return False


def cost_D_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_210() -> bool:
    return False


def hinge_status_210() -> dict:
    return {
        "eight_over_p_beats_need_proved": True,
        "implication_costD_ker_sc_proved": True,
        "cost_D_le_8_over_p_general": cost_D_bound_proved_general(),
        "ker_eq_sc_general": False,
        "residual_i_closed": residual_i_closed_via_210(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove cost_D≤8/p (or free-e_sc≤r·scheme) and ker=sc for all "
            "p≥5, or |μ|≤M_cand / K₄≤Wick_hi."
        ),
    }


def main() -> dict:
    a = theorem_eight_over_p_beats_need()
    b = theorem_implication_costD_ker_sc()
    # Reuse small prime set for speed in main; full census in certify when tested
    cert = certify_cost_D_budget([5, 7, 11])
    status = hinge_status_210()
    out = {
        "prop": "15.210",
        "title": "D(C) cost census; 8/p dual budget implication",
        "proved": {
            "eight_over_p_beats_need": a["proved"],
            "implication_costD_and_ker_sc": b["proved"],
        },
        "certified": {
            "cost_D_budget": cert,
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_210(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "soft_close": False,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15210.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {path}")
    print(f"  8/p < 2−α proved={a['proved']}")
    print(f"  implication proved={b['proved']}")
    print(f"  cost_D cert (p=5,7,11)={cert['certified']} max cost*p={cert['max_cost_times_p']:.4f}")
    print(f"  residual_i closed: {residual_i_closed_via_210()}")
    print(f"  gsum={gsum_disj_lb_proved_general()} e1={e1_closed_general()}")
    return out


if __name__ == "__main__":
    main()
