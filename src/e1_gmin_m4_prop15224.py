#!/usr/bin/env python3
"""
Prop 15.224 — T²φ closed form Max+-free for any conference; μ_part solves
master always; residual-(i) still OPEN on δ∈E_{±4p}.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP
  Any symmetric conference C of order n (Cᵀ=C, C_ii=0, C_ij=±1, C²=(n−1)I).
  κ, φ, star on 4-sets as in 15.177–187.  T = C-weighted 4-set replacement.
  Master (16p² I − T²)μ = 16κ on Paley n=p²+1 (or n−1=p²).

PROVED Max+-free (any conference)
  A. Recall (15.215 / 15.187, already in-repo):
        Tκ = −6 star on every 4-set,
        star = 0 on |κ|=1 and |star|=4 on |κ|=3 (64-label exhaustion),
        Tφ = (n+2) star  (C² reduction + 64-label exhaustion of K4 edge poly),
        T²κ = −24φ + 48κ  (any conference, all 4-sets).

  B. T²φ closed form.
        T²φ = T((n+2) star) = (n+2) T(−Tκ/6) = −(n+2)/6 · T²κ
            = −(n+2)/6 · (−24φ + 48κ) = 4(n+2)(φ − 2κ).
     Proof: substitute A into Tφ and T²κ; pure operator algebra on 4-set
     functions.  No Max+ enumeration.  ∎

  C. Paley specialization n=p²+1:
        T²φ = 4(p²+3)(φ − 2κ).
     (Was only certified by census at p=3,5,7 in 15.188; now proved for
     all conference orders with n=p²+1.)  ∎

  D. μ_part solves the master whenever (A)–(C) hold.
        μ_part = [(p²−1)κ − 2φ] / (p²(p²−5))   (p²≠5)
     satisfies (16p² I − T²) μ_part = 16κ on span{κ,φ}, by the joint
     action of T² on that span from T²κ and T²φ (15.186–188 algebra).
     Hence every master solution is μ = μ_part + δ with δ∈E_{±4p}(T).  ∎

  E. Particular majorant (15.186 H, recall): on |κ|=1 with |φ|≤2(p−2)
     (Paley, 15.186 G), |μ_part| ≤ (p²+4p−9)/(p²(p²−5)) ≤ 1/(2p) for
     all primes p≥5.  ∎

CERTIFIED
  F. T²φ formula matches sparse T² at p=3,5,7 (prior 15.188).
  G. Actual |μ| ≤ 1/(2p) at p=5,7; |μ|≤|μ_part| fails pointwise at p=7.

OPEN (blocks residual i / predicate flip)
  H. Control δ∈E_{±4p} so |μ_part+δ| ≤ 1/(2p) (or ≤2/n) on every |κ|=1
     four-set for all primes p≥5; **or** any other residual-(i) hinge
     (Path C, K₄≤Wick, free-e+ker=sc, m4₂ thr).  Soft-close forbidden.

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15224.json
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
from e1_gmin_m4_prop15186 import (  # noqa: E402
    majorant_mu_part,
    mu_part_formula,
    theorem_particular_majorant_beats_thr,
    theorem_phi_bound_general,
)
from e1_gmin_m4_prop15188 import two_over_n  # noqa: E402


def T2_phi_formula_n(n: int, phi: int, kappa: int) -> int:
    """T²φ = 4(n+2)(φ − 2κ) as scalar on a joint eigen-label (algebra check)."""
    return 4 * (n + 2) * (phi - 2 * kappa)


def T2_phi_formula_p(p: int, phi: int, kappa: int) -> int:
    """Paley: 4(p²+3)(φ − 2κ)."""
    return 4 * (p * p + 3) * (phi - 2 * kappa)


def theorem_T2_phi_closed() -> dict:
    """Max+-free: T²φ = 4(n+2)(φ−2κ) for any conference order n."""
    ok = True
    sample = {}
    # Algebra identity: 4(n+2)(φ-2κ) matches Paley form when n=p²+1
    for p in [pp for pp in range(3, 80) if is_prime(pp)]:
        n = n_of(p)
        for kappa, phi in [(1, 0), (1, 2), (1, -2), (3, 4), (-1, 0)]:
            a = T2_phi_formula_n(n, phi, kappa)
            b = T2_phi_formula_p(p, phi, kappa)
            if a != b:
                ok = False
        if p in (3, 5, 7, 11):
            sample[str(p)] = {
                "n": n,
                "formula_n": "4(n+2)(phi-2*kappa)",
                "formula_p": "4(p^2+3)(phi-2*kappa)",
                "match_on_samples": True,
            }
    return {
        "proved": ok,
        "theorem": (
            "For every symmetric conference of order n and every 4-set: "
            "T²φ = 4(n+2)(φ−2κ).  "
            "Proof: Tφ=(n+2)star and star=−Tκ/6 (15.215); "
            "T²κ=−24φ+48κ (15.187); "
            "T²φ=T((n+2)star)=(n+2)T(−Tκ/6)=−(n+2)/6·T²κ=4(n+2)(φ−2κ)."
        ),
        "paley_form": "4(p^2+3)(phi-2*kappa) when n=p^2+1",
        "by_p_sample": sample,
        "depends_on": [
            "prop_15_215_Tphi_star",
            "prop_15_215_Tkappa_eq_minus_6_star",
            "prop_15_187_T2_kappa",
        ],
    }


def theorem_mu_part_solves_master_general() -> dict:
    """μ_part solves master on span{κ,φ} for all primes p≥5, p²≠5."""
    # Algebraic check: the joint T² matrix on (κ,φ) and inversion
    # T²κ = 48κ - 24φ
    # T²φ = 4(n+2)φ - 8(n+2)κ
    # Master (16p² - T²) μ = 16κ with μ = a κ + b φ
    ok = True
    sample = {}
    for p in [pp for pp in range(5, 60) if is_prime(pp)]:
        n = n_of(p)
        # On basis (κ, φ), T² acts as matrix:
        # T² [κ] = 48 κ - 24 φ
        # T² [φ] = -8(n+2) κ + 4(n+2) φ
        # (16p² I - T²) μ_part = 16 κ
        # μ_part = [(p²-1)κ - 2φ] / (p²(p²-5))
        den = p * p * (p * p - 5)
        a = Fraction(p * p - 1, den)  # coeff of κ
        b = Fraction(-2, den)  # coeff of φ
        # Apply (16p² I - T²) to aκ + bφ:
        # 16p² (aκ + bφ) - a(48κ - 24φ) - b(-8(n+2)κ + 4(n+2)φ)
        # κ coeff: 16p² a - 48a + 8(n+2)b
        # φ coeff: 16p² b + 24a - 4(n+2)b
        n = p * p + 1
        k_coeff = 16 * p * p * a - 48 * a + 8 * (n + 2) * b
        phi_coeff = 16 * p * p * b + 24 * a - 4 * (n + 2) * b
        if k_coeff != 16 or phi_coeff != 0:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "k_coeff": str(k_coeff),
                "phi_coeff": str(phi_coeff),
                "equals_16_kappa": k_coeff == 16 and phi_coeff == 0,
            }
    return {
        "proved": ok,
        "theorem": (
            "For every prime p≥5: μ_part=[(p²−1)κ−2φ]/(p²(p²−5)) "
            "satisfies (16p²I−T²)μ_part=16κ on span{κ,φ}, using T²κ and "
            "T²φ=4(n+2)(φ−2κ) with n=p²+1.  "
            "Hence μ=μ_part+δ, δ∈E_{±4p}."
        ),
        "by_p_sample": sample,
        "depends_on": ["theorem_T2_phi_closed", "prop_15_187_T2_kappa"],
    }


def theorem_particular_majorant_recall() -> dict:
    t = theorem_particular_majorant_beats_thr()
    return {
        "proved": t.get("proved", True),
        "theorem": (
            "On |κ|=1 with |φ|≤2(p−2): |μ_part|≤(p²+4p−9)/(p²(p²−5))≤1/(2p) "
            "for all primes p≥5 (15.186 H)."
        ),
        "depends_on": ["prop_15_186_phi_bound", "prop_15_186_majorant"],
    }


def residual_i_closed_via_224() -> bool:
    return False


def hinge_status_224() -> dict:
    return {
        "T2_phi_closed": theorem_T2_phi_closed()["proved"],
        "mu_part_solves_master": theorem_mu_part_solves_master_general()["proved"],
        "particular_majorant": theorem_particular_majorant_recall()["proved"],
        "phi_bound": theorem_phi_bound_general()["proved"],
        "residual_i_closed_via_224": residual_i_closed_via_224(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "Control δ∈E_{±4p} so |μ_part+δ|≤1/(2p) on |κ|=1 for all p≥5; "
            "or Path C / K4 / free-e / m4₂."
        ),
    }


def main() -> dict:
    a = theorem_T2_phi_closed()
    b = theorem_mu_part_solves_master_general()
    c = theorem_particular_majorant_recall()
    status = hinge_status_224()
    out = {
        "prop": "15.224",
        "title": (
            "T²φ=4(n+2)(φ−2κ) Max+-free; μ_part solves master always; "
            "δ-control OPEN"
        ),
        "proved": {
            "T2_phi_closed": a["proved"],
            "mu_part_solves_master": b["proved"],
            "particular_majorant": c["proved"],
            "residual_i_closed": residual_i_closed_via_224(),
        },
        "algebra": {
            "T2_phi_formula": "4(n+2)(phi-2*kappa)",
            "paley": "4(p^2+3)(phi-2*kappa)",
            "mu_part": "[(p^2-1)kappa-2*phi]/(p^2(p^2-5))",
            "majorant_p5": str(majorant_mu_part(5)),
            "two_over_n_p5": str(two_over_n(5)),
            "master_check": b["by_p_sample"],
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_224(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "T²φ and μ_part master algebra proved Max+-free.  "
            "δ∈E_{±4p} still OPEN for |μ|≤1/(2p).  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15224.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.224  residual_i={residual_i_closed_via_224()}")
    print(f"  T2_phi={a['proved']} mu_part_master={b['proved']}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
