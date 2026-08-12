#!/usr/bin/env python3
"""
Prop 15.227 — Free-e dual sharp at μ_*: cost≡need; one-sided residual-(i)
reformulation; dual gap at −1/p; Aut-SOS hinge restated; residual OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.176, 15.218, 15.226)
  Row+ones free-e dual: W=a·row_e+b·1, 2a+b=1, a=min(½, 1/(2−G_min)),
  cost=α·sum_ne, need=2−α, α=6/(p n).
  Farkas thr μ_*=(6−4p)/(p(3p−2)); sharp |μ| thr μ_f=−μ_*/2 (15.226).
  Gsum_disj=2 m₁ μ₄ ⇒ G_min=−2 max_{|κ|=1}|μ₄|.

PROVED Max+-free
  A. Free-e dual is exactly calibrated to Farkas.
     Set G_min=μ_*.  Then a∈(0,½) and
        a_* = 1/(2−μ_*) = p(3p−2)/(6(p²−1)).
     Hence
        sum_ne = (n(n−1)−2)/2 + a_*·n·(2−n)
               = p(p²+1)/3 − 1,
        cost   = α·sum_ne = (2p³+2p−6)/(p³+p) = 2−α = need.
     Proof: substitute μ_*=(6−4p)/(p(3p−2)), n=p²+1, α=6/(p n) and clear
     positive denominators (Fraction identity for all primes p≥3 with
     p²≠1).  ∎
     Interpretation: the row+ones dual family has **zero slack** at the
     Farkas boundary.  cost<need requires G_min>μ_* (strict), equivalently
        max_{|κ|=1}|μ₄| < μ_f = (2p−3)/(p(3p−2)).
     (The convenient thr |μ|≤1/(2p) is a strict sublevel: −1/p>μ_*.)

  B. Dual gap at the convenient thr G_min=−1/p.
     With a=p/(2p+1),
        need − cost = (p−4)/(2p+1).
     For all primes p≥5: gap>0, so cost<need.  (Value 1/11 at p=5;
     increases to ½.)  Recovers 15.218 D with closed gap.  ∎

  C. One-sided residual-(i) free-e form.
     Dual-eq requires some κ∈ker(Gsum)∩box with κ_e=2−α.
     Free-e max κ_e ≤ cost(G_min).  Combined with (A)(B):
        G_min ≥ −1/p  ⇒  free-e max ≤ cost < 2−α  ⇒ dual-eq empty.
     Equivalently (15.176/218): max_{|κ|=1}|μ₄| ≤ 1/(2p) ⇒ G_min≥−1/p
     ⇒ residual-(i) dual-eq empty.  ∎
     One-sided sharpening: only the **negative** matching direction of
     Gsum binds Farkas (positive Gsum loosens the dual).  Writing
        μ_↓ := max{ −m₁ μ₄ : |κ|=1, m₁=C_e C_f for a disj pair of S},
     one has G_min = −2 μ_↓, and residual-(i) free-e/Farkas fires when
        μ_↓ ≤ 1/(2p)   (or the sharp μ_↓ < μ_f).
     On each |κ|=1 four-set the three matching signs give
     {m₁ μ₄} = {μ₄, μ₄, −μ₄} up to global sign, so μ_↓ = max|μ₄| on that
     set; globally μ_↓ = max_{|κ|=1}|μ₄|.  The “one-sided” language
     emphasises that the dual/Farkas obstruction is a lower bound on Gsum,
     not a two-sided Loewner constraint beyond PSD.  ∎

  D. Aut-SOS hinge (structure for E_{±4p}^G).
     By 15.225–226: μ=μ_part+δ with δ∈E_{±4p}^G, |μ_part|≤maj≤1/(2p)<μ_f,
     room_δ^f=μ_f−maj>room_δ=1/(2p)−maj>0, and ⟨δ,κ⟩=⟨δ,φ⟩=0.
     Residual-(i) closes if either
        (i)  |δ_S| ≤ room_δ  (or room_δ^f) for every |κ_S|=1, or
        (ii) max|μ_part+δ| ≤ 1/(2p) (or <μ_f) on |κ|=1,
     which is a finite-dimensional inequality on E_{±4p}^G (orbit coordinates).
     A parameter-dependent PSD/Schur certificate on that space (Gram of orbit
     indicators under T and the linear constraints ⟨δ,κ⟩=⟨δ,φ⟩=0) would
     prove (i) or (ii) uniformly in p.  No dim≤1 assumption (false for G
     at p=5,7).  ∎

  E. Quantitative dual gap monotonicity.
     For G∈[μ_*, 0], a(G)=1/(2−G) is decreasing in G, sum_ne is decreasing
     in a (since n(2−n)<0), hence cost is decreasing in G.  Therefore
        cost(G) ≤ cost(μ_*) = need    for G≥μ_*,
        cost(G) < need                 for G>μ_*,
     with equality only at G=μ_*.  ∎

CERTIFIED
  F. Identities A,B as exact Fraction for all primes 3≤p≤300.
  G. Census: max|μ|<μ_f and ≤1/(2p) at p=5,7; free-e dual cost matches.

OPEN (blocks residual i / predicate flip)
  H. Prove max_{|κ|=1}|μ₄| ≤ 1/(2p) (or <μ_f / ≤M_cand / ≤2/n) for all
     primes p≥5 — e.g. by Aut-SOS/Schur on E_{±4p}^G controlling δ —
     **or** Path C / K₄ / free-e+ker=sc by another dual.  Soft-close forbidden.

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15227.json
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
from e1_gmin_m4_prop15190 import alpha_particular  # noqa: E402
from e1_gmin_m4_prop15218 import sum_ne_row_ones_dual  # noqa: E402
from e1_gmin_m4_prop15226 import (  # noqa: E402
    mu_farkas_threshold,
    room_delta_farkas,
)


def a_star_at_mu_star(p: int) -> Fraction:
    """a_* = p(3p−2)/(6(p²−1)) = 1/(2−μ_*)."""
    return Fraction(p * (3 * p - 2), 6 * (p * p - 1))


def sum_ne_at_mu_star(p: int) -> Fraction:
    """sum_ne = p(p²+1)/3 − 1 at G_min=μ_*."""
    return Fraction(p * (p * p + 1), 3) - 1


def cost_at_mu_star(p: int) -> Fraction:
    """cost = need = (2p³+2p−6)/(p³+p)."""
    return Fraction(2 * p**3 + 2 * p - 6, p**3 + p)


def dual_gap_at_minus_1_over_p(p: int) -> Fraction:
    """need − cost at G_min=−1/p equals (p−4)/(2p+1)."""
    return Fraction(p - 4, 2 * p + 1)


def theorem_dual_sharp_at_mu_star() -> dict:
    ok = True
    sample = {}
    for p in range(3, 301):
        if not is_prime(p):
            continue
        G = mu_star_threshold(p)
        a_from_G = Fraction(1, 2 - G)
        a_closed = a_star_at_mu_star(p)
        sne = sum_ne_row_ones_dual(p, G)
        sne_closed = sum_ne_at_mu_star(p)
        al = alpha_particular(p)
        cost = al * sne
        need = 2 - al
        cost_closed = cost_at_mu_star(p)
        if not (
            a_from_G == a_closed
            and sne == sne_closed
            and cost == need
            and cost == cost_closed
        ):
            ok = False
            break
        if p in (3, 5, 7, 11, 13, 17, 19, 23, 31, 47, 101):
            sample[str(p)] = {
                "a_star": str(a_closed),
                "sum_ne": str(sne_closed),
                "cost": str(cost_closed),
                "need": str(need),
                "cost_eq_need": cost == need,
            }
    return {
        "proved": ok,
        "theorem": (
            "For every prime p≥3 with p²≠1: at G_min=μ_* the row+ones dual "
            "has a=p(3p−2)/(6(p²−1)), sum_ne=p(p²+1)/3−1, and "
            "cost=need=(2p³+2p−6)/(p³+p).  The dual family is sharp at Farkas."
        ),
        "by_p_sample": sample,
        "formulas": {
            "a_star": "p(3p-2)/(6(p^2-1))",
            "sum_ne": "p(p^2+1)/3 - 1",
            "cost_eq_need": "(2p^3+2p-6)/(p^3+p)",
        },
        "depends_on": ["prop_15_176_mu_star", "prop_15_218_row_ones_dual"],
    }


def theorem_dual_gap_minus_1_over_p() -> dict:
    ok = True
    sample = {}
    for p in range(5, 301):
        if not is_prime(p):
            continue
        G = Fraction(-1, p)
        al = alpha_particular(p)
        cost = al * sum_ne_row_ones_dual(p, G)
        need = 2 - al
        gap = need - cost
        gap_closed = dual_gap_at_minus_1_over_p(p)
        if gap != gap_closed or gap <= 0:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101):
            sample[str(p)] = {
                "gap": str(gap_closed),
                "cost": str(cost),
                "need": str(need),
            }
    return {
        "proved": ok,
        "theorem": (
            "For every prime p≥5: at G_min=−1/p the row+ones dual has "
            "need−cost=(p−4)/(2p+1)>0.  Hence cost<2−α ⇒ free-e max<2−α "
            "⇒ dual-eq empty whenever G_min≥−1/p."
        ),
        "by_p_sample": sample,
        "gap_formula": "(p-4)/(2p+1)",
        "depends_on": ["prop_15_218_row_ones_dual"],
    }


def theorem_onesided_residual_i_form() -> dict:
    return {
        "proved": True,
        "theorem": (
            "G_min=−2 max_{|κ|=1}|μ₄|.  Free-e/Farkas residual-(i) fires when "
            "G_min≥−1/p (⇔ max|μ₄|≤1/(2p)), which is a strict sublevel of the "
            "sharp dual threshold G_min>μ_* (⇔ max|μ₄|<μ_f).  The obstruction "
            "is one-sided (lower bound on Gsum); on each |κ|=1 four-set the "
            "matching signs force μ_↓=max|μ₄|."
        ),
        "depends_on": [
            "theorem_dual_sharp_at_mu_star",
            "theorem_dual_gap_minus_1_over_p",
            "prop_15_218_Gsum_eq_2_m1_mu4",
            "prop_15_226_mu_f",
        ],
        "sufficient_target": "max_{|κ|=1}|μ₄| ≤ 1/(2p)",
        "sharp_target": "max_{|κ|=1}|μ₄| < μ_f = (2p-3)/(p(3p-2))",
    }


def theorem_dual_gap_monotone() -> dict:
    """cost decreasing in G on [μ_*, 0] ⇒ cost<need for G>μ_*."""
    ok = True
    sample = {}
    for p in [pp for pp in range(5, 80) if is_prime(pp)]:
        mu_s = mu_star_threshold(p)
        # sample G values between μ_* and 0
        grid = [
            mu_s,
            (mu_s + Fraction(-1, p)) / 2,
            Fraction(-1, p),
            Fraction(-1, 2 * p),
            Fraction(0),
        ]
        costs = []
        al = alpha_particular(p)
        need = 2 - al
        prev = None
        mono = True
        for G in grid:
            if G < mu_s - Fraction(1, 10**9):
                continue
            c = al * sum_ne_row_ones_dual(p, G)
            costs.append((str(G), str(c)))
            if prev is not None and c > prev + Fraction(1, 10**12):
                # cost should be non-increasing as G increases
                mono = False
            prev = c
        if not mono:
            ok = False
            break
        # equality only at μ_*
        c_star = al * sum_ne_row_ones_dual(p, mu_s)
        c_thr = al * sum_ne_row_ones_dual(p, Fraction(-1, p))
        if c_star != need or c_thr >= need:
            ok = False
            break
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "cost_at_mu_star": str(c_star),
                "cost_at_minus_1_over_p": str(c_thr),
                "need": str(need),
                "strict_gap": str(need - c_thr),
            }
    return {
        "proved": ok,
        "theorem": (
            "On G∈[μ_*,0]: a=1/(2−G) decreases in G, sum_ne decreases in a "
            "(n(2−n)<0), hence cost decreases in G.  Therefore cost(G)<need "
            "for all G>μ_*, with equality only at G=μ_*."
        ),
        "by_p_sample": sample,
        "depends_on": ["theorem_dual_sharp_at_mu_star"],
    }


def theorem_Aut_SOS_hinge_restated() -> dict:
    return {
        "proved_structure": True,
        "bound_proved_general": False,
        "theorem": (
            "Residual-(i) free-e/Farkas closes for all primes p≥5 if a "
            "parameter-dependent PSD/Schur certificate on E_{±4p}^G proves "
            "max|μ_part+δ|≤1/(2p) (or |δ|≤room_δ / room_δ^f) on |κ|=1, "
            "using P±, ⟨δ,κ⟩=⟨δ,φ⟩=0, and Tδ=4p(δ_{+}−δ_{-}).  "
            "dim E_{4p}^G ∈ {2,7} at p=5,7 (not ≤1)."
        ),
        "depends_on": [
            "prop_15_225_projectors",
            "prop_15_226_delta_G_invariant",
            "prop_15_226_room_delta_farkas",
            "theorem_onesided_residual_i_form",
        ],
        "room_delta_farkas_p5": str(room_delta_farkas(5)),
        "room_delta_farkas_p7": str(room_delta_farkas(7)),
        "mu_f_p5": str(mu_farkas_threshold(5)),
        "mu_f_p7": str(mu_farkas_threshold(7)),
    }


def residual_i_closed_via_227() -> bool:
    return False


def hinge_status_227() -> dict:
    return {
        "dual_sharp_at_mu_star": theorem_dual_sharp_at_mu_star()["proved"],
        "dual_gap_minus_1_over_p": theorem_dual_gap_minus_1_over_p()["proved"],
        "onesided_form": theorem_onesided_residual_i_form()["proved"],
        "dual_gap_monotone": theorem_dual_gap_monotone()["proved"],
        "Aut_SOS_structure": theorem_Aut_SOS_hinge_restated()[
            "proved_structure"
        ],
        "Aut_SOS_bound": theorem_Aut_SOS_hinge_restated()[
            "bound_proved_general"
        ],
        "residual_i_closed_via_227": residual_i_closed_via_227(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "Aut-SOS/Schur on E_pm4p^G for max|μ|≤1/(2p); "
            "or Path C / K4 / free-e+ker=sc."
        ),
    }


def main() -> dict:
    a = theorem_dual_sharp_at_mu_star()
    b = theorem_dual_gap_minus_1_over_p()
    c = theorem_onesided_residual_i_form()
    d = theorem_dual_gap_monotone()
    e = theorem_Aut_SOS_hinge_restated()
    status = hinge_status_227()
    out = {
        "prop": "15.227",
        "title": (
            "Free-e dual sharp at μ_* (cost≡need); gap (p−4)/(2p+1) at −1/p; "
            "one-sided residual-(i); Aut-SOS hinge; OPEN"
        ),
        "proved": {
            "dual_sharp_at_mu_star": a["proved"],
            "dual_gap_minus_1_over_p": b["proved"],
            "onesided_residual_i_form": c["proved"],
            "dual_gap_monotone": d["proved"],
            "Aut_SOS_structure": e["proved_structure"],
            "Aut_SOS_bound": e["bound_proved_general"],
            "residual_i_closed": residual_i_closed_via_227(),
        },
        "algebra": {
            "a_star": "p(3p-2)/(6(p^2-1))",
            "sum_ne_at_mu_star": "p(p^2+1)/3 - 1",
            "cost_eq_need": "(2p^3+2p-6)/(p^3+p)",
            "gap_at_minus_1_over_p": "(p-4)/(2p+1)",
            "a_star_p5": str(a_star_at_mu_star(5)),
            "gap_p5": str(dual_gap_at_minus_1_over_p(5)),
            "gap_p7": str(dual_gap_at_minus_1_over_p(7)),
            "cost_at_mu_star_p5": str(cost_at_mu_star(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_227(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Row+ones dual is exactly sharp at Farkas (cost≡need at μ_*).  "
            "Convenient thr −1/p has closed gap (p−4)/(2p+1).  "
            "Aut-SOS bound on E^G still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15227.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.227  residual_i={residual_i_closed_via_227()}")
    print(
        f"  sharp={a['proved']} gap={b['proved']} "
        f"mono={d['proved']} Aut_SOS_bound={e['bound_proved_general']}"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
