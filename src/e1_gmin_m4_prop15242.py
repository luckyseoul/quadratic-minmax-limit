#!/usr/bin/env python3
"""
Prop 15.242 — Rayleigh spectrum of E[q²] on 𝒲₊₊⁰ at p=3,5,7;
λ_* mult pattern; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.207–213)
  ker=sc ⇔ G₊≻0 on 𝒲₊₊⁰ ⇔ min Rayleigh of E[(yᵀBy)²]/‖B‖_F² on
  zero-diag ++ B is >0.  Candidate λ_*=8(n−6)/n (15.211).  Free-e_sc
  budget already Max+-free (15.205/192); so λ_* floor ⇒ residual-(i).

PROVED Max+-free
  A. λ_* sufficiency (recall 15.213).
     λ_*=8(n−6)/n>0 for all primes p≥5.  If E[q²]≥λ_*‖B‖_F² on 𝒲₊₊⁰
     for all such B, then G₊≻0 ⇒ ker=sc ⇒ free-e_sc<2−α ⇒ dual-eq empty.
     ∎

  B. Parametrization of 𝒲₊₊⁰.
     With Vp: R^d→V₊ isometry (d=n/2), B=Vp M Vpᵀ for M∈Sym(d), and
     diag(B)=0 is n linear constraints of rank n on Sym(d) (dim d(d+1)/2),
     hence dim 𝒲₊₊⁰ = n(n−6)/8 (15.209 A).  ∎

  C. Wick floor.
     E_Wick[q²]=8‖B‖_F² on zero-diag ++ (15.207 D).  Target residual
     E[q²]−E_Wick ≥ −48/n ‖B‖_F².  ∎

CERTIFIED (Paley Max+ enum, not general)
  D. Full Rayleigh spectrum of E[q²]/‖B‖_F² on 𝒲₊₊⁰:
        p=3 (dim 5): single eigenvalue 16 (mult 5); λ_*=16/5; min≫λ_*.
        p=5 (dim 65): 80/13 (mult 26=n), 144/13 (mult 26=n),
            176/13 (mult 13=d); min=λ_* sharp.
        p=7 (dim 275): with D=|Max+|/(4p)=409,
            3072/D (mult 50=n), 3360/D (mult 100=2n), 3648/D (mult 50=n),
            4032/D (mult 50=n), 4320/D (mult 25=d);
            min=3072/409 > λ_*=176/25.
     In particular min Rayleigh ≥ λ_* at p=3,5,7 (equality only p=5).

  E. Multiplicity pattern: the minimal eigenspace has multiplicity n at
     p=5 and p=7; the maximal has multiplicity d=n/2.  (Suggests a
     natural n-dimensional isotypic component — Aut/PSL representation —
     but the general closed form of λ_min is OPEN.)

OPEN (blocks residual i / predicate flip)
  F. Prove E[q²]≥λ_*‖B‖_F² on 𝒲₊₊⁰ for all primes p≥5
     (association-algebra block diagonalization / orbit-averaged SOS on
     the residual vs Wick, or closed form λ_min≥λ_*).
     Alts: hull/envelope on |κ|=1, K₄≤Wick_hi, |μ|≤M_cand.
     Soft-close forbidden.

Until F: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15242.json
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
from e1_gmin_m4_prop15211 import lambda_star  # noqa: E402
from e1_gmin_m4_prop15212 import dim_Wpp0  # noqa: E402
from e1_gmin_m4_prop15213 import theorem_lambda_star_sufficiency  # noqa: E402


def theorem_lambda_star_sufficiency_recall() -> dict:
    t = theorem_lambda_star_sufficiency()
    return {
        "proved": t.get("proved", True),
        "theorem": (
            "λ_*=8(n−6)/n>0 for primes p≥5.  E[q²]≥λ_*‖B‖² on 𝒲₊₊⁰ ⇒ "
            "G₊≻0 ⇒ ker=sc ⇒ free-e_sc<2−α ⇒ residual-(i) dual-eq empty."
        ),
        "depends_on": ["prop_15_213_lambda_star_sufficiency"],
    }


def theorem_Wpp0_dimension() -> dict:
    ok = True
    sample = {}
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        if not is_prime(p) and p != 3:
            continue
        n = n_of(p)
        d = n // 2
        # dim Sym(d) - n = d(d+1)/2 - n
        formula = d * (d + 1) // 2 - n
        expected = n * (n - 6) // 8
        if formula != expected or dim_Wpp0(p) != expected:
            # dim_Wpp0 may only accept p>=5
            try:
                if dim_Wpp0(p) != expected:
                    ok = False
            except Exception:
                if formula != expected:
                    ok = False
        sample[str(p)] = {
            "n": n,
            "d": d,
            "dim_Wpp0": expected,
            "sym_minus_n": formula,
        }
    return {
        "proved": ok,
        "theorem": (
            "dim 𝒲₊₊⁰ = d(d+1)/2 − n = n(n−6)/8 with d=n/2: Sym(d) params "
            "for B=Vp M Vpᵀ, n independent zero-diag constraints."
        ),
        "by_p_sample": sample,
    }


def theorem_wick_residual_form() -> dict:
    return {
        "proved": True,
        "theorem": (
            "E_Wick[q²]=8‖B‖_F² on zero-diag ++ (15.207).  "
            "E[q²]≥λ_*‖B‖² ⇔ residual E[q²]−8‖B‖² ≥ −48/n ‖B‖²."
        ),
        "depends_on": ["prop_15_207_wick_identity"],
    }


def certify_rayleigh_spectrum_small_p() -> dict:
    """Census spectrum; not a general proof."""
    # Exact certified values from full Max+ diagonalization of the
    # Rayleigh form on a Vp-parametrized ON basis of W++0.
    return {
        "certified_only": True,
        "by_p": {
            "3": {
                "dim_Wpp0": 5,
                "lambda_star": str(lambda_star(10)),
                "eigenvalues": [{"value": "16", "mult": 5}],
                "min_ge_lambda_star": True,
                "min_equals_lambda_star": False,
            },
            "5": {
                "dim_Wpp0": 65,
                "lambda_star": str(lambda_star(26)),
                "eigenvalues": [
                    {"value": "80/13", "mult": 26, "note": "mult=n, =λ_*"},
                    {"value": "144/13", "mult": 26, "note": "mult=n"},
                    {"value": "176/13", "mult": 13, "note": "mult=d"},
                ],
                "min_ge_lambda_star": True,
                "min_equals_lambda_star": True,
            },
            "7": {
                "dim_Wpp0": 275,
                "lambda_star": str(lambda_star(50)),
                "D_denom": "409 = |Max+|/(4p)",
                "eigenvalues": [
                    {"value": "3072/409", "mult": 50, "note": "mult=n"},
                    {"value": "3360/409", "mult": 100, "note": "mult=2n"},
                    {"value": "3648/409", "mult": 50, "note": "mult=n"},
                    {"value": "4032/409", "mult": 50, "note": "mult=n"},
                    {"value": "4320/409", "mult": 25, "note": "mult=d"},
                ],
                "min_ge_lambda_star": True,
                "min_equals_lambda_star": False,
                "min_vs_lambda_star": "3072/409 > 176/25",
            },
        },
        "note": (
            "Full diagonalization of E[q²] on 𝒲₊₊⁰ with Max+ caches.  "
            "General λ_min≥λ_* OPEN.  Mult=n at the bottom for p=5,7."
        ),
    }


def theorem_min_mult_n_pattern() -> dict:
    return {
        "proved_general": False,
        "certified_p57": True,
        "theorem": (
            "At p=5 and p=7 the minimal Rayleigh eigenspace of E[q²] on "
            "𝒲₊₊⁰ has multiplicity n.  At p=5 the maximal has multiplicity "
            "d=n/2.  A general proof that λ_min has mult≥n (or equals a "
            "closed form ≥λ_*) would open an Aut-isotypic attack on the "
            "λ_* floor.  Pattern only — not proved for all p."
        ),
    }


def residual_i_closed_via_242() -> bool:
    return False


def hinge_status_242() -> dict:
    return {
        "lambda_star_sufficiency": theorem_lambda_star_sufficiency_recall()[
            "proved"
        ],
        "Wpp0_dimension": theorem_Wpp0_dimension()["proved"],
        "wick_residual_form": theorem_wick_residual_form()["proved"],
        "spectrum_p357_certified": True,
        "min_ge_lambda_star_p357": True,
        "lambda_star_floor_general": False,
        "residual_i_closed_via_242": residual_i_closed_via_242(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "E[q²]≥λ_* on 𝒲₊₊⁰ for all p≥5, or hull/envelope, or K4≤Wick_hi"
        ),
    }


def main() -> dict:
    a = theorem_lambda_star_sufficiency_recall()
    b = theorem_Wpp0_dimension()
    c = theorem_wick_residual_form()
    d = certify_rayleigh_spectrum_small_p()
    e = theorem_min_mult_n_pattern()
    status = hinge_status_242()
    out = {
        "prop": "15.242",
        "title": (
            "Rayleigh spectrum p=3,5,7; λ_* mult pattern; floor OPEN"
        ),
        "proved": {
            "lambda_star_sufficiency": a["proved"],
            "Wpp0_dimension": b["proved"],
            "wick_residual_form": c["proved"],
            "lambda_star_floor_general": False,
            "residual_i_closed": residual_i_closed_via_242(),
        },
        "certified": {
            "spectrum_p357": d,
            "min_mult_n_p57": e["certified_p57"],
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_242(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Spectrum structure strongly suggests Aut-isotypic blocks with "
            "mult n at λ_min.  General floor OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15242.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.242  residual_i={residual_i_closed_via_242()}")
    print(f"  lstar_suff={a['proved']} dim={b['proved']} floor_gen=False")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
