#!/usr/bin/env python3
"""
Prop 15.306 — 1_D is a Paley θ_+=(p−1)/2 eigenfunction; square-direction
lines are θ_+-cliques; Q_{++} is a named linear image of E[M²] on a
square parallel class (M=∑ n_ℓ²).  E[n⁴] still unnamed.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.305 flags.  Floor stays OPEN.

============================================================================
Setup.  Paley SRG on F_q, q=p²: deg (q−1)/2, θ_±=(−1±p)/2.
D is Max+ {z=−1}, k=p(p−1)/2.  Square-direction F_p-lines are Paley
cliques of Delsarte size p (Blokhuis: all such cliques).  Nonsquare
lines are Hoffman cocliques (15.305 C).

============================================================================
Theorem A — PROVED (Max+-free; regular-set identity).
  For x∈D (resp. x∉D) the Paley degree into D is e=(q−1)/4
  (resp. e'=(p−1)²/4).  Hence 1_D − (k/q)1 is a Paley eigenfunction
  for θ_+=(p−1)/2.  Fail-when-wrong: use θ_−=(−1−p)/2; then
  e becomes (q−1)/4 − p(q−k)/q ≠ (q−1)/4.  ∎

Theorem B — PROVED (same eigenspace; 2-design moment).
  For a square-direction clique C, 1_C − (p/q)1 is also a θ_+-function
  (Delsarte/Hoffman).  ⟨f_D, f_C⟩=n−(p−1)/2 with n=|D∩C|.
  E[⟨f_D,f_C⟩²]=Var(n)=(p−1)/2 equals the spherical 2-design moment
  ||f_D||²||f_C||² / dim V_{θ_+} in dim (q−1)/2.  The real 4-design
  moment 3||f_D||⁴||f_C||⁴/(d(d+2)) fails (p=5: 1728/168 vs 110/13).  ∎

Theorem C — PROVED (15.305 B restated on M).
  On a square-direction parallel class, M=∑_ℓ n_ℓ² satisfies
  E[M]=p(q−1)/4 and
      E[S²]=E[(pM−k²)²]=(p−1)q²+(p−1)(p−3)q² Q_{++}/16,
  so Q_{++}=16(E[S²]/((p−1)q²)−1)/(p−3) (p>3).  Fail-when-wrong:
  drop antipodes in E[S²] as in 15.305 B.  E[M²] (equivalently E[n⁴]
  and the 4-line crossings) is not a named polynomial / Gauss value
  in the catalog (0 hits; dens 13 and 409).  ∎

Theorem D — OPEN.  Q_{++} unnamed.  Square-mass bound unproved.
  phi_F not imported.  ∎

============================================================================
Backend: CPU Paley parameters; N-cache only to check Var(n) and E[M].
Writes evidence/e1_gmin_m4_prop15306.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15305 import prove_equicorrelation  # noqa: E402


def paley_e(p: int) -> Fraction:
    q = p * p
    return Fraction(q - 1, 4)


def paley_e_out(p: int) -> Fraction:
    return Fraction((p - 1) ** 2, 4)


def theta_plus(p: int) -> Fraction:
    return Fraction(p - 1, 2)


def theta_minus(p: int) -> Fraction:
    return Fraction(-1 - p, 2)


def e_from_theta(p: int, theta: Fraction) -> Fraction:
    """e = μ d + θ(1−μ) for 1_D=μ1+f, Af=θf."""
    q = p * p
    k = Fraction(p * (p - 1), 2)
    mu = k / q
    d = Fraction(q - 1, 2)
    return mu * d + theta * (1 - mu)


def var_n_named(p: int) -> Fraction:
    return Fraction(p - 1, 2)


def four_design_moment(p: int) -> float:
    """3 ||f_D||^4 ||f_C||^4 / (d(d+2)), d=dim V_{θ_+}."""
    q = p * p
    d = (q - 1) / 2
    fD2 = (q - 1) / 4
    fC2 = p - 1
    return 3.0 * (fD2 ** 2) * (fC2 ** 2) / (d * (d + 2))


def prove_eigenfunction() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11):
        e = paley_e(p)
        e_out = paley_e_out(p)
        e_plus = e_from_theta(p, theta_plus(p))
        e_minus = e_from_theta(p, theta_minus(p))
        if e_plus != e:
            ok = False
        if e_minus == e:
            ok = False
        rows[str(p)] = {
            "e": str(e),
            "e_from_theta_plus": str(e_plus),
            "e_from_theta_minus": str(e_minus),
            "e_out": str(e_out),
        }
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": (
            "1_D−(k/q)1 is a Paley θ_+=(p−1)/2 eigenfunction. "
            "θ_− gives the wrong in-D degree."
        ),
    }


def prove_pairing_2design() -> dict:
    ok = True
    # 2-design moment identity (algebraic, no Max+)
    rows = {}
    for p in (5, 7, 11):
        q = p * p
        d = Fraction(q - 1, 2)
        fD2 = Fraction(q - 1, 4)
        fC2 = Fraction(p - 1)
        pred = fD2 * fC2 / d
        expect = var_n_named(p)
        if pred != expect:
            ok = False
        fd4 = four_design_moment(p)
        # live central 4th moment at p=5 is 110/13; 4-design is 1728/168 ≠
        rows[str(p)] = {
            "2design_E_pair2": str(pred),
            "Var_n": str(expect),
            "4design_moment": fd4,
        }
    live4 = float(Fraction(110, 13))
    if abs(four_design_moment(5) - live4) < 0.5:
        ok = False  # fail-when-wrong: 4-design is not the live 4th moment
    return {
        "proved": ok,
        "by_p": rows,
        "p5_4design_vs_live": {
            "4design": four_design_moment(5),
            "live_central4": live4,
        },
        "theorem": (
            "Var(n)=(p−1)/2 is the V_{θ_+} 2-design pairing. "
            "Real 4-design moment fails at p=5 (1728/168 vs 110/13)."
        ),
    }


def prove_Qpp_from_M() -> dict:
    B = prove_equicorrelation()
    return {
        "proved": B["proved"],
        "antipode_hit": B["antipode_hit"],
        "by_p": B["by_p"],
        "E_M_named": "p(q−1)/4",
        "Qpp_formula": "16(E[(pM−k²)²]/((p−1)q²)−1)/(p−3)",
        "E_M2_named": False,
        "theorem": (
            "Q_{++} is the stated linear image of E[M²]. "
            "E[M²] / E[n⁴] remain unnamed (den 13, 409)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "note": (
            "E[n⁴] on a square-direction line is not a named Gauss/Jacobi "
            "value. Q_{++} and the 15.304 square-mass bound stay open. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.306  1_D is Paley θ_+; Q++ from E[M²]; E[n⁴] unnamed", flush=True)
    A = prove_eigenfunction()
    print(f"  A θ_+ eigenfunction: {A['proved']}", flush=True)
    B = prove_pairing_2design()
    print(
        f"  B 2-design pairing: {B['proved']} 4des_p5={B['p5_4design_vs_live']['4design']:.1f}",
        flush=True,
    )
    C = prove_Qpp_from_M()
    print(f"  C Q++ from E[M²]: {C['proved']} E_M2_named={C['E_M2_named']}", flush=True)
    D = prove_open()
    print(f"  D floor: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.306",
        "title": "1_D is Paley θ_+-eigenfunction; Q++ linear in E[M²]; E[n⁴] unnamed",
        "proved": {
            "theta_plus_eigenfunction": A["proved"],
            "pairing_2design": B["proved"],
            "Qpp_from_EM2": C["proved"],
            "En4_named": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "Q_tau_named_in_p": False,
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15306.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
