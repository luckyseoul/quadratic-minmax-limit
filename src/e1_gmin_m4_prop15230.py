#!/usr/bin/env python3
"""
Prop 15.230 — R_part max ≤ R̄₄-budget (Cy≡δ channel); residual-(i) OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.229, 15.225–226)
  On |κ|=1: (p⁴−1)μ + 2φ = R̄₄ (size-4 non-self Cy remainder).
  μ=μ_part+δ with μ_part=aκ+bφ, a=(p²−1)/den, b=−2/den, den=p²(p²−5).
  Define the particular Cy image
      R_part := (p⁴−1)μ_part + 2φ
               = (p⁴−1)a · κ + ((p⁴−1)b+2) · φ.
  Then R̄₄ = R_part + (p⁴−1)δ.  Budget for |R̄₄| sufficient for |μ|≤1/(2p):
      B(p) := (p⁴−1)/(2p) − 4(p−2)
  (15.229 G, under |φ|≤2(p−2)).

PROVED Max+-free
  A. Coefficients of R_part.
        (p⁴−1)a = (p⁴−1)(p²−1)/(p²(p²−5)) > 0,
        (p⁴−1)b+2 = 2(1−5p²)/(p²(p²−5)) < 0   for p≥5.  ∎

  B. Worst-case particular majorant of |R_part| on |κ|=1.
     On |κ|=1 with |φ|≤2(p−2):
        |R_part| ≤ R_part_max(p)
          := (p⁴−1)a + 2(p−2)·|(p⁴−1)b+2|
           = (p⁶ − p⁴ + 20p³ − 41p² − 4p + 9) / (p²(p²−5)).
     Proof: expand using (A) and |c_κ|+|c_φ|·2(p−2).  ∎

  C. R_part_max ≤ B for all primes p≥5.
        B − R_part_max
          = (p−1)(p+1)(p²+1)(p³−2p²−13p+18) / (2p²(p²−5)).
     Numerator: (p−1)(p+1)(p²+1)>0 and cubic p³−2p²−13p+18 has value 28
     at p=5 with derivative 3p²−4p−13>0 on [5,∞).  Denominator >0 for p≥5.
     Hence B − R_part_max > 0.  ∎
     (Same cubic as the particular majorant |μ_part|≤1/(2p) of 15.186 H.)

  D. Channel unification.
     R̄₄ = R_part + (p⁴−1)δ, so the Cy residual R̄₄−R_part is exactly
     the master residual δ scaled by (p⁴−1).  Moreover
        (B − R_part_max)/(p⁴−1) = room_δ := 1/(2p) − maj
     as rational identities.  Consequently the sufficient pointwise criteria
        |R̄₄| ≤ B    and    |δ| ≤ room_δ on |κ|=1
     are equivalent when |R_part| attains its majorant (worst φ), and in
     general |R̄₄|≤B is the Cy-channel form of residual-(i) Farkas via 15.229.
     The particular solution alone is R-safe (C); only the δ-correction is open.  ∎

CERTIFIED
  E. Identities B,C as exact Fraction for all primes 5≤p≤300.
  F. Census: p=5 max|R̄₄|=124/5=24.8 < B=252/5=50.4; p=3 has
     R_part_max=B_fail (maj>thr), consistent with residual scope p≥5.

OPEN (blocks residual i / predicate flip)
  G. Prove |R̄₄|≤B (equivalently control δ so |μ_part+δ|≤1/(2p)) on every
     |κ|=1 four-set for all primes p≥5 — e.g. by character-sum evaluation of
     size-4 maps (sum before abs), Aut-SOS+diamond, Path C, or K₄/free-e.
     Soft-close forbidden.

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15230.json
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
    residual_i_dual_eq_empty_proved_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import majorant_mu_part  # noqa: E402
from e1_gmin_m4_prop15225 import room_delta  # noqa: E402
from e1_gmin_m4_prop15229 import R_bar_budget_for_mu_thr  # noqa: E402


def mu_part_coeffs(p: int) -> tuple[Fraction, Fraction]:
    den = p * p * (p * p - 5)
    return Fraction(p * p - 1, den), Fraction(-2, den)


def R_part_coeffs(p: int) -> tuple[Fraction, Fraction]:
    """R_part = c_κ κ + c_φ φ."""
    a, b = mu_part_coeffs(p)
    D = p**4 - 1
    return D * a, D * b + 2


def R_part_max(p: int) -> Fraction:
    """Max |R_part| on |κ|=1, |φ|≤2(p−2)."""
    ck, cph = R_part_coeffs(p)
    return abs(ck) + abs(cph) * 2 * (p - 2)


def R_part_max_closed(p: int) -> Fraction:
    """(p⁶ − p⁴ + 20p³ − 41p² − 4p + 9)/(p²(p²−5))."""
    return Fraction(
        p**6 - p**4 + 20 * p**3 - 41 * p * p - 4 * p + 9,
        p * p * (p * p - 5),
    )


def budget_minus_R_part_max(p: int) -> Fraction:
    """(p−1)(p+1)(p²+1)(p³−2p²−13p+18)/(2p²(p²−5))."""
    return Fraction(
        (p - 1)
        * (p + 1)
        * (p * p + 1)
        * (p**3 - 2 * p * p - 13 * p + 18),
        2 * p * p * (p * p - 5),
    )


def theorem_R_part_coeffs() -> dict:
    ok = True
    sample = {}
    for p in range(5, 120):
        if not is_prime(p):
            continue
        ck, cph = R_part_coeffs(p)
        if ck <= 0 or cph >= 0:
            ok = False
            break
        # cph = 2(1-5p²)/(p²(p²-5))
        cph_closed = Fraction(2 * (1 - 5 * p * p), p * p * (p * p - 5))
        if cph != cph_closed:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23):
            sample[str(p)] = {"c_kappa": str(ck), "c_phi": str(cph)}
    return {
        "proved": ok,
        "theorem": (
            "R_part=(p⁴−1)a κ+((p⁴−1)b+2)φ with a,b the μ_part coefficients; "
            "(p⁴−1)a>0 and (p⁴−1)b+2=2(1−5p²)/(p²(p²−5))<0 for primes p≥5."
        ),
        "by_p_sample": sample,
    }


def theorem_R_part_max_closed() -> dict:
    ok = True
    sample = {}
    for p in range(5, 201):
        if not is_prime(p):
            continue
        if R_part_max(p) != R_part_max_closed(p):
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31):
            sample[str(p)] = str(R_part_max_closed(p))
    return {
        "proved": ok,
        "theorem": (
            "On |κ|=1 with |φ|≤2(p−2): |R_part|≤R_part_max with "
            "R_part_max=(p⁶−p⁴+20p³−41p²−4p+9)/(p²(p²−5))."
        ),
        "by_p_sample": sample,
        "depends_on": ["theorem_R_part_coeffs", "prop_15_186_phi_bound"],
    }


def theorem_R_part_max_le_budget() -> dict:
    ok = True
    sample = {}
    for p in range(5, 301):
        if not is_prime(p):
            continue
        B = R_bar_budget_for_mu_thr(p)
        rm = R_part_max_closed(p)
        d = budget_minus_R_part_max(p)
        if B - rm != d or d <= 0:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101):
            sample[str(p)] = {
                "R_part_max": str(rm),
                "budget": str(B),
                "slack": str(d),
            }
    # cubic positivity
    cubic_ok = all(
        (p**3 - 2 * p * p - 13 * p + 18) > 0
        for p in range(5, 301)
        if is_prime(p)
    )
    return {
        "proved": ok and cubic_ok,
        "theorem": (
            "For every prime p≥5: R_part_max ≤ B(p) with "
            "B−R_part_max=(p−1)(p+1)(p²+1)(p³−2p²−13p+18)/(2p²(p²−5))>0.  "
            "The cubic is the same as in |μ_part|≤1/(2p) (value 28 at p=5)."
        ),
        "by_p_sample": sample,
        "depends_on": ["theorem_R_part_max_closed", "prop_15_229_R_budget"],
    }


def theorem_channel_unification() -> dict:
    ok = True
    sample = {}
    for p in range(5, 120):
        if not is_prime(p):
            continue
        # (B - R_part_max)/(p^4-1) == room_delta
        slack = budget_minus_R_part_max(p)
        equiv = slack / (p**4 - 1)
        rd = room_delta(p)
        if equiv != rd:
            ok = False
            break
        if p in (5, 7, 11, 13, 17):
            sample[str(p)] = {
                "room_delta": str(rd),
                "slack_over_p4m1": str(equiv),
                "match": True,
            }
    return {
        "proved": ok,
        "theorem": (
            "R̄₄=R_part+(p⁴−1)δ and (B−R_part_max)/(p⁴−1)=room_δ as rational "
            "identities.  Hence the Cy-channel budget test and the δ-room "
            "test are equivalent at worst-case φ; the particular solution is "
            "R-safe for all primes p≥5, and residual-(i) reduces to controlling δ "
            "(or the full R̄₄) on |κ|=1."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_R_part_max_le_budget",
            "prop_15_225_mu_part_eq_mu_mn",
            "prop_15_229_functional_eq",
        ],
        "sufficient_for_residual_i": False,
        "open": "|R̄₄|≤B or |δ|≤room_δ on |κ|=1 for all primes p≥5",
    }


def residual_i_closed_via_230() -> bool:
    return False


def hinge_status_230() -> dict:
    return {
        "R_part_coeffs": theorem_R_part_coeffs()["proved"],
        "R_part_max_closed": theorem_R_part_max_closed()["proved"],
        "R_part_max_le_budget": theorem_R_part_max_le_budget()["proved"],
        "channel_unification": theorem_channel_unification()["proved"],
        "residual_i_closed_via_230": residual_i_closed_via_230(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "|R̄₄|≤B (Cy) or |δ|≤room_δ on |κ|=1; same open scalar after "
            "unification.  Alts: Path C, K4, free-e+ker=sc, M_cand."
        ),
    }


def main() -> dict:
    a = theorem_R_part_coeffs()
    b = theorem_R_part_max_closed()
    c = theorem_R_part_max_le_budget()
    d = theorem_channel_unification()
    status = hinge_status_230()
    out = {
        "prop": "15.230",
        "title": (
            "R_part_max≤R̄₄-budget for p≥5; Cy≡δ channel unification; "
            "residual OPEN"
        ),
        "proved": {
            "R_part_coeffs": a["proved"],
            "R_part_max_closed": b["proved"],
            "R_part_max_le_budget": c["proved"],
            "channel_unification": d["proved"],
            "residual_i_closed": residual_i_closed_via_230(),
        },
        "algebra": {
            "R_part_max": "(p^6-p^4+20p^3-41p^2-4p+9)/(p^2(p^2-5))",
            "budget_minus_R_part_max": (
                "(p-1)(p+1)(p^2+1)(p^3-2p^2-13p+18)/(2p^2(p^2-5))"
            ),
            "R_bar": "R_part + (p^4-1)*delta",
            "R_part_max_p5": str(R_part_max_closed(5)),
            "budget_p5": str(R_bar_budget_for_mu_thr(5)),
            "slack_p5": str(budget_minus_R_part_max(5)),
            "room_delta_p5": str(room_delta(5)),
            "maj_p5": str(majorant_mu_part(5)),
            "thr_p5": str(mu4_threshold(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_230(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Particular Cy image R_part is inside the R̄₄ budget for all "
            "primes p≥5.  Residual-(i) is exactly the δ-correction to R̄₄ "
            "(same room as thr−maj).  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15230.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.230  residual_i={residual_i_closed_via_230()}")
    print(
        f"  R_part_max_le_B={c['proved']} unify={d['proved']} "
        f"slack_p5={budget_minus_R_part_max(5)}"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
