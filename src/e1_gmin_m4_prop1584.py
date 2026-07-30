#!/usr/bin/env python3
"""
Prop 15.84 — GD ⇒ cand via S3 budget (Max+-free algebra).

From Prop 15.77: at star_a = +1 on a |κ|=1 set,
  star·(S1+S3) = p ρ - 2/p²,
and GD (star·S1 ≤ 0) implies
  p ρ ≤ 2/p² + S3    ⇒    ρ ≤ 2/p³ + S3/p.

Cand residual target (Prop 15.74 / 15.83):
  ρ_cand = (p² - 4p - 3) / (p² (2p+3)).

Hence GD + S3 ≤ B_cand(p)  ⇒  ρ ≤ ρ_cand  ⇒  |m4| ≤ M_cand on same-sign sets, where
  B_cand(p) := p ρ_cand - 2/p²
            = (p³ - 4p² - 7p - 6) / (p² (2p+3)).

Proved algebra (this module):
  1) Closed form B_cand as above.
  2) B_cand(5) = -16/325 < 0  (need strictly negative S3 at p=5).
  3) B_cand(p) > 0 for all primes p≥7  (cubic p³-4p²-7p-6 > 0).
  4) Naive degree bound |S3| ≤ d3 * max|ρ| with d3=p²-5 cannot close
     cand by absolute bootstrap when 4p < d1 (already Prop 15.76); recorded.
  5) Diagonal dominance 4p > d1 fails for all primes p≥5 (d1=3p²-7).

Does **not** prove GD or bound S3 for Max+. L = lim α_n remains OPEN.

Backend: CPU Fraction only (F20 GPU unused). No moduli keys (F19).
Writes evidence/e1_gmin_m4_prop1584.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    d = 3
    while d * d <= p:
        if p % d == 0:
            return False
        d += 2
    return True


def M_cand(p: int) -> Fraction:
    return Fraction(p - 2, p * (2 * p + 3))


def rho_cand(p: int) -> Fraction:
    return Fraction(p * p - 4 * p - 3, p * p * (2 * p + 3))


def B_cand(p: int) -> Fraction:
    """S3 upper budget under GD for cand: p*ρ_cand - 2/p²."""
    return Fraction(p**3 - 4 * p**2 - 7 * p - 6, p * p * (2 * p + 3))


def d1_kappa1(p: int) -> int:
    """One-center / set degree into |κ|=1 (Prop 15.72/76): 3p²-7."""
    return 3 * p * p - 7


def d3_kappa1(p: int) -> int:
    """Degree into |κ|=3: p²-5."""
    return p * p - 5


def prove_B_cand_formula(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        # Two expressions for B_cand
        b1 = p * rho_cand(p) - Fraction(2, p * p)
        b2 = B_cand(p)
        rows[str(p)] = {
            "B_cand": str(b2),
            "match": b1 == b2,
            "B_cand_float": float(b2),
            "rho_cand": str(rho_cand(p)),
            "M_cand": str(M_cand(p)),
        }
        if b1 != b2:
            ok = False
    return {
        "proved_closed_form": ok,
        "formula": "(p^3 - 4p^2 - 7p - 6) / (p^2 (2p+3))",
        "by_p": rows,
    }


def prove_B_cand_sign(primes: list[int]) -> dict:
    """Sign of cubic num = p³-4p²-7p-6."""
    rows = {}
    ok = True
    # Algebra: for integer p≥7, p³-4p²-7p-6 ≥ 7³-4*49-49-6=343-196-49-6=92>0
    # p=5: 125-100-35-6=-16<0
    # Derivative 3p²-8p-7 > 0 for p≥5 ⇒ increasing on [5,∞)
    for p in primes:
        num = p**3 - 4 * p**2 - 7 * p - 6
        b = B_cand(p)
        row = {
            "cubic_num": num,
            "B_cand_sign": 1 if b > 0 else (-1 if b < 0 else 0),
            "B_cand": str(b),
        }
        if p == 5:
            row["expect_negative"] = True
            if not (b < 0 and num == -16):
                ok = False
            row["B_cand_eq_neg_16_325"] = b == Fraction(-16, 325)
        elif p >= 7:
            row["expect_positive"] = True
            if not (b > 0 and num > 0):
                ok = False
        rows[str(p)] = row
    return {
        "proved_sign_pattern": ok,
        "statement": (
            "B_cand(5)=-16/325<0; B_cand(p)>0 for all primes p≥7 "
            "(cubic increasing on [5,∞) and positive at p=7)"
        ),
        "by_p": rows,
        "consequence": (
            "At p=5, GD alone forces S3 ≤ negative budget to reach cand "
            "(need strong cancellation). For p≥7, GD + S3≤B_cand (positive) "
            "suffices for cand if S3 is controlled."
        ),
    }


def prove_diagonal_dominance_fails(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        d1 = d1_kappa1(p)
        gap = 4 * p - d1  # 4p - (3p²-7) = -3p²+4p+7
        rows[str(p)] = {
            "d1": d1,
            "d3": d3_kappa1(p),
            "4p_minus_d1": gap,
            "diagonal_dominance_fails": gap < 0,
        }
        if gap >= 0:
            ok = False
    # Algebra: -3p²+4p+7 < 0 for p≥5 (disc 16+84=100, roots ( -4±10)/-6; positive root 7/3≈2.33)
    return {
        "proved_fails_for_p_ge_5": ok,
        "algebra": "4p-d1 = -3p²+4p+7 < 0 for all p≥3",
        "by_p": rows,
        "note": "Absolute |T| row-sum bootstrap cannot invert 4pI-T on κ1",
    }


def main() -> None:
    from io_atomic import write_json_atomic

    primes = [p for p in range(5, 80) if is_prime(p)]
    out = {
        "prop": "15.84",
        "backend": "CPU algebra (Fraction); GPU unused (F20)",
        "L_status": "OPEN",
        "B_cand_formula": prove_B_cand_formula(primes),
        "B_cand_sign": prove_B_cand_sign(primes),
        "diag_dom": prove_diagonal_dominance_fails(primes),
        "settlement_lemma": (
            "If GD (star·S1≤0) and S3≤B_cand on every same-sign |κ|=1 centre, "
            "then ρ≤ρ_cand ⇒ max|m4|≤M_cand ⇒ bi-tight empty for p≥5. "
            "GD and S3 bound remain OPEN for general primes."
        ),
    }
    out["proved"] = (
        out["B_cand_formula"]["proved_closed_form"]
        and out["B_cand_sign"]["proved_sign_pattern"]
        and out["diag_dom"]["proved_fails_for_p_ge_5"]
    )
    path = ROOT / "evidence" / "e1_gmin_m4_prop1584.json"
    write_json_atomic(path, out)
    print(f"Prop 15.84 proved={out['proved']} primes={len(primes)} wrote {path}")
    print(f"  B_cand(5)={B_cand(5)}  B_cand(7)={B_cand(7)}")
    print(f"  L OPEN — GD/S3 still unproved for ∀p≥5")


if __name__ == "__main__":
    main()
