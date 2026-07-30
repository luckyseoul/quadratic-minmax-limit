#!/usr/bin/env python3
"""
Prop 15.83 — Max+-free algebra: resolvent-budget hierarchy for M_cand vs L.

Proves (Fraction exact, no Max+ samples):
  1) Cascade M_cand < M_mid ≤ L_abs < T_abs for all primes p≥5.
  2) Same-sign residual targets:
        ρ_L   = (p-4)/(2 p²)          # for |m4|≤L_abs via Wick + residual
        ρ_cand= (p²-4p-3)/(p²(2p+3)) # for |m4|≤M_cand
     satisfy 0 < ρ_cand < ρ_L for primes p≥5.
  3) Resolvent gain targets (source amp 24/p² on |κ|=3):
        gain_L   = (p-4)/48
        gain_cand= (p²-4p-3)/(24(2p+3))
     satisfy 0 ≤ gain_cand < gain_L for all primes p≥5
     (equality never; limit ratio → 1 as p→∞).

Interpretation: any resolvent-gain bound ≤ gain_cand automatically gives
max|m4|≤M_cand≤L_abs and thus bi-tight empty (Prop 15.47) for p≥5.
This does **not** prove the gain bound; it only ranks the residual targets.
L = lim α_n remains OPEN.

F19: no moduli class keys. F20: pure algebra (CPU; GPU unused).
Writes evidence/e1_gmin_m4_prop1583.json
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


def M_mid(p: int) -> Fraction:
    return Fraction(p - 2, 2 * p * (p + 1))


def L_abs(p: int) -> Fraction:
    return Fraction(p - 2, 2 * p * p)


def T_abs(p: int) -> Fraction:
    return Fraction(p - 2, p * (2 * p - 1))


def rho_L(p: int) -> Fraction:
    """Same-sign residual thr for |m4| ≤ L_abs on |κ|=1 (Prop 15.66/67)."""
    return Fraction(p - 4, 2 * p * p)


def rho_cand(p: int) -> Fraction:
    """Same-sign residual thr for |m4| ≤ M_cand (Prop 15.74)."""
    # M_cand - 1/p² = (p-2)/(p(2p+3)) - 1/p²
    # = [(p-2)p - (2p+3)] / [p²(2p+3)] = (p² - 2p - 2p - 3) / ...
    # = (p² - 4p - 3) / [p²(2p+3)]
    return Fraction(p * p - 4 * p - 3, p * p * (2 * p + 3))


def gain_L(p: int) -> Fraction:
    """Resolvent gain budget for L (Prop 15.68/72): (p-4)/48."""
    return Fraction(p - 4, 48)


def gain_cand(p: int) -> Fraction:
    """Resolvent gain budget for cand (Prop 15.74): (p²-4p-3)/(24(2p+3))."""
    return Fraction(p * p - 4 * p - 3, 24 * (2 * p + 3))


def prove_cascade(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        Mc, Mm, La, Ta = M_cand(p), M_mid(p), L_abs(p), T_abs(p)
        row = {
            "M_cand": str(Mc),
            "M_mid": str(Mm),
            "L_abs": str(La),
            "T_abs": str(Ta),
            "Mc_lt_Mm": Mc < Mm,
            "Mm_le_La": Mm <= La,
            "La_lt_Ta": La < Ta,
        }
        if not (row["Mc_lt_Mm"] and row["Mm_le_La"] and row["La_lt_Ta"]):
            ok = False
        rows[str(p)] = row
    # Algebraic identities (p-independent rational functions)
    # M_mid - M_cand = (p-2)/(2p(p+1)) - (p-2)/(p(2p+3))
    # = (p-2)/p * [1/(2(p+1)) - 1/(2p+3)]
    # = (p-2)/p * [(2p+3) - 2(p+1)] / [2(p+1)(2p+3)]
    # = (p-2)/p * (2p+3-2p-2) / denom = (p-2)/p * 1 / [2(p+1)(2p+3)] > 0 for p>2
    # L - M_mid = (p-2)/(2p²) - (p-2)/(2p(p+1)) = (p-2)/(2p) [1/p - 1/(p+1)]
    # = (p-2)/(2p) * 1/(p(p+1)) > 0
    # T - L = (p-2)/p * [1/(2p-1) - 1/(2p)] = (p-2)/p * 1/(2p(2p-1)) > 0
    return {
        "proved_for_listed_primes": ok,
        "algebra_note": (
            "M_mid-M_cand=(p-2)/(2p(p+1)(2p+3))>0; "
            "L-M_mid=(p-2)/(2p²(p+1))>0; "
            "T-L=(p-2)/(2p²(2p-1))>0 for p>2"
        ),
        "by_p": rows,
    }


def prove_gain_hierarchy(primes: list[int]) -> dict:
    """Prove gain_cand < gain_L for all primes p≥5 by rational comparison."""
    rows = {}
    ok = True
    # Algebra: gain_L - gain_cand > 0
    # (p-4)/48 - (p²-4p-3)/(24(2p+3))
    # = (p-4)/(48) - 2(p²-4p-3)/(48(2p+3))
    # = [ (p-4)(2p+3) - 2(p²-4p-3) ] / [48(2p+3)]
    # num = (p-4)(2p+3) - 2p² + 8p + 6 = 2p² + 3p - 8p - 12 - 2p² + 8p + 6
    #     = (2p²-2p²) + (3p-8p+8p) + (-12+6) = 3p - 6 = 3(p-2)
    # So gain_L - gain_cand = 3(p-2) / [48(2p+3)] > 0 for p>2.
    alg_diff = "3(p-2)/(48(2p+3))"
    for p in primes:
        gL, gC = gain_L(p), gain_cand(p)
        rL, rC = rho_L(p), rho_cand(p)
        # ρ_cand > 0 iff p²-4p-3 > 0 iff (p-1)(p-3)>? p²-4p-3=(p-1)(p-3)-6+3 wait
        # roots of p²-4p-3: 2±√7; positive for p≥5
        row = {
            "gain_L": str(gL),
            "gain_cand": str(gC),
            "gain_diff": str(gL - gC),
            "gain_diff_closed": str(Fraction(3 * (p - 2), 48 * (2 * p + 3))),
            "gain_cand_lt_gain_L": gC < gL,
            "gain_cand_nonneg": gC >= 0,
            "rho_L": str(rL),
            "rho_cand": str(rC),
            "rho_cand_lt_rho_L": rC < rL,
            "rho_cand_pos": rC > 0,
            "ratio_gain_cand_over_L": float(gC / gL) if gL != 0 else None,
        }
        # exact match of closed form
        if gL - gC != Fraction(3 * (p - 2), 48 * (2 * p + 3)):
            ok = False
            row["closed_form_match"] = False
        else:
            row["closed_form_match"] = True
        if not (
            row["gain_cand_lt_gain_L"]
            and row["gain_cand_nonneg"]
            and row["rho_cand_lt_rho_L"]
            and row["rho_cand_pos"]
        ):
            ok = False
        rows[str(p)] = row
    return {
        "proved_algebra": True,
        "gain_L_minus_gain_cand": alg_diff,
        "proof_sketch": (
            "gain_L-gain_cand = (p-4)/48 - (p²-4p-3)/(24(2p+3)) "
            "= [(p-4)(2p+3) - 2(p²-4p-3)] / [48(2p+3)] "
            "= 3(p-2)/[48(2p+3)] > 0 for p>2. "
            "ρ_cand=M_cand-1/p²=(p²-4p-3)/(p²(2p+3)); "
            "ρ_L-ρ_cand = [(p-4)(2p+3) - 2(p²-4p-3)]/(2p²(2p+3)) "
            "= 3(p-2)/(2p²(2p+3)) > 0."
        ),
        "certified_primes": ok,
        "by_p": rows,
        "note": (
            "If resolvent gain ≤ gain_cand then same-sign |ρ|≤ρ_cand "
            "⇒ |m4|≤M_cand≤L_abs ⇒ g_min≥-M_cand>T(p) ⇒ bi-tight empty "
            "(Props 15.47, 15.66–15.74). Gain bound itself remains OPEN for ∀p≥5."
        ),
    }


def main() -> None:
    from io_atomic import write_json_atomic

    primes = [p for p in range(5, 80) if is_prime(p)]
    out = {
        "prop": "15.83",
        "backend": "CPU algebra (Fraction); GPU unused (F20)",
        "L_status": "OPEN",
        "cascade": prove_cascade(primes),
        "gain_hierarchy": prove_gain_hierarchy(primes),
    }
    # Single boolean for tests
    out["proved"] = (
        out["cascade"]["proved_for_listed_primes"]
        and out["gain_hierarchy"]["proved_algebra"]
        and out["gain_hierarchy"]["certified_primes"]
    )
    path = ROOT / "evidence" / "e1_gmin_m4_prop1583.json"
    write_json_atomic(path, out)
    print(f"Prop 15.83 proved={out['proved']} primes={len(primes)} wrote {path}")
    print(f"  gain_L-gain_cand = {out['gain_hierarchy']['gain_L_minus_gain_cand']}")
    print(f"  L remains OPEN — gain bound not proved")


if __name__ == "__main__":
    main()
