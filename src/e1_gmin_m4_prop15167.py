#!/usr/bin/env python3
"""Prop 15.167 — spectral majorization algebra; bi-tight conclusion retracted.

The Fraction identities in Theorems A--B are valid: conditional on a spectral
floor lambda_min(Phi)>=6, majorization gives L_*<2d and hence the cycle
spectrum lies below d.  The old final arrow through Proposition 15.55 was
false.  If R=G-(n/2)P_1, then

    ker R = span{1} + ker G,

not merely span{1}; Proposition 15.56 itself supplies n-2 star-difference
vectors in ker G.  Thus a tight centered indicator can lie in ker G even
when the top eigenvalue n/2 is simple.  This module now records the valid
conditional spectral algebra and returns False for the retracted bi-tight
claim. Proposition 15.720 supplies the valid discrete replacement for the
required levels 2 and 3; bi-tight level 4 is only a corollary and does not
exclude one-sided tightness.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7, dim_Z  # noqa: E402


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


def L_star_bulk6(p: int) -> Fraction:
    """Majorization UB: mult≥d−1, λ_min≥6, tr=n(n−2) ⇒ λ_max ≤ L_*."""
    n = n_of(p)
    d = d_of(p)
    m = dim_Z(p)
    k = d - 1
    T = n * (n - 2)
    return Fraction(T - 6 * (m - k), k)


def L_star_closed(p: int) -> Fraction:
    """Closed form (p⁴ + 24 p² − 1)/(2(p² − 1))."""
    return Fraction(p**4 + 24 * p * p - 1, 2 * (p * p - 1))


def two_d_minus_L_star(p: int) -> Fraction:
    """2d − L_* = (p⁴ − 24 p² − 1)/(2(p² − 1)) > 0 for p≥5."""
    return Fraction(p**4 - 24 * p * p - 1, 2 * (p * p - 1))


def lambda_cycle_ub_from_L(L: Fraction | float) -> Fraction:
    """λ_cycle = λ_max(Φ)/2 (15.61 normalization)."""
    return Fraction(L) / 2


def bitight_from_majorization(p: int) -> dict:
    """
    Record the valid arithmetic and the retracted final implication.
    """
    ok_p = p >= 5 and is_prime(p)
    L = L_star_closed(p)
    gap = two_d_minus_L_star(p)
    d = d_of(p)
    cycle_ub = lambda_cycle_ub_from_L(L)
    spectral_gap_conditional = bool(ok_p and gap > 0 and cycle_ub < d)
    return {
        "p": p,
        "prime_ge_5": ok_p,
        "L_star": str(L),
        "two_d": 2 * d,
        "two_d_minus_L_star": str(gap),
        "L_star_lt_2d": gap > 0,
        "lambda_cycle_ub": str(cycle_ub),
        "d": d,
        "lambda_cycle_ub_lt_d": cycle_ub < d,
        "spectral_gap_conditional_on_floor": spectral_gap_conditional,
        "bitight_empty": False,
        "retracted": True,
        "theorem": (
            "15.167 A-B: conditional spectral gap arithmetic is valid; "
            "bi-tight does not follow because ker(G-(n/2)P1) contains ker G."
        ),
    }


def bitight_empty_for_all_primes_ge_5(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 200) if is_prime(p)]
    rows = {str(p): bitight_from_majorization(p) for p in primes}
    ok = False
    return {
        "proved": ok,
        "n_checked": len(primes),
        "by_p_sample": {k: rows[k] for k in list(rows)[:6]},
        "all_bitight_empty": ok,
    }


def prove_theorem_A(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 100) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        a = L_star_bulk6(p)
        b = L_star_closed(p)
        rows[str(p)] = {
            "L_star_bulk6": str(a),
            "L_star_closed": str(b),
            "match": a == b,
        }
        if a != b:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "L_* = (T − 6(m − (d−1)))/(d−1) = (p⁴+24p²−1)/(2(p²−1)); "
            "majorization UB on λ_max under mult≥d−1 and λ_min≥6."
        ),
        "by_p_sample": {k: rows[k] for k in list(rows)[:5]},
        "n_checked": len(primes),
    }


def prove_theorem_B(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 200) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        gap = two_d_minus_L_star(p)
        # poly p^4 - 24p^2 - 1
        poly = p**4 - 24 * p * p - 1
        rows[str(p)] = {
            "gap": str(gap),
            "poly": poly,
            "gap_positive": gap > 0,
        }
        if not (gap > 0 and poly > 0):
            ok = False
    # p=5 base + monotonicity of f(x)=x^2-24x-1 for x>=25
    f5 = 25**2 - 24 * 25 - 1  # x=p^2
    return {
        "proved": ok and f5 == 24,
        "theorem": (
            "2d−L_*=(p⁴−24p²−1)/(2(p²−1)); numerator >0 for primes p≥5 "
            "(value 24 at p=5; f(x)=x²−24x−1 increasing on x≥25)."
        ),
        "f_at_p2_5": f5,
        "by_p_sample": {k: rows[k] for k in list(rows)[:5]},
        "n_checked": len(primes),
    }


def prove_theorem_C() -> dict:
    bt = bitight_empty_for_all_primes_ge_5()
    return {
        "proved": False,
        "retracted": True,
        "theorem": (
            "The old arrow λ_cycle<d ⇒ bi-tight empty is invalid: "
            "ker(G-(n/2)P1)=span{1}+ker G, and ker G is nontrivial."
        ),
        "bitight_check": bt,
    }


def certify_census_spectrum() -> dict:
    """Census: actual λ_max(Φ) ≤ L_* and cycle < d at p=5,7."""
    out = {}
    for p, spec in ((5, SPECTRUM_P5), (7, SPECTRUM_P7)):
        lam = max(spec)
        mult = spec[lam]
        d = d_of(p)
        L = L_star_closed(p)
        cycle = Fraction(lam) / 2
        out[str(p)] = {
            "lambda_max_Phi": str(lam),
            "mult": mult,
            "d_minus_1": d - 1,
            "mult_ge_d_minus_1": mult >= d - 1,
            "lambda_min": str(min(spec)),
            "lambda_min_ge_6": min(spec) >= 6,
            "L_star": str(L),
            "lambda_max_le_L_star": lam <= L,
            "lambda_cycle": str(cycle),
            "d": d,
            "lambda_cycle_lt_d": cycle < d,
            "spectral_conditions": cycle < d and mult >= d - 1 and min(spec) >= 6,
            "bitight_empty": False,
        }
    return {
        "certified_spectral_conditions": all(r["spectral_conditions"] for r in out.values()),
        "certified": False,
        "by_p": out,
    }


def prove_open() -> dict:
    return {
        "bi_tight_empty_for_all_p_ge_5": False,
        "residual_closed_general": False,
        "sixteen_N_for_all_p": False,
        "Es4_star_for_all_p": False,
        "E1_closed": False,  # E(1) is Prop 15.168, not this module
        "L_status": "OPEN",  # this module alone does not close L; see 15.168
        "note": (
            "15.167 bi-tight implication retracted. Required levels 2 and 3 "
            "are instead closed by the discrete degree congruence in 15.720."
        ),
    }


def main() -> dict:
    A = prove_theorem_A()
    B = prove_theorem_B()
    C = prove_theorem_C()
    cert = certify_census_spectrum()
    open_ = prove_open()
    out = {
        "title": (
            "Prop 15.167 spectral majorization algebra; bi-tight conclusion retracted"
        ),
        "L_status": "OPEN",
        "proved": {
            "L_star_closed_form": A["proved"],
            "L_star_lt_2d_for_p_ge_5": B["proved"],
            "bitight_empty_for_all_p_ge_5": C["proved"],
            "census_spectral_conditions_p5_p7": cert["certified_spectral_conditions"],
            "residual_closed_general": False,
            "sixteen_N_for_all_p": False,
            "E1_closed": False,
        },
        "algebra": {"A": A, "B": B, "C": C, "open": open_},
        "census": cert,
        "bitight_predicate": bitight_from_majorization(5),
        "F3": "no soft-close of L (E(1)/Main still open)",
        "F19": "no class_key thrash",
        "F20": "CPU Fraction only for general-p claim",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15167.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.167 spectral algebra (bi-tight implication retracted)")
    print(f"  A L_* closed form: {A['proved']}")
    print(f"  B L_* < 2d for p≥5: {B['proved']}")
    print(f"  C bi-tight empty all p≥5: {C['proved']} (retracted)")
    print(f"  census p=5,7 spectral conditions: {cert['certified_spectral_conditions']}")
    for p, r in cert["by_p"].items():
        print(
            f"    p={p} λmax={r['lambda_max_Phi']} ≤ L_*={r['L_star']} "
            f"cycle={r['lambda_cycle']} < d={r['d']}"
        )
    print(f"  residual_closed_general={open_['residual_closed_general']}")
    print(f"  bi_tight_empty_for_all_p_ge_5={open_['bi_tight_empty_for_all_p_ge_5']}")
    print(f"  L_status={open_['L_status']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
