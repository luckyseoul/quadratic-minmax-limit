#!/usr/bin/env python3
"""
Prop 15.167 — Bi-tight empty for all primes p≥5 via mult≥d−1 + λ_min≥6
majorization (no residual / 16N required).

Bypasses Path-C residual for the bi-tight link. Residual / 16N / E[s⁴]≤Es4_*
remain OPEN (honest). L OPEN until E(1)/deep ND + Main.

============================================================================
Setup. Φ on Z = zero-diag traceless forms on V₊, dim m=d(d−3)/2,
d=n/2, n=p²+1.  tr(Φ)=n(n−2).  λ_min(Φ)≥6 (15.160.D / residual_h_close
Thm D: Q₄≥0).  mult(λ_max(Φ))≥d−1 (15.162.B / 15.98: PSL min irrep).

From 15.61: λ_cycle = λ_max(Φ)/2  (same Rayleigh: Q_max = N λ_max(Φ)‖B‖²,
λ_cycle = 2 d²/N · λ₂(W), and λ₂(W)=N λ_max(Φ)/(4 d²)).

From 15.55–15.56: λ_max(G)=max(n/2, λ_cycle) and
  λ_max(G)=n/2  ⇔  λ₂(P⊙P) ≤ d/(2N)  ⇔  λ_cycle ≤ d.
When λ_cycle < n/2 the n/2-eigenspace is simple (star differences in ker G;
cycle space strictly below), so bi-tight covers empty (15.55).

============================================================================
Theorem A (majorization UB on λ_max) — PROVED Fraction.
  Let T=n(n−2), k=d−1, m=d(d−3)/2, ℓ=6.  Any spectrum with sum T,
  at least k coordinates equal to the top value L, and all coordinates ≥ℓ,
  satisfies
      L ≤ L_*(p) := (T − ℓ(m−k))/k
  with equality only for the two-level spectrum (k copies of L, rest =ℓ).
  Closed form:
      L_*(p) = (p⁴ + 24 p² − 1) / (2(p² − 1)).
  Certified identity for primes p=5..97.  ∎

Theorem B (L_* < 2d for primes p≥5) — PROVED Fraction.
  2d − L_* = (p⁴ − 24 p² − 1)/(2(p²−1)).
  Numerator p⁴−24p²−1 at p=5 equals 24>0; derivative of f(x)=x²−24x−1
  (x=p²) is 2x−24>0 for x≥25, so f(p²)>0 for all primes p≥5.
  Hence L_*(p) < 2d for all primes p≥5.  ∎

Theorem C (bi-tight empty for all primes p≥5) — PROVED.
  mult≥d−1 + λ_min≥6 + tr=T ⇒ λ_max(Φ) ≤ L_*(p) < 2d.
  Therefore λ_cycle = λ_max(Φ)/2 < d = n/2.
  Hence λ₂(P⊙P) < d/(2N), λ_max(G)=n/2 simple, bi-tight empty (15.55).
  **Does not use residual, Es4_*, or 16N.**  ∎

Theorem D (status, honest).
  bi_tight_empty_for_all_p_ge_5 = True (this prop).
  residual_closed_general = False (16N / Es4_* / H still open; not needed for bi-tight).
  E(1)_closed = False; Main / L OPEN until deep ND + denseness wiring.
  Census: λ_max(Φ) ≤ L_* and λ_cycle < d at p=5,7 (spectrum Fraction).  ∎

============================================================================
Shipped API: L_star_bulk6, two_d_minus_L_star, bitight_from_majorization,
bitight_empty_for_all_primes_ge_5.
Writes evidence/e1_gmin_m4_prop15167.json
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
    Predicate: prime p≥5 ⇒ bi-tight empty via L_* < 2d.
    Pure Fraction; no Max+ census required.
    """
    ok_p = p >= 5 and is_prime(p)
    L = L_star_closed(p)
    gap = two_d_minus_L_star(p)
    d = d_of(p)
    cycle_ub = lambda_cycle_ub_from_L(L)
    bitight = bool(ok_p and gap > 0 and cycle_ub < d)
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
        "bitight_empty": bitight,
        "theorem": "15.167.C mult≥d−1 + λ_min≥6 ⇒ bi-tight empty (no residual)",
    }


def bitight_empty_for_all_primes_ge_5(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 200) if is_prime(p)]
    rows = {str(p): bitight_from_majorization(p) for p in primes}
    ok = all(r["bitight_empty"] for r in rows.values())
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
        "proved": bt["proved"],
        "theorem": (
            "For every prime p≥5: λ_max(Φ)≤L_*<2d ⇒ λ_cycle<d=n/2 "
            "⇒ λ_max(G)=n/2 simple ⇒ bi-tight empty (15.55). "
            "No residual / Es4_* / 16N used."
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
            "bitight_empty": cycle < d and mult >= d - 1 and min(spec) >= 6,
        }
    return {
        "certified": all(r["bitight_empty"] for r in out.values()),
        "by_p": out,
    }


def prove_open() -> dict:
    return {
        "bi_tight_empty_for_all_p_ge_5": True,
        "residual_closed_general": False,
        "sixteen_N_for_all_p": False,
        "Es4_star_for_all_p": False,
        "E1_closed": False,  # E(1) is Prop 15.168, not this module
        "L_status": "OPEN",  # this module alone does not close L; see 15.168
        "note": (
            "Bi-tight closed by majorization. Residual/16N still open (optional). "
            "E(1)+L closed in Prop 15.168 + denseness Prop 6.2."
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
            "Prop 15.167 bi-tight empty for all primes p≥5 via "
            "mult≥d−1 + λ_min≥6 majorization (no residual)"
        ),
        "L_status": "OPEN",
        "proved": {
            "L_star_closed_form": A["proved"],
            "L_star_lt_2d_for_p_ge_5": B["proved"],
            "bitight_empty_for_all_p_ge_5": C["proved"],
            "census_spectrum_p5_p7": cert["certified"],
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
    print("Prop 15.167 bi-tight via majorization (no residual)")
    print(f"  A L_* closed form: {A['proved']}")
    print(f"  B L_* < 2d for p≥5: {B['proved']}")
    print(f"  C bi-tight empty all p≥5: {C['proved']}")
    print(f"  census p=5,7: {cert['certified']}")
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
