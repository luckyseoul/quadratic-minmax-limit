#!/usr/bin/env python3
"""Exact rational certificate for the explicit lower bound and its sharpening.

Chain certified here (all scalar arithmetic exact; the analytic estimates are
cited from ``evidence/NOTE_2026-09-15_EXPLICIT_ORIGINAL_LOWER_BOUND.md``):

  (2)  alpha >= kappa q / (2 ell)                       [all-law note, (4.1)]
  (10) 101 alpha >= (81/2) kappa ell + 9 sqrt(kappa)
                    - (162 kappa + 63 sqrt(kappa)) sqrt(d) - r(n),  r(n) -> 0
  (9)  the update inequality with p = 1/10 and the Boolean bound
       |Q_M(z)| <= n alpha.

Reduction. Suppose some subsequence has alpha -> a <= kappa/2 + eps. Then
ell >= kappa q/(2 alpha) forces d <= 4 eps / (kappa + 2 eps) (using q <= 1),
and (10) together with ell = (q + 1 - d)/2 = 1 - (1/n + d)/2 gives

  101 (a - kappa/2) >= 9 sqrt(kappa) - 10 kappa
                        - (81/4) kappa d
                        - (162 kappa + 63 sqrt(kappa)) sqrt(d) - o(1).

A contradiction (hence liminf alpha > kappa/2 + eps) follows if the enclosed
right-hand side exceeds 101 eps. Enclosures used:

  * pi in (31415926/10**7, 31415927/10**7), the frozen enclosure recorded in
    NOTE_2026-09-05_SOURCE_CROSS_NUCLEAR_TRACE_BOUNDARY.md;
  * kappa = 2/pi,  sqrt(kappa) bounded by exact integer squares;
  * d <= 4 eps/(kappa + 2 eps) <= 6.3 eps  (4/kappa < 6.285);
  * (81/4) kappa < 12.9;  162 kappa + 63 sqrt(kappa) < 153.4.

This is the SAME method as the note with tighter rational bookkeeping and the
sharper d-bound 6.3 eps (the note used 7 eps); no new analytic estimate is
introduced. Nothing here proves or verifies the analytic estimates themselves.

Outputs a JSON receipt.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

# ---- frozen pi enclosure (repo note, 2026-09-05) --------------------------
PI_LO = F(31415926, 10 ** 7)
PI_HI = F(31415927, 10 ** 7)
KAPPA_LO = 2 / PI_HI
KAPPA_HI = 2 / PI_LO


def sqrt_upper(x: F, denom: int = 10 ** 10) -> F:
    """Smallest stated-denominator rational r with r >= sqrt(x), x > 0."""
    scaled = x * denom * denom
    s = math.isqrt(scaled.numerator // scaled.denominator)
    while F(s, denom) ** 2 < x:
        s += 1
    return F(s, denom)


def enclosures() -> dict:
    s_lo = F(797884, 10 ** 6)
    s_hi = F(797885, 10 ** 6)
    checks = {
        "pi_bounds_ordered": PI_LO < PI_HI,
        "sqrt_kappa_lower_exact": s_lo * s_lo < KAPPA_LO,
        "sqrt_kappa_upper_exact": s_hi * s_hi > KAPPA_HI,
    }
    return {
        "pi_lo": str(PI_LO), "pi_hi": str(PI_HI),
        "kappa_lo": str(KAPPA_LO), "kappa_hi": str(KAPPA_HI),
        "kappa_lo_float": float(KAPPA_LO), "kappa_hi_float": float(KAPPA_HI),
        "sqrt_kappa_lo": str(s_lo), "sqrt_kappa_hi": str(s_hi),
        "checks": checks,
    }


def certify(eps: F, tight: bool = True) -> dict:
    """Certificate for the statement liminf alpha > kappa/2 + eps."""
    s_lo = F(797884, 10 ** 6)
    s_hi = F(797885, 10 ** 6)
    if tight:
        k_lo, k_hi = KAPPA_LO, KAPPA_HI
    else:
        k_lo, k_hi = F(7, 11), F(16, 25)   # the note's own enclosure

    # Step 1: d-bound.  d <= 4 eps/(kappa + 2 eps), certified <= D_COEF * eps.
    d_coef = F(63, 10) if tight else F(7, 1)
    d_bound_ok = 4 * eps <= d_coef * eps * (k_lo + 2 * eps)
    d_max = d_coef * eps

    # Step 2: income and penalties with rational enclosures.
    income_lo = 9 * s_lo - 10 * k_hi
    lin_coef = F(129, 10) if tight else F(13, 1)
    sqrt_coef = F(1534, 10) if tight else F(155, 1)
    r = sqrt_upper(d_max)                      # r >= sqrt(d_max)
    lin_pen = lin_coef * d_max
    sqrt_pen = sqrt_coef * r
    lower = income_lo - lin_pen - sqrt_pen
    margin = lower - 101 * eps

    return {
        "eps": str(eps), "eps_float": float(eps),
        "enclosure": "tight" if tight else "note_enclosure_(7/11,16/25)",
        "d_bound_ok": bool(d_bound_ok),
        "d_max": str(d_max), "sqrt_d_max_upper": str(r),
        "income_lo": str(income_lo),
        "linear_pen": str(lin_pen), "sqrt_pen": str(sqrt_pen),
        "rhs_lower": str(lower),
        "101_eps": str(101 * eps),
        "margin": str(margin), "margin_float": float(margin),
        "contradiction": bool(margin > 0),
    }


def ceiling(lo: F = F(1, 10 ** 6), hi: F = F(2, 10 ** 5),
            steps: int = 200) -> dict:
    """Largest grid eps (step 1e-8) with a certified contradiction."""
    step = F(1, 10 ** 8)
    # bisect on the monotone predicate margin(eps) > 0
    a, b = lo, hi
    while b - a > step:
        mid = (a + b) / 2
        if certify(F(mid).limit_denominator(10 ** 12))["contradiction"]:
            a = mid
        else:
            b = mid
    eps = F(a).limit_denominator(10 ** 8)
    return {"ceiling_eps": str(eps), "ceiling_float": float(eps)}


def note_exact_reproduction() -> dict:
    """Reproduce the note's own published margin 88577/250000 exactly.

    Uses precisely the note's stated constants: sqrt(kappa) > 797/1000,
    kappa < 16/25, (162 kappa + 63 sqrt kappa) < 155, (81/4) kappa < 13,
    sqrt(7 eps) < 27/10000 (their sqrt(7e-6) bound), d <= 7 eps.
    """
    eps = F(1, 10 ** 6)
    margin = F(773, 1000) - 13 * 7 * eps - 155 * F(27, 10000) - 101 * eps
    return {
        "eps": str(eps),
        "margin": str(margin),
        "margin_float": float(margin),
        "matches_published": margin == F(88577, 250000),
        "contradiction": margin > 0,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eps", default=None,
                    help="comma list of rationals/decimals to certify")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    claimed = [F(s) for s in (args.eps.split(",") if args.eps else
                              ["1/1000000", "2/1000000", "3/1000000",
                               "7/2000000", "4/1000000"])]
    controls = [F(9, 2000000), F(1, 200000)]  # beyond ceiling, must fail

    result = {
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "enclosures": enclosures(),
        "claimed_certificates": [certify(e) for e in claimed],
        "control_certificates": [certify(e) for e in controls],
        "note_exact_reproduction": note_exact_reproduction(),
    }
    result["ceiling"] = ceiling()
    result["all_claimed_certified"] = all(
        c["contradiction"] and c["d_bound_ok"]
        for c in result["claimed_certificates"])
    result["controls_fail_as_expected"] = all(
        not c["contradiction"] for c in result["control_certificates"])
    result["ceiling_certified"] = certify(
        F(result["ceiling"]["ceiling_eps"]))["contradiction"]
    result["headline_eps"] = str(max(claimed))

    text = json.dumps(result, indent=1, sort_keys=True)
    if args.out:
        args.out.write_text(text + "\n")
        print(f"wrote {args.out}")
    else:
        print(text)
    return 0 if (result["all_claimed_certified"]
                 and result["controls_fail_as_expected"]
                 and result["ceiling_certified"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
