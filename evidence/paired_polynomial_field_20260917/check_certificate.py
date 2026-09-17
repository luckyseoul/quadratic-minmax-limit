"""Exact scalar certificate for the new paired-field bound; no census.

This is not a formal verification of the analytic theorem or its imports.
Run once on a mesh CPU after reviewing the proof and this implementation.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import platform

import sympy as sp


def main():
    a, b, q, y = sp.symbols("a b q y", real=True)
    v = a + b*y*y
    # After the first squaring, only v < 2a needs a second squaring.
    radicand_gap = sp.expand(v*v - 4*b*b*q*y*y - (2*a-v)**2)
    polynomial_certificate = 4*b*y*y*(a-b*q)
    k_lo = F(7, 11)
    field_lo = F(791, 1000)
    field_squared_lo = k_lo*(F(2, 3)+k_lo/2)
    p = F(1, 10)
    bound_lo = ((1-p)**2*k_lo/2 + p*(1-p)*field_lo)/(1+p*p)
    margin = bound_lo-F(13, 40)
    checks = {
        "paired_root_polynomial_identity": sp.expand(
            radicand_gap-polynomial_certificate) == 0,
        "field_squared_lower": field_squared_lo == F(455, 726),
        "field_lower_squared": field_lo**2 < field_squared_lo,
        "update_alpha_coefficient": 100*(1+p*p) == 101,
        "update_baseline_coefficient": 100*(1-p)**2 == 81,
        "update_field_coefficient": 100*p*(1-p) == 9,
        "strict_global_bound": margin > 0,
        # f''(3/4)^2 = 2304/343 < 9; magnitude is monotone.
        "arcsine_second_derivative_bound": F(2304, 343) < 9,
        # q>=1/2 implies (1-q)/(1+q)<=1/3.
        "positive_completed_square_cushion": F(2, 3)-F(16, 25)/3 > 0,
    }
    result = {
        "classification": "exact scalar and polynomial corroboration",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "polynomial_certificate": str(sp.factor(radicand_gap)),
        "polynomial_sign_hypotheses": "a>=b*q; b,q,y>=0",
        "field_squared_lower": str(field_squared_lo),
        "conservative_bound_lower": str(bound_lo),
        "margin_above_13_over_40": str(margin),
        "checks": checks,
        "all_passed": all(checks.values()),
        "signings_enumerated": 0,
        "convergence_proved": False,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
