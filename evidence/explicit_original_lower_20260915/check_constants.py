"""Exact arithmetic for the explicit original lower bound.

No signings, spectral profiles, or spin states are enumerated. This does
not verify the analytic proof or the imported local Gaussianization lemma.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import platform


def main():
    eps = F(1, 10**6)
    k_lo, k_hi = F(7, 11), F(16, 25)
    s_lo, s_hi = F(797, 1000), F(4, 5)
    r_hi = F(27, 10000)
    p = F(1, 10)
    checks = {
        "update_alpha_coefficient": 100 * (1 + p*p) == 101,
        "update_baseline_coefficient": 100 * (1-p)**2 == 81,
        "update_field_coefficient": 100 * p * (1-p) == 9,
        "sqrt_kappa_lower": s_lo*s_lo < k_lo,
        "sqrt_kappa_upper": s_hi*s_hi == k_hi,
        "flat_gain_lower": 9*s_lo - 10*k_hi == F(773, 1000),
        "defect_linear_coefficient": F(81, 4)*k_hi < 13,
        "defect_sqrt_coefficient": 162*k_hi + 63*s_hi < 155,
        "nuclear_defect_factor": 4/k_lo < 7,
        "defect_sqrt_upper": 7*eps < r_hi*r_hi,
    }
    margin = F(773, 1000) - 13*7*eps - 155*r_hi - 101*eps
    checks["strict_contradiction_margin"] = margin == F(88577, 250000) > 0
    result = {
        "classification": "exact scalar corroboration; analytic theorem not formalized",
        "python": platform.python_version(),
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "epsilon": str(eps),
        "margin": str(margin),
        "checks": checks,
        "all_passed": all(checks.values()),
        "convergence_proved": False,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
