#!/usr/bin/env python3
"""Five NEW exact comparisons for cubic-frame template alignment.

The Gaussian enclosures and the squared sqrt(2) endpoint comparison
are reused from the previously verified 28-check baseline certificate.
This script does not recompute that certificate. No floating-point
arithmetic, optimization, matrix construction, or search is used.
Written for one remote verification run; initially UNRUN.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path


checks = []


def check(label, lhs, relation, rhs):
    lhs, rhs = F(lhs), F(rhs)
    passed = {"<": lhs < rhs, ">": lhs > rhs}[relation]
    record = {
        "label": label,
        "lhs": str(lhs),
        "relation": relation,
        "rhs": str(rhs),
        "passed": passed,
    }
    checks.append(record)
    if not passed:
        raise AssertionError(json.dumps(record, sort_keys=True))


# Reused strict bounds: P > p_lo, phi > phi_lo, 0 <= R < r_hi.
# The prior certificate also verifies positivity of 2 P^2 - 1 + 2 phi
# and (283/200)^2 > 2. None of its 28 checks is repeated here.
p_lo = F(68268, 100000)
phi_lo = F(24197, 100000)
r_hi = F(844, 10000)
q_lo = F(12, 5)
x_cap = F(83, 200)
d_coarse = F(4824, 10000)
b_coarse = F(63, 100)
penalty_cap = F(67, 1000)

d_lo = q_lo * (2 * p_lo ** 2 - 1 + 2 * phi_lo) - 1 + 2 * phi_lo
check("D lower threshold at q=12/5", d_lo, ">", d_coarse)

b_lo = 4 * q_lo * (q_lo - 1) * phi_lo ** 2 / (3 * x_cap)
check("cubic coefficient lower threshold", b_lo, ">", b_coarse)

check("maximal penalty cubed after clearing positive denominator",
      27 * r_hi ** 2, "<", 1024 * b_coarse * penalty_cap ** 3)

check("strict contradiction margin", d_coarse - penalty_cap, ">", x_cap)

check("weak-Dirac barrier squared", 17 ** 2 * 2, ">", 24 ** 2)

print(json.dumps({
    "status": "PASS",
    "arithmetic": "fractions.Fraction only",
    "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    "check_count": len(checks),
    "checks": checks,
    "reused_baseline_not_rerun": {
        "proof_filename": "original_mo_scalar_template_gamma_bound.md",
        "proof_sha256": "bd5997203c52895744a078048e206241996c46ef485e8975d7955b73be41f1c6",
        "checker_sha256": "d3af3d3bac9ba4d73a7589ba9ed4ff6261fde3263c64d04de36da7f36a1c65d3",
        "verified_result_sha256": "fbc10c4760d963f9364dca586cca3d8df5692ab786cd155634a651fac3a62d9d",
        "prior_check_count": 28,
        "statement": "Provenance of reused bounds; this script does not reverify the baseline.",
    },
}, indent=2, sort_keys=True))
