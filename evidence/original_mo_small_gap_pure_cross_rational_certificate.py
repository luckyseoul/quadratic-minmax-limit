#!/usr/bin/env python3
"""Fixed rational certificate for one analytic pure-cross upper bound.

This script is not a signing search, SDP, optimization, or parameter scan.
It checks predetermined rational enclosures using exact Fraction arithmetic.
Run only on the authorized offload host, not on the local coordinator.
"""

from fractions import Fraction as F
import json


checks = {}


def check(name, predicate):
    checks[name] = bool(predicate)
    if not predicate:
        raise AssertionError(name)


# Reuse the already-proved baseline pi interval; do not repeat Machin.
# Baseline checker SHA256:
# d3af3d3bac9ba4d73a7589ba9ed4ff6261fde3263c64d04de36da7f36a1c65d3
# Baseline result SHA256:
# fbc10c4760d963f9364dca586cca3d8df5692ab786cd155634a651fac3a62d9d
pi_lower = F(31415926, 10**7)
pi_upper = F(31415927, 10**7)
kappa_lower = 2 / pi_upper
kappa_upper = 2 / pi_lower
check("kappa_coarse_interval", F(3, 5) < kappa_lower < kappa_upper < F(2, 3))

t = F(3, 5)
r = t * t
s_lower = kappa_lower**2
s_upper = kappa_upper**2
m_lower = s_lower / 2
m_upper = s_upper / 2

# A_s(m) = [1-(21/25+s/5)m+(9/25)s m^2]/[1-(9/25)m]^2.
a_numerator_upper = (
    1 - (F(21, 25) + s_lower / 5) * m_lower
    + F(9, 25) * s_upper * m_upper**2
)
a_denominator_lower = (1 - r * m_upper)**2
a_upper = a_numerator_upper / a_denominator_lower
check("a_interval_positive", a_numerator_upper > 0 and a_denominator_lower > 0)

# B_s(1)-1=(297-375s)/128, so convexity gives this mean upper.
b_upper = 1 + m_upper * (297 - 375 * s_lower) / 128
check("b_chord_coefficient_positive", 297 - 375 * s_upper > 0)

first_square_upper = t**2 * (1 - kappa_lower) * a_upper
second_square_upper = (1-t)**2 * kappa_upper * b_upper
first_bound = F(35317, 100000)
second_bound = F(35391, 100000)
total_bound = first_bound + second_bound
check("first_square_root_bound", first_square_upper < first_bound**2)
check("second_square_root_bound", second_square_upper < second_bound**2)
check("total_is_70708_over_100000", total_bound == F(70708, 100000))
check("strictly_below_inverse_sqrt_two", total_bound**2 < F(1, 2))

# Coarse exact constants used in the analytic curvature/monotonicity proof.
check("a_curvature_endpoint_positive", F(5076, 15625) - F(7020, 15625)*F(2, 3) > 0)
check("b_curvature_endpoint_positive", F(1206, 625) - F(1770, 625)*F(2, 3) > 0)
check("a_along_curve_decreases", F(3, 25) - F(33, 250)*F(2, 3) > 0)
check("upper_function_decreases", F(3, 16) > F(189, 1280))

print(json.dumps({
    "status": "PASS",
    "check_count": len(checks),
    "checks": checks,
    "first_bound": str(first_bound),
    "second_bound": str(second_bound),
    "total_bound": str(total_bound),
    "square_margin_to_half": str(F(1, 2) - total_bound**2),
    "scope": "one fixed eta=3/5 rational certificate, not an optimization",
    "pi_interval": "31415926/10000000 < pi < 31415927/10000000",
    "pi_interval_source_checker_sha256": "d3af3d3bac9ba4d73a7589ba9ed4ff6261fde3263c64d04de36da7f36a1c65d3",
    "pi_interval_source_result_sha256": "fbc10c4760d963f9364dca586cca3d8df5692ab786cd155634a651fac3a62d9d",
}, indent=2, sort_keys=True))
