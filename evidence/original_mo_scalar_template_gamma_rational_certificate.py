#!/usr/bin/env python3
"""Exact rational checks printed in the finite-template Gamma proof.

No floating-point arithmetic, matrix construction, numerical integration,
optimization, or search is used. Written for one remote verification run.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import factorial
from pathlib import Path


def T(m, z):
    return sum(((-1) ** j * z ** (2 * j + 1) / (2 * j + 1)
                for j in range(m + 1)), F(0))


def E(m):
    return sum((F((-1) ** j, 2 ** j * factorial(j))
                for j in range(m + 1)), F(0))


def I(m):
    return sum((F((-1) ** j, 2 ** j * factorial(j) * (2 * j + 1))
                for j in range(m + 1)), F(0))


checks = []


def check(label, lhs, relation, rhs):
    lhs, rhs = F(lhs), F(rhs)
    passed = {"<": lhs < rhs, ">": lhs > rhs, "=": lhs == rhs}[relation]
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


pi_lo, pi_hi = F("3.1415926"), F("3.1415927")
machin_lo = 16 * T(5, F(1, 5)) - 4 * T(2, F(1, 239))
machin_hi = 16 * T(6, F(1, 5)) - 4 * T(1, F(1, 239))
check("Machin lower endpoint", pi_lo, "<", machin_lo)
check("Machin finite interval ordered", machin_lo, "<", machin_hi)
check("Machin upper endpoint", machin_hi, "<", pi_hi)

e_lo, e_hi = F("0.6065306"), F("0.6065307")
check("exponential lower endpoint", e_lo, "<", E(9))
check("exponential finite interval ordered", E(9), "<", E(10))
check("exponential upper endpoint", E(10), "<", e_hi)

i_lo, i_hi = F("0.8556243"), F("0.8556245")
check("integral lower endpoint", i_lo, "<", I(7))
check("integral finite interval ordered", I(7), "<", I(8))
check("integral upper endpoint", I(8), "<", i_hi)

z_lo, z_hi = F("0.398942"), F("0.398943")
check("normal density scale lower bound", z_lo ** 2 * 2 * pi_hi, "<", 1)
check("normal density scale upper bound", z_hi ** 2 * 2 * pi_lo, ">", 1)

phi_lo, phi_hi = F("0.24197"), F("0.24198")
p_lo, p_hi = F("0.68268"), F("0.68270")
check("phi lower product", phi_lo, "<", z_lo * e_lo)
check("phi upper product", z_hi * e_hi, "<", phi_hi)
check("P lower product", p_lo, "<", 2 * z_lo * i_lo)
check("P upper product", 2 * z_hi * i_hi, "<", p_hi)

v_lo = 1 - 2 * phi_hi
k_lo = p_lo - 2 * phi_hi
fourth_hi = 1 + 2 * p_hi - 8 * phi_lo
check("positive lower variance factor", v_lo, ">", 0)
check("positive lower quadratic factor", k_lo, ">", 0)
check("printed lower v", v_lo, "=", F("0.51604"))
check("printed lower k", k_lo, "=", F("0.19872"))
check("printed upper fourth moment", fourth_hi, "=", F("0.42964"))
check("D slope positive", 2 * p_lo ** 2 - 1 + 2 * phi_lo, ">", 0)

d_lo = 5 * p_lo ** 2 - F(7, 2) + 7 * phi_lo
check("printed D lower arithmetic", d_lo, "=", F("0.524049912"))
check("D lower coarse threshold", d_lo, ">", F("0.524"))
r_hi = fourth_hi - v_lo ** 2 - 2 * k_lo ** 2
check("printed R upper arithmetic", r_hi, "=", F("0.0843634416"))
check("R upper coarse threshold", r_hi, "<", F("0.0844"))

check("sqrt2 upper endpoint squared", F("1.415") ** 2, ">", 2)
check("completion fluctuation squared", F("0.0844") * F("0.415") / 3,
      "<", F("0.109") ** 2)
check("contradiction margin arithmetic", F("0.524") - F("0.109"),
      "=", F("0.415"))

print(json.dumps({
    "status": "PASS",
    "arithmetic": "fractions.Fraction only",
    "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    "check_count": len(checks),
    "checks": checks,
}, indent=2, sort_keys=True))
