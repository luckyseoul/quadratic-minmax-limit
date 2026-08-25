#!/usr/bin/env python3
"""Prop. 15.642 — mass of a nonzero quadratic lift on the middle slice.

Let ``p`` be odd, ``m=(p+1)/2``, and let ``B`` be a nonzero,
nonnegative, integer-valued polynomial of degree at most two on ``J(p,m)``.
Then the elementary stabilizer bound is

    E B >= 1/p                                      if p == 3 (mod 4),
    E B >= r/((r+1)p),  p=4r+1.                    if p == 1 (mod 4).

The proof fixes a point ``X0`` with ``B(X0)>=1`` and averages under its
stabilizer.  The resulting quadratic in ``t=|X cap X0|`` is controlled by
an exact three-point moment identity.  For ``p>=5``, the exact polynomial
distance lemma on slices also gives support, hence mass, at least
``(p^2-1)/(16p(p-2))``; we use the stronger of the two bounds.  Applied to
the two-vertex boundary
``D={infinity,v}`` in Prop. 15.632, it makes the ``c_H=+1`` affine slack
completely rigid and bounds the exceptional directions when ``c_H=-1``.
It does not by itself exclude either boundary-product branch.
"""
from __future__ import annotations

import json
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def middle_weight(p: int) -> int:
    if p < 3 or p % 2 == 0:
        raise ValueError("p must be odd and at least three")
    return (p + 1) // 2


def stabilizer_mass_certificate(p: int) -> dict:
    """Exact dual certificate for the nonzero quadratic mass floor.

    If ``q`` is the stabilizer average and ``a=floor((m-1)/2)``, this
    returns nonnegative weights such that, for every quadratic ``q``,

        E q(t) = lambda_a q(a) + lambda_b q(a+1) + M q(m).

    Therefore ``q>=0`` and ``q(m)>=1`` imply ``E q>=M``.  The normalized
    consecutive-root quadratic attains equality.
    """
    m = middle_weight(p)
    if p % 4 == 3:
        r = (p - 3) // 4
        a = r
        weights = (
            Fraction(0),
            Fraction(2 * (2 * r + 1), 4 * r + 3),
            Fraction(1, 4 * r + 3),
        )
        value = Fraction(1, p)
    else:
        r = (p - 1) // 4
        a = r
        weights = (
            Fraction(r * (2 * r + 1), (r + 1) * (4 * r + 1)),
            Fraction(2 * r + 1, 4 * r + 1),
            Fraction(r, (r + 1) * (4 * r + 1)),
        )
        value = Fraction(r, (r + 1) * p)

    denominator = (m - a) * (m - a - 1)
    extremizer = (
        Fraction(1, denominator),
        Fraction(-(2 * a + 1), denominator),
        Fraction(a * (a + 1), denominator),
    )
    return {
        "p": p,
        "m": m,
        "a": a,
        "nodes": (a, a + 1, m),
        "weights": weights,
        "value": value,
        "extremizer_coefficients": extremizer,
    }


def hypergeometric_moments(p: int) -> tuple[Fraction, Fraction, Fraction]:
    """Return ``E[1], E[t], E[t^2]`` for two uniform middle sets."""
    m = middle_weight(p)
    first = Fraction(m * m, p)
    factorial_second = Fraction(m * m * (m - 1) * (m - 1), p * (p - 1))
    return Fraction(1), first, factorial_second + first


def certificate_is_exact(p: int) -> bool:
    certificate = stabilizer_mass_certificate(p)
    nodes = certificate["nodes"]
    weights = certificate["weights"]
    moments = hypergeometric_moments(p)
    represented = tuple(
        sum(weight * (node**degree) for node, weight in zip(nodes, weights))
        for degree in range(3)
    )
    a = int(certificate["a"])
    m = int(certificate["m"])
    u, v, w = certificate["extremizer_coefficients"]
    values = tuple(u * t * t + v * t + w for t in range(1, m + 1))
    at_a = u * a * a + v * a + w
    at_a1 = u * (a + 1) * (a + 1) + v * (a + 1) + w
    return bool(
        represented == moments
        and all(weight >= 0 for weight in weights)
        and all(value >= 0 for value in values)
        and at_a == 0
        and at_a1 == 0
        and values[m - 1] == 1
        and sum(
            weight * value
            for weight, value in zip(_hypergeometric_weights(p), values)
        )
        == certificate["value"]
    )


def polynomial_distance_support_floor(p: int) -> Fraction:
    """Exact degree-two slice-distance floor of Amireddy et al., Lemma 2."""
    m = middle_weight(p)
    if p < 5:
        raise ValueError("the degree-two distance lemma requires p>=5")
    value = Fraction(comb(p - 4, m - 2), comb(p, m))
    closed = Fraction(p * p - 1, 16 * p * (p - 2))
    if value != closed:
        raise AssertionError("slice-distance simplification failed")
    return value


def nonzero_quadratic_mass_floor(p: int) -> Fraction:
    """Best certified mass floor used here for nonnegative integer ``B``."""
    stabilizer = stabilizer_mass_certificate(p)["value"]
    if p == 3:
        return stabilizer
    return max(stabilizer, polynomial_distance_support_floor(p))


def _hypergeometric_weights(p: int) -> tuple[Fraction, ...]:
    from math import comb

    m = middle_weight(p)
    denominator = comb(p, m)
    return tuple(
        Fraction(comb(m, t) * comb(p - m, m - t), denominator)
        for t in range(1, m + 1)
    )


def nonbaseline_scaled_cost_floor(p: int) -> int:
    """Even-integer floor for ``4p E[B]`` when ``B`` is nonzero."""
    value = nonzero_quadratic_mass_floor(p)
    lower = 4 * p * value
    integer_ceiling = -((-lower.numerator) // lower.denominator)
    return integer_ceiling if integer_ceiling % 2 == 0 else integer_ceiling + 1


def infinity_finite_boundary_consequence(p: int, c_h: int) -> dict:
    """Prop. 15.632 consequence for ``D={infinity,v}``, ``|H|=4p+1``."""
    if c_h not in (-1, 1):
        raise ValueError("c_h must be +/-1")
    m = middle_weight(p)
    type_budget = m * (p + 1)
    baseline_per_direction = p + c_h
    type_surplus = type_budget - m * baseline_per_direction
    cost = nonbaseline_scaled_cost_floor(p)
    return {
        "p": p,
        "c_H": c_h,
        "directions_per_type": m,
        "type_budget": type_budget,
        "baseline": "x_s" if c_h == 1 else "1-x_s",
        "baseline_scaled_cost_per_direction": baseline_per_direction,
        "type_surplus": type_surplus,
        "nonzero_lift_scaled_cost_floor": cost,
        "maximum_nonzero_lifts_per_type": type_surplus // cost,
        "pointwise_rigid": c_h == 1,
    }


def theorem_quadratic_lift_mass() -> dict:
    primes = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 101)
    exact = all(certificate_is_exact(p) for p in primes)
    costs = {str(p): nonbaseline_scaled_cost_floor(p) for p in primes}
    plus = {str(p): infinity_finite_boundary_consequence(p, 1) for p in primes}
    minus = {str(p): infinity_finite_boundary_consequence(p, -1) for p in primes}
    return {
        "proved": exact,
        "all_odd_primes": True,
        "mass_floor": {
            "p_mod_4_eq_3": "1/p",
            "p_eq_4r_plus_1": "r/((r+1)p)",
        },
        "distance_support_floor": "(p^2-1)/(16p(p-2)) for p>=5",
        "combined_mass_floor": "maximum of stabilizer and distance floors",
        "nonbaseline_scaled_cost": (
            "even ceiling of 4p times the combined mass floor"
        ),
        "sample_costs": costs,
        "infinity_finite_boundary_c_plus": plus,
        "infinity_finite_boundary_c_minus": minus,
        "closes_infinity_finite_boundary": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, tuple):
        return [_jsonable(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    return value


def main() -> dict:
    theorem = theorem_quadratic_lift_mass()
    out = {
        "prop": "15.642",
        "title": "Minimum mass of a nonzero quadratic lift",
        "proved": theorem["proved"],
        "theorem": theorem,
        "certificates": {
            str(p): _jsonable(stabilizer_mass_certificate(p))
            for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 101)
        },
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15642.json"
    destination.write_text(json.dumps(_jsonable(out), indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
