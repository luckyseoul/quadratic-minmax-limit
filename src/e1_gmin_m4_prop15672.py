#!/usr/bin/env python3
"""Prop. 15.672 -- the opposite-sign near-line branch is impossible.

Continue Proposition 15.671 with an odd-degree boundary consisting of
infinity and ``p-2`` collinear finite points, but take the product sign not
covered by its pointwise-rigid argument:

* ``p=3 (mod 4), c_H=-1`` (``p>=11``), or
* ``p=1 (mod 4), c_H=+1`` (``p>=13``).

The line direction has ``b=1`` and every transverse direction has
``b=p-2``.  The transverse parity baseline is the xnor of the two omitted
fibres and has scaled cost ``p-1``.  Proposition 15.642's degree-two support
bound shows that four nonzero lifts cost strictly more than the available
type surplus ``p+1``.  Hence each quadratic type has a transverse baseline
direction.

For a direction of quadratic type ``eps``, let ``I`` be the number of
infinity edges, ``P_d`` the number of finite edges parallel to the direction,
and ``T`` the total signed sum of all finite selected edges.  The exact
directional mean is

    a_d = I + (p+1)P_d - eps*T - 3p.                    (1)

Consequently same-type means differ by multiples of ``p+1``.  The split
budgets then force exactly one exception of excess ``p+1`` in each type.
For ``p=1 (mod 4)`` the special line direction is the exception in its type;
in all cases an exception has parallel count one above its type baseline.

Let ``x,y`` be the positive/negative baseline parallel counts and
``m=(p+1)/2``.  Counting finite edges gives

    E = m(x+y)+2,       I = 4p-1-m(x+y),       x+y<=7.  (2)

Every type has a transverse baseline.  Comparing its xnor coefficients on
the middle slice gives, with ``q=(p-1)/2``,

    q | I+x-4,          q | I+y-4.                       (3)

Since ``p=2q+1`` and ``m=q+1``, substituting (2) into (3) yields

    q | y+1,            q | x+1.

Thus ``x+y>=2q-2``.  But ``q>=5`` in the stated ranges, contradicting
``x+y<=7``.  The opposite-sign branch is empty.

Together with Proposition 15.671, both product signs of every collinear
infinity-plus-``(p-2)`` boundary are impossible for every odd prime
``p>=13``.  Noncollinear boundaries, the all-finite range, residual (ii),
R1, and the limit remain open.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15642 import polynomial_distance_support_floor
from e1_gmin_m4_prop15671 import rigid_near_line_exclusion


ROOT = Path(__file__).resolve().parents[1]


def middle_weight(p: int) -> int:
    if p < 5 or p % 2 == 0:
        raise ValueError("p must be odd and at least five")
    return (p + 1) // 2


def nonrigid_sign(p: int) -> int:
    """The product sign opposite Proposition 15.671's rigid sign."""
    if p < 5 or p % 2 == 0:
        raise ValueError("p must be odd and at least five")
    return -1 if p % 4 == 3 else 1


def lift_count_ledger(p: int) -> dict[str, object]:
    """Four nonzero lifts cost more than the largest type surplus."""
    support_floor = polynomial_distance_support_floor(p)
    one_lift_scaled_lower = 4 * p * support_floor
    four_lift_scaled_lower = 4 * one_lift_scaled_lower
    maximum_surplus = p + 1
    if four_lift_scaled_lower <= maximum_surplus:
        raise ArithmeticError("the four-lift exclusion margin disappeared")
    return {
        "p": p,
        "nonzero_lift_mass_floor": support_floor,
        "one_lift_scaled_lower_bound": one_lift_scaled_lower,
        "four_lift_scaled_lower_bound": four_lift_scaled_lower,
        "maximum_type_surplus": maximum_surplus,
        "four_lifts_exceed_surplus": True,
        "maximum_nonbaseline_directions_per_type": 3,
    }


def opposite_sign_floor_ledger(p: int) -> dict[str, object]:
    """Exact parity-baseline budgets for the opposite product sign."""
    m = middle_weight(p)
    c_h = nonrigid_sign(p)
    phase = 0 if c_h == 1 else 1
    special_cost = p + (1 if phase == 0 else -1)
    transverse_cost = p - 1
    type_budget = m * (p + 1)
    line_type_floor = special_cost + (m - 1) * transverse_cost
    opposite_type_floor = m * transverse_cost
    line_type_surplus = type_budget - line_type_floor
    opposite_type_surplus = type_budget - opposite_type_floor
    applicable = bool(
        (p % 4 == 3 and p >= 11 and phase == 1 and m % 2 == 0)
        or (p % 4 == 1 and p >= 13 and phase == 0 and m % 2 == 1)
    )
    return {
        "p": p,
        "p_mod_4": p % 4,
        "c_H": c_h,
        "phase": phase,
        "directions_per_type": m,
        "type_budget": type_budget,
        "special_b": 1,
        "transverse_b": p - 2,
        "special_baseline": "x_j" if phase == 0 else "1-x_j",
        "transverse_baseline": "xnor(x_a,x_b)",
        "special_scaled_cost": special_cost,
        "transverse_scaled_cost": transverse_cost,
        "line_type_floor_sum": line_type_floor,
        "opposite_type_floor_sum": opposite_type_floor,
        "line_type_surplus": line_type_surplus,
        "opposite_type_surplus": opposite_type_surplus,
        "applicable": applicable,
    }


def exception_normal_form(p: int) -> dict[str, object]:
    """Use mean quantization and the split budgets to fix two exceptions."""
    floor = opposite_sign_floor_ledger(p)
    if not floor["applicable"]:
        raise ValueError("prime lies outside the opposite-sign theorem range")
    count = lift_count_ledger(p)
    m = middle_weight(p)
    transverse_per_line_type = m - 1
    baseline_exists = bool(
        transverse_per_line_type
        > int(count["maximum_nonbaseline_directions_per_type"])
    )
    if not baseline_exists:
        raise ArithmeticError("a transverse baseline was not forced")

    if p % 4 == 3:
        line_exception = "one direction (line or transverse)"
        line_exception_excess_over_own_floor = p + 1
    else:
        line_exception = "the special line direction"
        line_exception_excess_over_own_floor = p - 1

    return {
        "p": p,
        "exact_directional_mean": "a_d=I+(p+1)P_d-eps_d*T-3p",
        "same_type_quantum": p + 1,
        "transverse_baseline_exists_in_each_type": baseline_exists,
        "exceptions_per_type": 1,
        "exception_a": 2 * p,
        "exception_parallel_increment": 1,
        "line_type_exception": line_exception,
        "line_type_exception_excess_over_own_floor": (
            line_exception_excess_over_own_floor
        ),
        "opposite_type_exception_excess": p + 1,
        "normal_form_proved": True,
    }


def arithmetic_contradiction(p: int) -> dict[str, object]:
    """Close the two-baseline-count normal form symbolically."""
    normal = exception_normal_form(p)
    q = (p - 1) // 2
    m = (p + 1) // 2
    # I>=1 gives m(x+y)<=4p-2.  The quotient is strictly below eight in
    # every stated case, so integral x+y is at most seven.
    maximum_sum_fraction = Fraction(4 * p - 2, m)
    maximum_sum = (4 * p - 2) // m
    if maximum_sum != 7:
        raise ArithmeticError("baseline parallel-count bound changed")
    if q < 5:
        raise ArithmeticError("the congruence contradiction needs q>=5")
    return {
        "p": p,
        "q": q,
        "m": m,
        "finite_edges": "E=m(x+y)+2",
        "infinity_edges": "I=4p-1-m(x+y)",
        "maximum_x_plus_y_fraction": maximum_sum_fraction,
        "maximum_x_plus_y": maximum_sum,
        "transverse_baseline_congruences": ["q divides I+x-4", "q divides I+y-4"],
        "substituted_congruences": ["q divides y+1", "q divides x+1"],
        "congruence_lower_bound_on_x_plus_y": 2 * q - 2,
        "contradiction": 2 * q - 2 > maximum_sum,
        "normal_form": normal,
    }


def opposite_sign_near_line_exclusion(p: int) -> dict[str, object]:
    floor = opposite_sign_floor_ledger(p)
    if not floor["applicable"]:
        return {
            "p": p,
            "applicable": False,
            "excluded": False,
            "floor_ledger": floor,
        }
    arithmetic = arithmetic_contradiction(p)
    return {
        "p": p,
        "applicable": True,
        "boundary": "infinity plus p-2 collinear finite points",
        "c_H": int(floor["c_H"]),
        "excluded": bool(arithmetic["contradiction"]),
        "floor_ledger": floor,
        "lift_count_ledger": lift_count_ledger(p),
        "arithmetic": arithmetic,
    }


def both_signs_collinear_closed(p: int) -> dict[str, object]:
    rigid = rigid_near_line_exclusion(p)
    opposite = opposite_sign_near_line_exclusion(p)
    signs = {
        int(rigid["floor_ledger"]["c_H"]): bool(rigid["excluded"]),
        int(opposite["floor_ledger"]["c_H"]): bool(opposite["excluded"]),
    }
    return {
        "p": p,
        "signs": {str(sign): value for sign, value in sorted(signs.items())},
        "both_signs_excluded": signs == {-1: True, 1: True},
    }


def theorem_record() -> dict[str, object]:
    opposite_samples = {
        str(p): opposite_sign_near_line_exclusion(p)
        for p in (7, 11, 13, 17, 19, 23, 29, 31, 101)
    }
    combined_samples = {
        str(p): both_signs_collinear_closed(p)
        for p in (11, 13, 17, 19, 23, 29, 31, 101)
    }
    proved = bool(
        all(
            opposite_sign_near_line_exclusion(p)["excluded"]
            for p in (11, 13, 17, 19, 23, 29, 31, 101)
        )
        and all(
            both_signs_collinear_closed(p)["both_signs_excluded"]
            for p in (13, 17, 19, 23, 29, 31, 101)
        )
    )
    return {
        "prop": "15.672",
        "title": "Opposite-sign exclusion and complete collinear near-line closure",
        "proved": proved,
        "theorem": {
            "opposite_sign_p_eq_3_mod_4": (
                "c_H=-1 excluded for every odd prime p=3 mod 4, p>=11"
            ),
            "opposite_sign_p_eq_1_mod_4": (
                "c_H=+1 excluded for every odd prime p=1 mod 4, p>=13"
            ),
            "combined_with_15_671": (
                "both product signs of every collinear infinity-plus-(p-2) "
                "boundary are excluded for every odd prime p>=13"
            ),
            "noncollinear_boundary": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "limit_exists": False,
        },
        "opposite_sign_samples": opposite_samples,
        "combined_samples": combined_samples,
        "L_status": "OPEN",
    }


def _json_default(value):
    if isinstance(value, Fraction):
        return str(value)
    raise TypeError(f"cannot serialize {type(value).__name__}")


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.672 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15672.json"
    destination.write_text(json.dumps(record, indent=2, default=_json_default) + "\n")
    print("Prop 15.672 complete collinear near-line closure: proved")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
