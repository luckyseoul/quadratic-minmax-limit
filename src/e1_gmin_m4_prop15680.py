#!/usr/bin/env python3
"""Prop. 15.680 -- close the p=37 next all-finite endpoint.

At ``p=37`` the second even all-finite boundary above ``3(p-1)/4`` has
size ``s=30``.  Exact quotient arithmetic leaves phase-zero residues
``u=2,3,4,5``.  Each forces a quotient-zero direction with ``b=0`` and

    2u = 4p E[B],

where ``B`` is a nonzero nonnegative integer-valued quadratic on
``J(37,19)``.  Proposition 15.642 excludes ``u<=4``.  This proposition
excludes the formerly sharp ``u=5`` row.

The new ingredient is an exact small-mass argument.  Stabilizer averaging
forces a hypothetical mass-ten lift to be ``{0,1,2}``-valued.  The exact
degree-two and degree-four slice-distance bounds exclude the value two, so
the lift would be Boolean.  A self-contained paired-cube restriction then
shows that every nonzero Boolean quadratic on ``J(p,(p+1)/2)`` has density
at least ``(p-3)/(4p)``.  At ``p=37`` this is ``17/74``, whereas the lift
has mean ``5/74``.

This closes only the ``p=37,s=30`` endpoint.  The same boundary at
``p=17,19,23,29,31,41``, later all-finite sizes, strict infinity-plus-p,
residual (ii), R1, global QVAR, Type I, and the limit remain open.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15642 import (
    nonbaseline_scaled_cost_floor,
    polynomial_distance_support_floor,
)
from e1_gmin_m4_prop15669 import full_symbolic_floor


ROOT = Path(__file__).resolve().parents[1]
P = 37
M = 19
S = 30
PERIOD = 38


def slice_distance_support_floor(p: int, degree: int) -> Fraction:
    """Lemma 2 of Amireddy et al. for the middle slice."""
    if p < 2 * degree + 1 or p % 2 == 0 or degree < 0:
        raise ValueError("need odd p with degree at most (p-1)/2")
    middle = (p + 1) // 2
    return Fraction(comb(p - 2 * degree, middle - degree), comb(p, middle))


def paired_cube_boolean_quadratic_floor(p: int) -> dict[str, object]:
    """Density floor for a nonzero Boolean quadratic on ``J(p,(p+1)/2)``.

    Put ``p=2m-1``.  Given a middle set ``X``, choose a uniformly random
    leftover point of ``X``, pair the other ``m-1`` points bijectively with
    the complement, and choose one endpoint from every pair.  This gives a
    uniform point in a Boolean ``(m-1)``-cube through ``X``.

    The resulting Markov operator ``T`` satisfies, on the monomial basis
    through degree two,

        T 1 = 1,
        T x_i = 1/2 + rho x_i,
        T x_i x_j = 1/4 + rho x_i x_j,

    where ``rho=1/(2m)=1/(p+1)``.  Equivalently,

        T f = rho f + (1-rho) E[f]

    for every quadratic ``f``.  If Boolean ``f`` is nonzero and ``f(X)=1``,
    its restriction to every paired cube through ``X`` is a nonzero
    degree-two cube polynomial, hence has support density at least ``1/4``.
    Solving ``rho+(1-rho)E[f]>=1/4`` gives the claimed floor.
    """
    if p < 5 or p % 2 == 0:
        raise ValueError("need odd p>=5")
    middle = (p + 1) // 2
    rho = Fraction(1, 2 * middle)
    coordinate_mean = Fraction(middle, p)
    pair_mean = Fraction(middle * (middle - 1), p * (p - 1))
    coordinate_constant = (1 - rho) * coordinate_mean
    pair_constant = (1 - rho) * pair_mean
    floor_value = (Fraction(1, 4) - rho) / (1 - rho)
    closed = Fraction(p - 3, 4 * p)
    if coordinate_constant != Fraction(1, 2):
        raise ArithmeticError("paired-cube coordinate transition changed")
    if pair_constant != Fraction(1, 4):
        raise ArithmeticError("paired-cube pair transition changed")
    if floor_value != closed:
        raise ArithmeticError("paired-cube density simplification changed")
    return {
        "p": p,
        "middle_weight": middle,
        "cube_dimension": middle - 1,
        "rho": rho,
        "transition_on_one": "T(1)=1",
        "transition_on_coordinate": "T(x_i)=1/2+rho*x_i",
        "transition_on_pair": "T(x_i*x_j)=1/4+rho*x_i*x_j",
        "degree_two_operator": "T(f)=rho*f+(1-rho)*E[f]",
        "cube_distance_floor": Fraction(1, 4),
        "boolean_quadratic_density_floor": floor_value,
        "proved": True,
    }


def p37_mass_ten_exclusion() -> dict[str, object]:
    """Exclude a nonnegative integral quadratic with ``4p E[B]=10``."""
    target_mean = Fraction(10, 4 * P)

    # At p=4r+1 with r=9, the stabilizer identity has endpoint weight
    # 9/370.  Applying it at every point where B=h gives h<=2.
    endpoint_weight = Fraction(9, 370)
    max_value = (target_mean / endpoint_weight).numerator // (
        target_mean / endpoint_weight
    ).denominator

    degree_two_floor = polynomial_distance_support_floor(P)
    if degree_two_floor != slice_distance_support_floor(P, 2):
        raise ArithmeticError("degree-two distance floors disagree")
    value_two_density_upper = target_mean - degree_two_floor
    degree_four_floor = slice_distance_support_floor(P, 4)

    boolean = paired_cube_boolean_quadratic_floor(P)
    boolean_floor = Fraction(boolean["boolean_quadratic_density_floor"])
    proved = bool(
        max_value == 2
        and value_two_density_upper == Fraction(2, 1295)
        and degree_four_floor == Fraction(1938, 441595)
        and degree_four_floor > value_two_density_upper
        and boolean_floor == Fraction(17, 74)
        and boolean_floor > target_mean
    )
    return {
        "p": P,
        "middle_weight": M,
        "scaled_mass": 10,
        "target_mean": target_mean,
        "stabilizer_endpoint_weight": endpoint_weight,
        "maximum_point_value": max_value,
        "value_range": [0, 1, 2],
        "degree_two_support_floor": degree_two_floor,
        "value_two_density_upper": value_two_density_upper,
        "degree_four_support_floor_for_B_times_B_minus_1": degree_four_floor,
        "degree_four_gap": degree_four_floor - value_two_density_upper,
        "value_two_excluded": degree_four_floor > value_two_density_upper,
        "therefore_B_is_boolean": True,
        "paired_cube_boolean_floor": boolean_floor,
        "boolean_gap": boolean_floor - target_mean,
        "mass_ten_excluded": proved,
        "paired_cube_certificate": boolean,
        "proved": proved,
    }


def _exact_profile_dp(phase: int) -> list[dict[str, object]]:
    """Independent exact quotient/floor replay at ``p=37,s=30``."""
    if phase not in (0, 1):
        raise ValueError("phase must be zero or one")
    rows = []
    for u in range(M):
        quotient_sum = M - u
        best_by_quotient: dict[int, tuple[int, int]] = {}
        for b in range(0, S + 1, 2):
            floor_value = full_symbolic_floor(P, b, phase)
            for quotient in range(quotient_sum + 1):
                excess = 2 * u + PERIOD * quotient - floor_value
                if excess < 0 or excess == 2:
                    continue
                candidate = (S - b, b)
                old = best_by_quotient.get(quotient)
                if old is None or candidate[0] < old[0]:
                    best_by_quotient[quotient] = candidate

        states: dict[int, tuple[int, tuple[int, ...]]] = {0: (0, ())}
        for _ in range(M):
            next_states: dict[int, tuple[int, tuple[int, ...]]] = {}
            for used, (deficit, profile) in states.items():
                for quotient, (added, b) in best_by_quotient.items():
                    new_used = used + quotient
                    if new_used > quotient_sum:
                        continue
                    candidate = (deficit + added, profile + (b,))
                    old = next_states.get(new_used)
                    if old is None or candidate[0] < old[0]:
                        next_states[new_used] = candidate
            states = next_states
        if quotient_sum in states:
            deficit, profile = states[quotient_sum]
            rows.append(
                {
                    "u": u,
                    "quotient_sum": quotient_sum,
                    "minimum_deficit": deficit,
                    "profile": dict(sorted(Counter(profile).items())),
                }
            )
    return rows


def p37_pair_and_lift_ledger() -> dict[str, object]:
    """Close every exact pair-surviving residue at ``p=37,s=30``."""
    phase_one = _exact_profile_dp(1)
    expected_phase_one = [
        {
            "u": 18,
            "quotient_sum": 1,
            "minimum_deficit": 504,
            "profile": {2: 18, 30: 1},
        }
    ]
    if phase_one != expected_phase_one:
        raise ArithmeticError("p=37 phase-one profile changed")
    phase_zero = _exact_profile_dp(0)
    phase_one_deficit = 504
    pair_budget = S * (S - 1)
    pair_rows = []
    for row in phase_zero:
        required = int(row["minimum_deficit"]) + phase_one_deficit
        pair_rows.append(
            {
                **row,
                "required_total_deficit": required,
                "pair_budget": pair_budget,
                "pair_slack": pair_budget - required,
                "survives_pair_budget": required <= pair_budget,
            }
        )
    survivors = [row for row in pair_rows if row["survives_pair_budget"]]
    surviving_residues = [int(row["u"]) for row in survivors]
    expected_survivors = {
        2: (328, 38, {0: 10, 2: 1, 30: 8}),
        3: (330, 36, {0: 11, 30: 8}),
        4: (358, 8, {0: 11, 2: 1, 30: 7}),
        5: (360, 6, {0: 12, 30: 7}),
    }
    for row in survivors:
        expected = expected_survivors[int(row["u"])]
        if (
            int(row["minimum_deficit"]),
            int(row["pair_slack"]),
            row["profile"],
        ) != expected:
            raise ArithmeticError("p=37 pair-survivor ledger changed")

    positive_b_floor = min(
        full_symbolic_floor(P, b, 0) for b in range(2, S + 1, 2)
    )
    old_lift_floor = nonbaseline_scaled_cost_floor(P)
    mass_ten = p37_mass_ten_exclusion()
    lift_rows = []
    for u in surviving_residues:
        scaled_mean = 2 * u
        lift_rows.append(
            {
                "u": u,
                "quotient_sum": M - u,
                "forces_quotient_zero": M - u < M,
                "scaled_mean_in_zero_quotient_direction": scaled_mean,
                "minimum_positive_b_floor": positive_b_floor,
                "therefore_b_zero": scaled_mean < positive_b_floor,
                "lift_form": "A=2B with B nonzero nonnegative integral quadratic",
                "excluded_by": (
                    "Proposition 15.642"
                    if scaled_mean < old_lift_floor
                    else "mass-ten paired-cube argument"
                ),
                "excluded": (
                    scaled_mean < old_lift_floor
                    or (scaled_mean == 10 and bool(mass_ten["mass_ten_excluded"]))
                ),
            }
        )

    proved = bool(
        surviving_residues == [2, 3, 4, 5]
        and positive_b_floor == PERIOD
        and old_lift_floor == 10
        and all(row["excluded"] for row in lift_rows)
    )
    return {
        "p": P,
        "s": S,
        "phase_one": phase_one,
        "phase_zero": phase_zero,
        "pair_budget": pair_budget,
        "pair_rows": pair_rows,
        "surviving_residues": surviving_residues,
        "old_nonzero_lift_floor": old_lift_floor,
        "lift_rows": lift_rows,
        "mass_ten_exclusion": mass_ten,
        "endpoint_excluded": proved,
        "proved": proved,
    }


def theorem_record() -> dict[str, object]:
    ledger = p37_pair_and_lift_ledger()
    return {
        "prop": "15.680",
        "title": "The p=37 next all-finite endpoint is impossible",
        "proved": bool(ledger["proved"]),
        "theorem": {
            "p": P,
            "boundary_size": S,
            "statement": "the second even all-finite boundary above 3(p-1)/4 is impossible",
            "remaining_smaller_endpoints": [17, 19, 23, 29, 31, 41],
            "remaining_smaller_endpoints_status": "OPEN_AT_THIS_BOUNDARY_SIZE",
            "later_all_finite_sizes": "OPEN",
            "infinity_present_remainder": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "global_QVAR": False,
            "type_I": False,
            "limit_exists": False,
        },
        "paired_cube_boolean_quadratic_lemma": paired_cube_boolean_quadratic_floor(P),
        "p37_pair_and_lift_ledger": ledger,
        "L_status": "OPEN",
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.680 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15680.json"
    destination.write_text(json.dumps(_jsonable(record), indent=2) + "\n")
    print("Prop 15.680 p=37,s=30 next all-finite endpoint: excluded")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
