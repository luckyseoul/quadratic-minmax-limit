#!/usr/bin/env python3
"""Prop. 15.688 -- sharp integral quadratic-lift floor and p=19 reduction.

Let ``B`` be a nonzero nonnegative integer-valued quadratic on the middle
slice ``J(p,(p+1)/2)``.  If ``H=max B``, the paired-cube operator gives

    4p E[B] >= p-3                         when H=1,
    4p E[B] >= 2(p+1)-4H                  when H>=2.

The second line uses that the mean of an integer-valued quadratic on a
Boolean cube lies in ``(1/4) Z``: support density is at least ``1/4``, and
equality of the mean with ``1/4`` is impossible when the cube contains the
value ``H>=2``.  Combining it with Proposition 15.642's exact stabilizer
bound shows that ``H>=2`` costs at least ``p+1`` for ``p=3 mod 4`` and at
least ``p-1`` for ``p=1 mod 4``.  Hence, sharply,

    4p E[B] >= p-3.

The Boolean quadratic ``(1-x_i)(1-x_j)`` attains equality.

At the live p=19 second all-finite boundary, corrected minimum row arithmetic
leaves phase-zero residues ``0,2,3,4,6,7`` paired with the unique phase-one
residue 9.  Every positive phase-zero residue forces a quotient-zero, b=0
direction with scaled lift mass in ``{4,6,8,12,14}``, all below the new
floor 16.  In residue zero the minimum profiles have inadmissible pair slack 34,
but they are not the whole exact row: completion-bounded enumeration and
the slack congruence leave 143 phase-labelled profiles (75 global shapes),
with slack histogram ``{0:54,4:37,8:25,12:13,16:7,20:4,24:1,28:1,32:1}``.
The endpoint and top-level gates therefore stay open.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15681 import exact_type_rows, second_even_boundary
from e1_gmin_m4_prop15723 import floor_excess_admissible


ROOT = Path(__file__).resolve().parents[1]
P = 19
M = 10
S = 16
PERIOD = 20
PAIR_DEFICIT_BUDGET = S * (S - 1)


def _histogram(profile: tuple[int, ...]) -> dict[int, int]:
    return dict(sorted(Counter(profile).items()))


@lru_cache(maxsize=None)
def _profile_rows(
    phase: int, u: int, deficit_cap: int
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Enumerate every exact p=19 profile within a deficit cap."""
    if phase not in (0, 1) or not 0 <= u < M:
        raise ValueError("phase must be 0/1 and 0<=u<10")
    target = M - u
    options: list[tuple[int, int, int]] = []
    for b in range(0, S + 1, 2):
        floor_value = full_symbolic_floor(P, b, phase)
        for quotient in range(target + 1):
            excess = 2 * u + PERIOD * quotient - floor_value
            if floor_excess_admissible(P, b, phase, excess):
                options.append((quotient, S - b, b))

    infinity = deficit_cap + S * M + 1
    completion = [[infinity] * (target + 1) for _ in range(M + 1)]
    completion[0][0] = 0
    for count in range(1, M + 1):
        for quotient_sum in range(target + 1):
            completion[count][quotient_sum] = min(
                (
                    added + completion[count - 1][quotient_sum - quotient]
                    for quotient, added, _b in options
                    if quotient <= quotient_sum
                ),
                default=infinity,
            )

    states: set[tuple[int, int, tuple[int, ...]]] = {(0, 0, ())}
    for count in range(M):
        next_states: set[tuple[int, int, tuple[int, ...]]] = set()
        for used, deficit, profile in states:
            for quotient, added, b in options:
                new_used = used + quotient
                new_deficit = deficit + added
                remaining_count = M - count - 1
                remaining_sum = target - new_used
                if (
                    new_used <= target
                    and new_deficit <= deficit_cap
                    and new_deficit
                    + completion[remaining_count][remaining_sum]
                    <= deficit_cap
                ):
                    next_states.add(
                        (new_used, new_deficit, tuple(sorted(profile + (b,))))
                    )
        states = next_states
    return tuple(
        sorted(
            (deficit, profile)
            for used, deficit, profile in states
            if used == target
        )
    )


@lru_cache(maxsize=1)
def p19_residue_zero_profiles() -> dict[str, object]:
    """Complete arithmetic census after pair-slack divisibility."""
    phase_zero = _profile_rows(0, 0, PAIR_DEFICIT_BUDGET - 126)
    phase_one = _profile_rows(1, 9, PAIR_DEFICIT_BUDGET - 80)
    candidates = []
    for deficit_zero, profile_zero in phase_zero:
        for deficit_one, profile_one in phase_one:
            slack = PAIR_DEFICIT_BUDGET - deficit_zero - deficit_one
            if slack < 0 or slack % 4:
                continue
            global_profile = Counter(profile_zero)
            global_profile.update(profile_one)
            candidates.append(
                {
                    "phase_profiles_b": {
                        "0": _histogram(profile_zero),
                        "1": _histogram(profile_one),
                    },
                    "phase_deficits": {
                        "0": deficit_zero,
                        "1": deficit_one,
                    },
                    "pair_slack": slack,
                    "global_b_profile": dict(sorted(global_profile.items())),
                    "undetermined_directions": global_profile[S],
                    "arc": slack == 0,
                }
            )

    slack_histogram = dict(sorted(Counter(row["pair_slack"] for row in candidates).items()))
    undetermined_by_slack: dict[int, dict[int, int]] = {}
    for slack in slack_histogram:
        undetermined_by_slack[slack] = dict(
            sorted(
                Counter(
                    int(row["undetermined_directions"])
                    for row in candidates
                    if int(row["pair_slack"]) == slack
                ).items()
            )
        )
    global_shapes = {
        (int(row["pair_slack"]), tuple(row["global_b_profile"].items()))
        for row in candidates
    }
    expected_slack = {0: 54, 4: 37, 8: 25, 12: 13, 16: 7, 20: 4, 24: 1, 28: 1, 32: 1}
    if len(phase_zero) != 60 or len(phase_one) != 9:
        raise ArithmeticError("p=19 exact profile row counts changed")
    if len(candidates) != 143 or slack_histogram != expected_slack:
        raise ArithmeticError("p=19 residue-zero census changed")
    if len(global_shapes) != 75:
        raise ArithmeticError("p=19 global shape count changed")
    return {
        "p": P,
        "boundary_size": S,
        "pair_deficit_budget": PAIR_DEFICIT_BUDGET,
        "phase_zero_row_count": len(phase_zero),
        "phase_one_row_count": len(phase_one),
        "phase_labelled_profile_count": len(candidates),
        "global_shape_count": len(global_shapes),
        "pair_slack_histogram": slack_histogram,
        "undetermined_direction_histogram_by_slack": undetermined_by_slack,
        "profiles": candidates,
        "proved": True,
    }


def sharp_integral_quadratic_lift_floor(p: int) -> dict[str, object]:
    """Return the symbolic sharp floor for odd ``p>=5``.

    The code checks the algebra at the intersection of the paired-cube and
    stabilizer inequalities.  The support and quarter-integrality facts are
    the proof inputs recorded in the module docstring and evidence note.
    """
    if p < 5 or p % 2 == 0:
        raise ValueError("need odd p>=5")

    h_one_floor = p - 3
    if p % 4 == 3:
        stabilizer_coefficient = Fraction(4)
        intersection_h = Fraction(p + 1, 4)
        h_ge_two_floor = p + 1
    else:
        r = (p - 1) // 4
        stabilizer_coefficient = Fraction(4 * r, r + 1)
        intersection_h = Fraction(r + 1)
        h_ge_two_floor = p - 1

    paired_at_intersection = (
        2 * (p + 1) - 4 * intersection_h
    )
    stabilizer_at_intersection = stabilizer_coefficient * intersection_h
    if paired_at_intersection != h_ge_two_floor:
        raise ArithmeticError("paired-cube intersection changed")
    if stabilizer_at_intersection != h_ge_two_floor:
        raise ArithmeticError("stabilizer intersection changed")
    if h_ge_two_floor <= h_one_floor:
        raise ArithmeticError("H>=2 branch no longer separates")

    equality_mass = Fraction(p - 3, 4 * p)
    equality_scaled_mass = 4 * p * equality_mass
    if equality_scaled_mass != h_one_floor:
        raise ArithmeticError("equality example normalization changed")

    return {
        "p": p,
        "paired_cube_identity": "T B(X)=(B(X)+p E[B])/(p+1)",
        "cube_degree_two_support_floor": Fraction(1, 4),
        "integer_quadratic_cube_mean_lattice": "(1/4) Z",
        "H_equals_one_scaled_floor": h_one_floor,
        "H_at_least_two_paired_floor": "2(p+1)-4H",
        "H_at_least_two_stabilizer_coefficient": stabilizer_coefficient,
        "H_at_least_two_intersection": intersection_h,
        "H_at_least_two_scaled_floor": h_ge_two_floor,
        "sharp_scaled_floor": h_one_floor,
        "sharp_mass_floor": equality_mass,
        "equality_example": "B=(1-x_i)(1-x_j)",
        "equality_example_scaled_mass": equality_scaled_mass,
        "equality_rigidity": (
            "B is Boolean and every paired cube through a support point "
            "has minimum degree-two support 1/4"
        ),
        "proved": True,
    }


def p19_second_boundary_reduction() -> dict[str, object]:
    """Apply the sharp floor to the live ``p=19,s=16`` endpoint."""
    p = P
    s = second_even_boundary(p)
    m = (p + 1) // 2
    pair_budget = s * (s - 1)
    phase_zero = exact_type_rows(p, 0)
    phase_one = exact_type_rows(p, 1)

    pair_survivors: list[dict[str, object]] = []
    for zero in phase_zero:
        for one in phase_one:
            required = int(zero["minimum_deficit"]) + int(
                one["minimum_deficit"]
            )
            if required <= pair_budget:
                pair_survivors.append(
                    {
                        "u0": int(zero["u"]),
                        "u1": int(one["u"]),
                        "phase_zero_minimum_deficit": int(
                            zero["minimum_deficit"]
                        ),
                        "phase_one_minimum_deficit": int(
                            one["minimum_deficit"]
                        ),
                        "pair_slack": pair_budget - required,
                        "phase_zero_profile": zero["profile"],
                        "phase_one_profile": one["profile"],
                    }
                )

    expected_residues = [0, 2, 3, 4, 6, 7]
    if [int(row["u0"]) for row in pair_survivors] != expected_residues:
        raise ArithmeticError("p=19 pair-surviving residues changed")
    if {int(row["u1"]) for row in pair_survivors} != {9}:
        raise ArithmeticError("p=19 phase-one row changed")

    least_positive_b_floor = min(
        full_symbolic_floor(p, b, 0) for b in range(2, s + 1, 2)
    )
    sharp_floor = int(
        sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"]
    )
    positive_rows: list[dict[str, object]] = []
    for row in pair_survivors:
        u = int(row["u0"])
        if u == 0:
            continue
        scaled_mean = 2 * u
        quotient_sum = m - u
        forces_quotient_zero = quotient_sum < m
        forces_b_zero = scaled_mean < least_positive_b_floor
        excluded = (
            forces_quotient_zero
            and forces_b_zero
            and scaled_mean < sharp_floor
        )
        positive_rows.append(
            {
                "u0": u,
                "quotient_sum": quotient_sum,
                "forces_quotient_zero": forces_quotient_zero,
                "zero_quotient_scaled_mean": scaled_mean,
                "least_positive_b_floor": least_positive_b_floor,
                "therefore_b_zero": forces_b_zero,
                "sharp_scaled_lift_floor": sharp_floor,
                "excluded": excluded,
            }
        )

    residue_zero_minimum = [row for row in pair_survivors if int(row["u0"]) == 0]
    if len(residue_zero_minimum) != 1:
        raise ArithmeticError("p=19 residue-zero minimum changed")
    if not all(bool(row["excluded"]) for row in positive_rows):
        raise ArithmeticError("a positive p=19 residue survived")
    census = p19_residue_zero_profiles()

    return {
        "proposition": "15.688",
        "p": p,
        "boundary_size": s,
        "pair_deficit_budget": pair_budget,
        "pair_survivors_before_new_floor": pair_survivors,
        "positive_residue_rows": positive_rows,
        "newly_restored_residue": 7,
        "positive_residues_all_excluded": True,
        "residue_zero_minimum_row": residue_zero_minimum[0],
        "residue_zero_minimum_pair_slack": residue_zero_minimum[0]["pair_slack"],
        "residue_zero_minimum_rejected_modulo_four": (
            int(residue_zero_minimum[0]["pair_slack"]) % 4 != 0
        ),
        "residue_zero_exact_census": census,
        "p19_second_all_finite_endpoint_closed": False,
        "remaining_same_boundary_primes": [17, 19, 23],
        "top_level_gates_changed": False,
        "proved": True,
    }


def theorem_sharp_lift_and_p19() -> dict[str, object]:
    """Combined theorem record used by tests and evidence generation."""
    sample_primes = (5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 101)
    floors = {
        str(p): sharp_integral_quadratic_lift_floor(p)
        for p in sample_primes
    }
    p19 = p19_second_boundary_reduction()
    return {
        "proposition": "15.688",
        "sharp_floor_all_odd_p_at_least_five": True,
        "sharp_floor": "4p E[B] >= p-3",
        "sample_floors": floors,
        "p19_reduction": p19,
        "closes_p19_endpoint": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "closes_type_I": False,
        "L_status": "OPEN",
        "proved": True,
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def main() -> None:
    theorem = theorem_sharp_lift_and_p19()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15688.json"
    target.write_text(
        json.dumps(_jsonable(theorem), indent=2, sort_keys=True) + "\n"
    )
    remainder = theorem["p19_reduction"]["residue_zero_exact_census"]
    print(
        "Prop. 15.688: sharp lift floor 4p E[B]>=p-3; "
        "p=19 positive residues excluded; "
        f"exact residue-zero profiles={remainder['phase_labelled_profile_count']}"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
