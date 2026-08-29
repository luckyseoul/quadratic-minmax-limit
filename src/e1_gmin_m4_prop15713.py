#!/usr/bin/env python3
"""Prop. 15.713 -- direction-theorem cut for p=7 infinity plus seven.

In the positive-product branch both four-direction types use phase-zero
floors. The exhaustive projection to odd-fibre-count multisets gives 35
profiles per type and 1,217 ordered pairs satisfying the finite-boundary
pair-deficit budget. These are projected b-profiles, not a count of
residue/quotient-labelled arithmetic states.

If at least four of the eight directions have ``b_d=7``, the seven-point
affine boundary determines at most four directions.  Szőnyi's theorem forces
it to be collinear.  Only the two labelled line profiles survive, excluding
208 arithmetic profiles and leaving 1,009 in this branch.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import product
from pathlib import Path

from e1_gmin_m4_prop15632 import scaled_direction_floor


ROOT = Path(__file__).resolve().parents[1]
P = 7
FINITE_BOUNDARY_SIZE = 7
DIRECTIONS_PER_TYPE = 4
TYPE_BUDGET = 32
PAIR_DEFICIT_BUDGET = FINITE_BOUNDARY_SIZE * (FINITE_BOUNDARY_SIZE - 1)
ALLOWED_B = (1, 3, 5, 7)


def _count_profiles(total: int, parts: int) -> list[tuple[int, ...]]:
    return [row for row in product(range(total + 1), repeat=parts) if sum(row) == total]


def p7_positive_infinity_plus_seven_direction_reduction() -> dict[str, object]:
    """Apply Szőnyi's theorem to the exact positive-branch count ledger."""
    floors = {b: scaled_direction_floor(P, b, 0) for b in ALLOWED_B}
    if floors != {1: 8, 3: 8, 5: 8, 7: 0}:
        raise ArithmeticError("p7 positive infinity-plus-seven floors changed")

    one_type = _count_profiles(DIRECTIONS_PER_TYPE, len(ALLOWED_B))
    if len(one_type) != 35:
        raise ArithmeticError("p7 one-type count-profile census changed")
    # Every count profile completes to budget 32: its minimum is 8 times the
    # number of non-b7 directions and the remaining multiple of eight can be
    # assigned as quotient increments.
    if any(
        sum(count * floors[b] for count, b in zip(row, ALLOWED_B)) > TYPE_BUDGET
        or (TYPE_BUDGET - sum(count * floors[b] for count, b in zip(row, ALLOWED_B))) % 8
        for row in one_type
    ):
        raise ArithmeticError("p7 type-budget completion changed")

    ledger = []
    z_histogram: Counter[int] = Counter()
    for left in one_type:
        for right in one_type:
            deficit = sum(
                (left[index] + right[index]) * (FINITE_BOUNDARY_SIZE - b)
                for index, b in enumerate(ALLOWED_B)
            )
            if deficit > PAIR_DEFICIT_BUDGET:
                continue
            undetermined = left[-1] + right[-1]
            record = {
                "type_zero_counts": dict(zip(ALLOWED_B, left)),
                "type_one_counts": dict(zip(ALLOWED_B, right)),
                "pair_deficit": deficit,
                "undetermined_direction_count": undetermined,
            }
            ledger.append(record)
            z_histogram[undetermined] += 1
    expected_z = {0: 217, 1: 300, 2: 280, 3: 210, 4: 126, 5: 56, 6: 21, 7: 6, 8: 1}
    if len(ledger) != 1217 or z_histogram != expected_z:
        raise ArithmeticError("p7 positive count-profile ledger changed")

    collinear_keys = {
        ((1, 1), (7, 3), (7, 4)),
        ((7, 4), (1, 1), (7, 3)),
    }

    def compact(profile: dict[int, int]) -> tuple[tuple[int, int], ...]:
        return tuple((b, count) for b, count in profile.items() if count)

    excluded = []
    survivors = []
    collinear_survivors = []
    for row in ledger:
        if int(row["undetermined_direction_count"]) < 4:
            survivors.append(row)
            continue
        key = compact(row["type_zero_counts"]) + compact(row["type_one_counts"])
        if key in collinear_keys:
            survivors.append(row)
            collinear_survivors.append(row)
        else:
            excluded.append(row)

    excluded_z = Counter(int(row["undetermined_direction_count"]) for row in excluded)
    surviving_z = Counter(int(row["undetermined_direction_count"]) for row in survivors)
    if (
        len(collinear_survivors) != 2
        or len(excluded) != 208
        or len(survivors) != 1009
        or excluded_z != {4: 126, 5: 56, 6: 21, 7: 4, 8: 1}
        or surviving_z != {0: 217, 1: 300, 2: 280, 3: 210, 7: 2}
    ):
        raise ArithmeticError("p7 Szőnyi reduction ledger changed")

    return {
        "proposition": "15.713",
        "p": P,
        "boundary": "infinity plus seven finite points",
        "product_sign": "positive",
        "common_phase": 0,
        "phase_zero_floors": floors,
        "directions_per_type": DIRECTIONS_PER_TYPE,
        "type_budget": TYPE_BUDGET,
        "one_type_profile_count": len(one_type),
        "ordered_profile_count_before_pair_budget": len(one_type) ** 2,
        "projected_b_profile_count_before": len(ledger),
        "undetermined_direction_histogram_before": dict(sorted(z_histogram.items())),
        "szonyi_noncollinear_minimum_directions": 5,
        "four_undetermined_directions_force_collinear": True,
        "labelled_collinear_profiles": [
            {"type_zero": {1: 1, 7: 3}, "type_one": {7: 4}},
            {"type_zero": {7: 4}, "type_one": {1: 1, 7: 3}},
        ],
        "projected_b_profiles_excluded_here": len(excluded),
        "excluded_undetermined_direction_histogram": dict(sorted(excluded_z.items())),
        "projected_b_profile_count_after": len(survivors),
        "remaining_undetermined_direction_histogram": dict(sorted(surviving_z.items())),
        "positive_p7_infinity_plus_seven_closed": False,
        "negative_p7_infinity_plus_seven_changed": False,
        "top_level_gates_changed": False,
        "uses_solver": False,
        "counts_residue_quotient_labelled_states": False,
        "proved_analytically": True,
    }


def main() -> None:
    theorem = p7_positive_infinity_plus_seven_direction_reduction()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15713.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.713: p7 positive infinity-plus-seven projected b-profiles "
        f"{theorem['projected_b_profile_count_before']} -> "
        f"{theorem['projected_b_profile_count_after']}"
    )


if __name__ == "__main__":
    main()
