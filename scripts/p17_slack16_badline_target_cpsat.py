#!/usr/bin/env python3
"""Test the final phase-blind p=17 slack-sixteen boundary profiles.

Slack-sixteen equality restricts a finite sixteen-point boundary to one of
three global affine-line occupancy patterns.  In terms of the numbers of
3- and 4-secants these are

    (n3, n4) = (0, 2), (2, 1), or (4, 0).

This model asks directly whether one such boundary can have one of the seven
remaining directional odd-fibre histograms.  GL(2,17), acting as PGL(2,17)
on directions, is 3-transitive; hence all directions with odd-fibre count
fourteen may be put in canonical positions when there are at most three of
them.  Such a direction has a unique fibre of occupancy two or three.
Translations normalize the exceptional fibres of two canonical directions;
the residual common scalar normalizes the third exceptional fibre to zero or
one.  When fewer than two such directions exist, a boundary point supplies
the remaining translation normalization.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from ortools.sat.python import cp_model

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import field_direction_data, projective_directions


TARGETS = (
    {0: 4, 2: 11, 14: 3},
    {0: 5, 2: 9, 4: 1, 14: 3},
    {0: 6, 2: 8, 6: 1, 14: 3},
    {0: 6, 2: 8, 8: 1, 12: 1, 14: 2},
    {0: 6, 2: 8, 10: 1, 12: 2, 14: 1},
    {0: 6, 2: 8, 10: 2, 14: 2},
    {0: 6, 2: 8, 12: 4},
)
BAD_LINE_PATTERNS = ((0, 2), (2, 1), (4, 0))


def indicator(model: cp_model.CpModel, value, target: int, name: str):
    out = model.new_bool_var(name)
    model.add(value == target).only_enforce_if(out)
    model.add(value != target).only_enforce_if(out.Not())
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=int, choices=range(len(TARGETS)), required=True)
    parser.add_argument(
        "--pattern", type=int, choices=range(len(BAD_LINE_PATTERNS)), required=True
    )
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    target = TARGETS[args.case]
    required_n3, required_n4 = BAD_LINE_PATTERNS[args.pattern]
    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"x_{point}") for point in range(17**2)]
    model.add(sum(selected) == 16)

    odd_counts = []
    occupancies = []
    direction_three_counts = []
    direction_four_counts = []
    for direction_index, direction in enumerate(projective_directions(17)):
        _epsilon, labels = field_direction_data(17, direction)
        odd_fibres = []
        direction_occupancies = []
        direction_three_secants = []
        direction_four_secants = []
        for fibre in range(17):
            occupancy = model.new_int_var(0, 4, f"occ_{direction_index}_{fibre}")
            model.add(
                occupancy
                == sum(
                    selected[point]
                    for point, label in enumerate(labels)
                    if label == fibre
                )
            )
            direction_occupancies.append(occupancy)
            parity = model.new_bool_var(f"parity_{direction_index}_{fibre}")
            model.add_modulo_equality(parity, occupancy, 2)
            odd_fibres.append(parity)
            direction_three_secants.append(
                indicator(
                    model,
                    occupancy,
                    3,
                    f"is_three_{direction_index}_{fibre}",
                )
            )
            direction_four_secants.append(
                indicator(
                    model,
                    occupancy,
                    4,
                    f"is_four_{direction_index}_{fibre}",
                )
            )

        odd_count = model.new_int_var(0, 16, f"odd_count_{direction_index}")
        model.add(odd_count == sum(odd_fibres))
        model.add_allowed_assignments([odd_count], [(value,) for value in target])
        number_three = model.new_int_var(0, 5, f"number_three_{direction_index}")
        number_four = model.new_int_var(0, 4, f"number_four_{direction_index}")
        model.add(number_three == sum(direction_three_secants))
        model.add(number_four == sum(direction_four_secants))
        # If n_j is the number of j-point fibres in this direction, then
        # n_1+n_3=b and n_1+2n_2+3n_3+4n_4=16.  Exposing the resulting
        # tiny table makes the global bad-line budget propagate by direction.
        allowed_direction_types = []
        for b_value in target:
            for n3_value in range(6):
                for n4_value in range(5):
                    n1_value = b_value - n3_value
                    remainder = 16 - n1_value - 3 * n3_value - 4 * n4_value
                    if n1_value < 0 or remainder < 0 or remainder % 2:
                        continue
                    n2_value = remainder // 2
                    if n1_value + n2_value + n3_value + n4_value <= 17:
                        allowed_direction_types.append(
                            (b_value, n3_value, n4_value)
                        )
        model.add_allowed_assignments(
            [odd_count, number_three, number_four], allowed_direction_types
        )
        odd_counts.append(odd_count)
        occupancies.append(direction_occupancies)
        direction_three_counts.append(number_three)
        direction_four_counts.append(number_four)

    model.add(sum(direction_three_counts) == required_n3)
    model.add(sum(direction_four_counts) == required_n4)
    for value, count in target.items():
        model.add(
            sum(
                indicator(model, odd_count, value, f"is_b_{direction}_{value}")
                for direction, odd_count in enumerate(odd_counts)
            )
            == count
        )

    # PGL(2,17) is 3-transitive on the eighteen projective directions.
    b14_count = target.get(14, 0)
    canonical_directions = {
        1: (17,),
        2: (0, 17),
        3: (0, 1, 17),
    }
    for direction in canonical_directions.get(b14_count, ()):
        model.add(odd_counts[direction] == 14)

    if b14_count >= 2:
        # Directions 0 and 17 are x and y.  Independent translations put
        # their unique multi-point fibres on x=0 and y=0.
        model.add(occupancies[0][0] >= 2)
        model.add(occupancies[17][0] >= 2)
        if b14_count == 3:
            # Direction 1 is x+y.  A common nonzero scalar fixes x=0 and
            # y=0 and sends its exceptional label to 1 unless it is 0.
            third_at_zero = model.new_bool_var("third_exceptional_at_zero")
            model.add(occupancies[1][0] >= 2).only_enforce_if(third_at_zero)
            model.add(occupancies[1][1] >= 2).only_enforce_if(third_at_zero.Not())
    elif b14_count == 1:
        # Put the unique exceptional line at y=0, then translate along it so
        # that one of its selected points is the origin.
        model.add(occupancies[17][0] >= 2)
        model.add(selected[0] == 1)
    else:
        model.add(selected[0] == 1)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.seconds
    solver.parameters.num_search_workers = args.workers
    solver.parameters.symmetry_level = 2
    started = time.time()
    status = solver.solve(model)
    result = {
        "case": args.case,
        "target": target,
        "pattern": args.pattern,
        "required_three_secants": required_n3,
        "required_four_secants": required_n4,
        "status": solver.status_name(status),
        "seconds": time.time() - started,
        "branches": solver.num_branches,
        "conflicts": solver.num_conflicts,
        "wall_time": solver.wall_time,
        "workers": args.workers,
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        result["boundary"] = [
            point for point, variable in enumerate(selected) if solver.value(variable)
        ]
        result["directional_odd_counts"] = [
            solver.value(variable) for variable in odd_counts
        ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
