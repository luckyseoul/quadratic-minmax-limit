#!/usr/bin/env python3
"""Exact p=17 first-survivor arc-profile diagnostic.

The remaining ``u_0=0`` ledger has total pair deficit equal to the pair
budget.  Hence every line contains at most two boundary points.  Encoding
that arc condition directly replaces weak modulo-two fibre constraints by
line occupancies in ``{0,1,2}`` and exact secant counts per direction.

This is diagnostic until every normalized case is proved infeasible.  It
does not assert a residual theorem merely from a timeout.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)


PROFILES = {
    "A": {
        0: {7: 6, 0: 3},
        1: {6: 8, 1: 1},
    },
    "B": {
        0: {7: 6, 1: 1, 0: 2},
        1: {6: 8, 0: 1},
    },
}


def solve(case: str, c_h: int, seconds: float, workers: int) -> dict[str, object]:
    from ortools.sat.python import cp_model

    if case not in PROFILES or c_h not in (-1, 1):
        raise ValueError("invalid case or c_H")
    p = 17
    boundary_size = 14
    m = 9
    started = time.time()
    model = cp_model.CpModel()
    points = [model.new_bool_var(f"x_{v}") for v in range(p * p)]
    model.add(sum(points) == boundary_size)

    # Translation and multiplication by a nonzero F_{p^2} scalar normalize
    # any ordered pair of distinct selected points to 0 and 1.  Running both
    # c_H values covers the possible quadratic-type swap.
    model.add(points[0] == 1)
    model.add(points[1] == 1)

    profile = PROFILES[case]
    by_phase: dict[int, list[tuple[object, dict[int, object], tuple[int, int]]]] = {
        0: [],
        1: [],
    }
    normalized_pair_direction: tuple[int, int] | None = None
    directions = projective_directions(p)
    for direction_index, direction in enumerate(directions):
        eps, labels = field_direction_data(p, direction)
        secants = []
        for fibre in range(p):
            on_line = [
                points[v] for v, label in enumerate(labels) if label == fibre
            ]
            occupancy = model.new_int_var(0, 2, f"occ_{direction_index}_{fibre}")
            model.add(occupancy == sum(on_line))
            secant = model.new_bool_var(f"sec_{direction_index}_{fibre}")
            model.add(occupancy == 2).only_enforce_if(secant)
            model.add(occupancy <= 1).only_enforce_if(secant.Not())
            secants.append(secant)

        count = model.new_int_var(0, boundary_size // 2, f"t_{direction_index}")
        model.add(count == sum(secants))
        phase = 0 if -eps * c_h == 1 else 1

        # Choose the normalized pair from one of the six t=7 directions.
        # Scalar multiplication sends its difference to 1, so the direction
        # whose fibres contain both normalized points must itself have seven
        # secants.  Its epsilon is +1, hence it must be phase zero and the
        # normalized orientation has c_H=-1.
        if labels[0] == labels[1]:
            if normalized_pair_direction is not None:
                raise ArithmeticError("the normalized pair has two directions")
            normalized_pair_direction = direction
            model.add(count == 7)
        indicators = {}
        for target in profile[phase]:
            indicator = model.new_bool_var(f"is_{direction_index}_{target}")
            model.add(count == target).only_enforce_if(indicator)
            model.add(count != target).only_enforce_if(indicator.Not())
            indicators[target] = indicator
        model.add(sum(indicators.values()) == 1)
        by_phase[phase].append((count, indicators, direction))

    for phase in (0, 1):
        if len(by_phase[phase]) != m:
            raise ArithmeticError("quadratic direction split changed")
        for secant_count, target_count in profile[phase].items():
            model.add(
                sum(row[1][secant_count] for row in by_phase[phase])
                == target_count
            )

    if normalized_pair_direction is None:
        raise ArithmeticError("normalized pair direction was not found")

    # The profile itself has exactly C(14,2)=91 secants, but retaining this
    # redundant equality materially strengthens propagation.
    model.add(
        sum(row[0] for phase_rows in by_phase.values() for row in phase_rows)
        == boundary_size * (boundary_size - 1) // 2
    )

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = 15677017 + (case == "B") + (c_h > 0)
    solver.parameters.symmetry_level = 2
    status = solver.solve(model)
    result: dict[str, object] = {
        "experiment": "p17_first_survivor_arc_cpsat",
        "case": case,
        "c_H": c_h,
        "profile_by_phase_secant_count": profile,
        "normalization": {
            "selected_points": [0, 1],
            "pair_chosen_in_seven_secant_direction": True,
            "direction": normalized_pair_direction,
            "forced_normalized_c_H": -1,
        },
        "solver_status": solver.status_name(status),
        "feasible": status in (cp_model.OPTIMAL, cp_model.FEASIBLE),
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time": solver.wall_time,
        "seconds": time.time() - started,
        "rigorous_only_if_infeasible_or_witness_rechecked": True,
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        boundary = [v for v, variable in enumerate(points) if solver.value(variable)]
        result["boundary"] = boundary
        result["boundary_coordinates"] = [[v % p, v // p] for v in boundary]
        result["direction_rows"] = [
            {
                "phase": phase,
                "direction": list(direction),
                "secants": solver.value(count),
                "odd_fibres": boundary_size - 2 * solver.value(count),
            }
            for phase, phase_rows in by_phase.items()
            for count, _indicators, direction in phase_rows
        ]
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=sorted(PROFILES), required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = solve(args.case, args.c_h, args.seconds, args.workers)
    rendered = json.dumps(result, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
