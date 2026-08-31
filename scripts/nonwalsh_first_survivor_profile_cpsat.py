#!/usr/bin/env python3
"""Exact diagnostic for one quantized first-survivor boundary profile.

This is not a proposition.  It asks whether an all-finite subset of
``F_p^2`` realizes prescribed odd-fibre-count histograms in the phase-zero
and phase-one quadratic direction types.  Translation and field-scalar
symmetry fix points 0 and 1; both ``c_H`` values should be tested because a
nonsquare scalar swaps the quadratic direction types.
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


def parse_profile(text: str) -> dict[int, int]:
    out: dict[int, int] = {}
    for item in text.split(","):
        b_text, count_text = item.split(":", 1)
        out[int(b_text)] = int(count_text)
    return out


def solve(
    p: int,
    boundary_size: int,
    c_h: int,
    phase_zero_profile: dict[int, int],
    phase_one_profile: dict[int, int],
    time_limit: float,
    workers: int,
) -> dict[str, object]:
    from ortools.sat.python import cp_model

    started = time.time()
    m = (p + 1) // 2
    if sum(phase_zero_profile.values()) != m:
        raise ValueError("phase-zero profile must contain m directions")
    if sum(phase_one_profile.values()) != m:
        raise ValueError("phase-one profile must contain m directions")
    model = cp_model.CpModel()
    points = [model.new_bool_var(f"point_{index}") for index in range(p * p)]
    model.add(sum(points) == boundary_size)
    model.add(points[0] == 1)
    model.add(points[1] == 1)

    by_phase: dict[int, list[tuple[object, dict[int, object]]]] = {0: [], 1: []}
    directions = projective_directions(p)
    for direction_index, direction in enumerate(directions):
        eps, labels = field_direction_data(p, direction)
        odd_fibres = []
        for fibre in range(p):
            parity = model.new_bool_var(f"r_{direction_index}_{fibre}")
            fibre_points = [
                points[index] for index, label in enumerate(labels) if label == fibre
            ]
            model.add_modulo_equality(parity, sum(fibre_points), 2)
            odd_fibres.append(parity)
        b_var = model.new_int_var(0, p, f"b_{direction_index}")
        model.add(b_var == sum(odd_fibres))
        phase = 0 if -eps * c_h == 1 else 1
        profile = phase_zero_profile if phase == 0 else phase_one_profile
        indicators = {}
        for b in profile:
            indicator = model.new_bool_var(f"is_{direction_index}_{b}")
            model.add(b_var == b).only_enforce_if(indicator)
            model.add(b_var != b).only_enforce_if(indicator.Not())
            indicators[b] = indicator
        model.add(sum(indicators.values()) == 1)
        by_phase[phase].append((b_var, indicators))

    for phase, profile in ((0, phase_zero_profile), (1, phase_one_profile)):
        if len(by_phase[phase]) != m:
            raise ArithmeticError("quadratic direction split changed")
        for b, target in profile.items():
            model.add(
                sum(indicators[b] for _value, indicators in by_phase[phase])
                == target
            )

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = 15675000 + p + (1 if c_h > 0 else 0)
    solver.parameters.symmetry_level = 2
    status = solver.solve(model)
    result: dict[str, object] = {
        "experiment": "nonwalsh_first_survivor_profile_cpsat",
        "status": "boundary_profile_only_not_an_edge_lift",
        "p": p,
        "boundary_size": boundary_size,
        "c_H": c_h,
        "phase_zero_profile": phase_zero_profile,
        "phase_one_profile": phase_one_profile,
        "fixed_points": [0, 1],
        "solver_status": solver.status_name(status),
        "feasible": status in (cp_model.OPTIMAL, cp_model.FEASIBLE),
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time": solver.wall_time,
        "seconds": time.time() - started,
        "not_a_full_residual_certificate": True,
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        result["boundary"] = [
            index for index, variable in enumerate(points) if solver.value(variable)
        ]
        result["boundary_coordinates"] = [
            [index % p, index // p] for index in result["boundary"]
        ]
        result["direction_rows"] = [
            {
                "direction": list(direction),
                "eps": field_direction_data(p, direction)[0],
                "b": solver.value(by_phase[phase][offset][0]),
                "phase": phase,
            }
            for phase in (0, 1)
            for offset, direction in enumerate(
                [
                    d
                    for d in directions
                    if (0 if -field_direction_data(p, d)[0] * c_h == 1 else 1)
                    == phase
                ]
            )
        ]
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--boundary-size", type=int, required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--phase-zero-profile", required=True)
    parser.add_argument("--phase-one-profile", required=True)
    parser.add_argument("--time-limit", type=float, default=120.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = solve(
        args.p,
        args.boundary_size,
        args.c_h,
        parse_profile(args.phase_zero_profile),
        parse_profile(args.phase_one_profile),
        args.time_limit,
        args.workers,
    )
    rendered = json.dumps(result, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
