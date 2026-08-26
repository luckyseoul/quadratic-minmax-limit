#!/usr/bin/env python3
"""Exact coefficient model for saturated p=7 four-finite boundaries.

This handles exactly the no-infinity boundary profiles whose Proposition
15.632 floor sum is 32 in both quadratic direction types.  Every directional
slack is then minimum-mean.  The b=0 and b=2 slacks are pointwise forced;
phase-zero b=4 is uniquely ``(t-2)^2``; and phase-one b=4 has the complete
36-element catalog from ``p7_size_four_slack_classify.py``.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from p7_size_four_slack_classify import (  # noqa: E402
    classify_four_odd_fibres_phase_one,
)
from residual_boundary_four_lift_cpsat import geometry, verify_witness  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(
        f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp"
    )
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def mapped_phase_one_b4_targets(B: set[int]) -> tuple[tuple[int, ...], ...]:
    """Map the canonical 36 target coefficient vectors to an arbitrary B."""
    if len(B) != 4:
        raise ValueError("B must have four fibres")
    catalog = classify_four_odd_fibres_phase_one()
    if not catalog["proved"]:
        raise AssertionError("phase-one b=4 catalog failed its exact audit")
    canonical_pairs = tuple(itertools.combinations(range(7), 2))
    actual_B = sorted(B)
    actual_complement = sorted(set(range(7)) - B)
    permutation = dict(
        zip(range(4), actual_B)
    ) | dict(zip(range(4, 7), actual_complement))
    actual_pairs = tuple(itertools.combinations(range(7), 2))
    actual_pair_index = {pair: index for index, pair in enumerate(actual_pairs)}
    rows = []
    for item in catalog["catalog"]:
        canonical = item["target_coefficients"]
        constant = int(canonical[0])
        linear = [0] * 7
        pair = [0] * 21
        for s in range(7):
            linear[permutation[s]] = int(canonical[1 + s])
        for index, (s, t) in enumerate(canonical_pairs):
            mapped_pair = tuple(sorted((permutation[s], permutation[t])))
            pair[actual_pair_index[mapped_pair]] = int(canonical[8 + index])
        rows.append((constant, *linear, *pair))
    if len(rows) != 36 or len(set(rows)) != 36:
        raise AssertionError("mapped phase-one catalog lost a target")
    return tuple(rows)


def fixed_target_options(b: int, phase: int, B: set[int]) -> tuple[tuple[int, ...], ...]:
    """Return rows `(constant, seven linear, 21 pair coefficients)`."""
    pairs = tuple(itertools.combinations(range(7), 2))
    linear = [0] * 7
    pair = [0] * len(pairs)
    if b == 0:
        constant = 3 if phase == 0 else 5
    elif b == 2:
        constant = 4
        unique_pair = tuple(sorted(B))
        pair[pairs.index(unique_pair)] = -1 if phase == 0 else 1
    elif b == 4 and phase == 0:
        constant = 5
        for index, endpoints in enumerate(pairs):
            pair[index] = int(set(endpoints) <= B)
    elif b == 4 and phase == 1:
        return mapped_phase_one_b4_targets(B)
    else:
        raise ValueError(f"unsupported saturated target b={b}, phase={phase}")
    return ((constant, *linear, *pair),)


def solve_case(
    c_h: int,
    fixed_boundary: tuple[int, ...],
    seconds: float,
    workers: int,
    seed: int,
) -> dict:
    from ortools.sat.python import cp_model

    if c_h not in (-1, 1):
        raise ValueError("c_h must be +/-1")
    fixed_boundary = tuple(sorted(fixed_boundary))
    if len(fixed_boundary) != 4 or 0 in fixed_boundary:
        raise ValueError("need four finite boundary vertices")
    started = time.time()
    data = geometry(7, "affine")
    C = data["C"]
    edges = data["edges"]
    signs = data["edge_signs"]
    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"edge_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == 29)
    model.add(selected[edges.index((0, 1))] == 1)

    fixed_set = set(fixed_boundary)
    for vertex in range(50):
        incident = [
            selected[j] for j, (a, b) in enumerate(edges) if vertex in (a, b)
        ]
        if vertex in fixed_set:
            model.add_bool_xor(incident)
        else:
            model.add_bool_xor([~incident[0], *incident[1:]])
    negative = [selected[j] for j, sign in enumerate(signs) if sign == -1]
    if c_h == -1:
        model.add_bool_xor(negative)
    else:
        model.add_bool_xor([~negative[0], *negative[1:]])

    infinity_edges = [selected[edges.index((0, u + 1))] for u in range(49)]
    infinity_count = model.new_int_var(0, 29, "infinity_count")
    model.add(infinity_count == sum(infinity_edges))
    pair_order = tuple(itertools.combinations(range(7), 2))
    type_floors = {-1: 0, 1: 0}
    direction_rows = []
    coefficient_constraints = 0
    for d, direction in enumerate(projective_directions(7)):
        eps, labels = field_direction_data(7, direction)
        boundary_counts = [0] * 7
        for vertex in fixed_boundary:
            boundary_counts[labels[vertex - 1]] += 1
        B = {s for s, value in enumerate(boundary_counts) if value & 1}
        b = len(B)
        phase = int(-eps * c_h == -1)
        floor = scaled_direction_floor(7, b, phase)
        type_floors[eps] += floor
        options = fixed_target_options(b, phase, B)
        if len(options) == 1:
            target_constant = options[0][0]
            target_linear = options[0][1:8]
            target_pairs = options[0][8:]
            option_variable = None
        else:
            option_variable = model.new_int_var(0, len(options) - 1, f"option_{d}")
            target_constant = model.new_int_var(-8, 12, f"target_constant_{d}")
            target_linear = [
                model.new_int_var(-8, 12, f"target_linear_{d}_{s}")
                for s in range(7)
            ]
            target_pairs = [
                model.new_int_var(-8, 12, f"target_pair_{d}_{s}_{t}")
                for s, t in pair_order
            ]
            model.add_allowed_assignments(
                [option_variable, target_constant, *target_linear, *target_pairs],
                [(index, *row) for index, row in enumerate(options)],
            )

        star_counts = [
            sum(infinity_edges[u] for u, label in enumerate(labels) if label == s)
            for s in range(7)
        ]
        parallel = sum(
            selected[j]
            for j, (a, edge_b) in enumerate(edges)
            if a != 0 and labels[a - 1] == labels[edge_b - 1]
        )
        k_d = model.new_int_var(-20, 30, f"coefficient_kernel_{d}")
        model.add(
            parallel
            == target_constant + sum(target_linear) + 3 * k_d - infinity_count
        )
        coefficient_constraints += 1
        for pair_index, (s, t) in enumerate(pair_order):
            signed_cross = sum(
                eps * int(C[a, edge_b]) * selected[j]
                for j, (a, edge_b) in enumerate(edges)
                if a != 0 and {labels[a - 1], labels[edge_b - 1]} == {s, t}
            )
            model.add(
                signed_cross
                == target_pairs[pair_index]
                + k_d
                - star_counts[s]
                - star_counts[t]
                + target_linear[s]
                + target_linear[t]
            )
            coefficient_constraints += 1
        direction_rows.append(
            {
                "direction": direction,
                "eps": eps,
                "b": b,
                "phase": phase,
                "floor": floor,
                "target_option_count": len(options),
                "option_variable": option_variable,
            }
        )
    if type_floors != {-1: 32, 1: 32}:
        raise ValueError(f"boundary is not doubly saturated: {type_floors}")

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    out = {
        "experiment": "p7_no_infinity_saturated_cpsat",
        "status": "exact_saturated_coefficient_edge_model",
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(fixed_boundary),
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "type_floor_sums": {str(key): value for key, value in type_floors.items()},
        "coefficient_constraints": coefficient_constraints,
        "phase_one_b4_catalog_size": 36,
        "direction_rows": [
            {key: value for key, value in row.items() if key != "option_variable"}
            for row in direction_rows
        ],
        "workers": workers,
        "seed": seed,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
    }
    if feasible:
        chosen = [
            list(edge) for edge, variable in zip(edges, selected) if solver.value(variable)
        ]
        out["chosen_edges_H"] = chosen
        out["selected_target_options"] = [
            None
            if row["option_variable"] is None
            else solver.value(row["option_variable"])
            for row in direction_rows
        ]
        out["witness_audit"] = verify_witness(
            7, c_h, chosen, 0, fixed_boundary, "affine"
        )
        if not out["witness_audit"]["valid"]:
            raise AssertionError("saturated coefficient witness failed audit")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs=4, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15653001)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve_case(
        args.c_h,
        tuple(args.fixed_boundary),
        args.seconds,
        args.workers,
        args.seed,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
