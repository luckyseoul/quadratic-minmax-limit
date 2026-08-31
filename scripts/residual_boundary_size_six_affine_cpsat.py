#!/usr/bin/env python3
"""Search an exact size-six residual lift in all affine eigenshell rows.

This is an exact necessary-condition model for one prime, Paley edge-product
sign, and infinity-boundary bit.  It selects ``H`` directly, fixes the
distinguished edge, imposes ``|H|=4p+1``, all affine Max+/Max- margin-three
inequalities, the six-vertex odd-degree boundary, the product sign, and
Proposition 15.632's split parity-floor budgets.

INFEASIBLE is a finite certificate for the requested branch.  FEASIBLE is
only an affine witness and still has to be tested against the complete
Boolean eigenshells.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from residual_affine_johnson_milp import (  # noqa: E402
    affine_shell,
    feature_rows,
    unique_rows,
)


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def solve(
    p: int,
    c_h: int,
    infinity_value: int,
    seconds: float,
    workers: int,
) -> dict:
    from ortools.sat.python import cp_model

    if p < 3 or p % 2 == 0:
        raise ValueError("p must be an odd prime")
    if c_h not in (-1, 1) or infinity_value not in (0, 1):
        raise ValueError("need c_h in {+-1} and infinity_value in {0,1}")
    started = time.time()
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    y_plus = affine_shell(p, 1, C)
    y_minus = affine_shell(p, -1, C)
    edges, f_plus = feature_rows(y_plus, C)
    edges_minus, f_minus = feature_rows(y_minus, C)
    if edges_minus != edges:
        raise AssertionError("affine shell edge orders differ")
    edge_index = {edge: index for index, edge in enumerate(edges)}
    distinguished = edge_index[(0, 1)]
    h = 4 * p + 1

    model = cp_model.CpModel()
    selected = [
        model.new_bool_var(f"edge_{a}_{b}") for a, b in edges
    ]
    model.add(sum(selected) == h)
    model.add(selected[distinguished] == 1)

    boundary = [model.new_bool_var(f"boundary_{v}") for v in range(len(C))]
    for vertex in range(len(C)):
        incident = [
            selected[index]
            for index, edge in enumerate(edges)
            if vertex in edge
        ]
        model.add_modulo_equality(boundary[vertex], sum(incident), 2)
    model.add(sum(boundary) == 6)
    model.add(boundary[0] == infinity_value)

    negative = [
        selected[index]
        for index, (a, b) in enumerate(edges)
        if int(C[a, b]) == -1
    ]
    negative_parity = model.new_bool_var("negative_edge_parity")
    model.add_modulo_equality(negative_parity, sum(negative), 2)
    model.add(negative_parity == int(c_h == -1))

    type_budget = (p + 1) ** 2 // 2
    costs_by_type: dict[int, list] = {-1: [], 1: []}
    slacks_by_type: dict[int, list] = {-1: [], 1: []}
    direction_records = []
    for direction_index, direction in enumerate(projective_directions(p)):
        eps, labels = field_direction_data(p, direction)
        fibre_parities = []
        for fibre in range(p):
            parity = model.new_bool_var(
                f"fibre_parity_{direction_index}_{fibre}"
            )
            model.add_modulo_equality(
                parity,
                sum(
                    boundary[1 + value]
                    for value, label in enumerate(labels)
                    if label == fibre
                ),
                2,
            )
            fibre_parities.append(parity)
        odd_fibres = model.new_int_var(0, p, f"odd_fibres_{direction_index}")
        model.add(odd_fibres == sum(fibre_parities))
        floor = model.new_int_var(0, 2 * p, f"parity_floor_{direction_index}")
        table = []
        for b in range(p + 1):
            sign = -eps * c_h
            if infinity_value:
                sign *= eps
            if b & 1:
                sign *= -1
            phase = int(sign == -1)
            table.append([b, scaled_direction_floor(p, b, phase)])
        model.add_allowed_assignments([odd_fibres, floor], table)
        directional_contributions = []
        for index, (a, b) in enumerate(edges):
            if a == 0:
                coefficient = 1
            elif labels[a - 1] == labels[b - 1]:
                coefficient = p
            else:
                coefficient = -eps * int(C[a, b])
            directional_contributions.append(coefficient * selected[index])
        slack = model.new_int_var(
            0, type_budget, f"direction_slack_{direction_index}"
        )
        model.add(slack == sum(directional_contributions) - 3 * p)
        model.add(slack >= floor)
        model.add_modulo_equality(0, slack, 2)
        costs_by_type[eps].append(floor)
        slacks_by_type[eps].append(slack)
        parity_sign = -eps * c_h
        if infinity_value:
            parity_sign *= eps
        if (6 - infinity_value) & 1:
            parity_sign *= -1
        phase = int(parity_sign == -1)
        direction_records.append(
            (
                direction,
                eps,
                labels,
                odd_fibres,
                floor,
                slack,
                fibre_parities,
                phase,
            )
        )
    common_residues = {}
    for eps in (-1, 1):
        # For two directions of the same quadratic type, every selected
        # edge has congruent directional coefficient modulo p+1.  A finite
        # edge parallel to a direction of its own Paley type contributes
        # p == -1 (mod p+1), exactly as it does transversely in every other
        # direction of that type; an edge of the opposite type contributes
        # +1 throughout.  Infinity edges also contribute +1 throughout.
        # Hence all scaled directional slacks of one type share a residue.
        residue = model.new_int_var(0, p, f"common_residue_{eps}")
        for slack in slacks_by_type[eps]:
            model.add_modulo_equality(residue, slack, p + 1)
        common_residues[eps] = residue
        model.add(sum(costs_by_type[eps]) <= type_budget)
        model.add(sum(slacks_by_type[eps]) == type_budget)
    if p == 7 and c_h == -1 and infinity_value == 1:
        # Odd b and phase one give floors 6 (b=1,5) or 14 (b=3).
        # The common residue and total 32 force residue six, hence exactly
        # one mean-14 direction and three mean-six directions in each type.
        for eps in (-1, 1):
            model.add(common_residues[eps] == 6)
            for slack in slacks_by_type[eps]:
                model.add_allowed_assignments([slack], [[6], [14]])

    plus_rows = unique_rows(f_plus)
    minus_rows = unique_rows(f_minus)
    if len(plus_rows) != len(f_plus) or len(minus_rows) != len(f_minus):
        raise AssertionError("affine construction unexpectedly has duplicate rows")

    pointwise_parity_equalities = 0
    p7_infinity_unique_slack_rows = 0
    m = (p + 1) // 2
    for shell_eps, shell, features in (
        (1, y_plus, f_plus),
        (-1, y_minus, f_minus),
    ):
        candidate_directions = [
            record for record in direction_records if record[1] == shell_eps
        ]
        for row_index, (y, feature_row) in enumerate(zip(shell, features)):
            matches = []
            for record in candidate_directions:
                labels = record[2]
                fibre_signs = []
                valid = True
                for fibre in range(p):
                    values = {
                        int(y[1 + value])
                        for value, label in enumerate(labels)
                        if label == fibre
                    }
                    if len(values) != 1:
                        valid = False
                        break
                    fibre_signs.append(values.pop())
                if valid and fibre_signs.count(1) == m:
                    matches.append((record, tuple(fibre_signs)))
            if len(matches) != 1:
                raise AssertionError("affine row did not identify one direction")
            record, fibre_signs = matches[0]
            fibre_parities = record[6]
            phase = int(record[7])
            chosen_fibres = [
                fibre for fibre, sign in enumerate(fibre_signs) if sign == 1
            ]
            intersection = model.new_int_var(
                0, 4, f"boundary_intersection_{shell_eps}_{row_index}"
            )
            model.add(
                intersection
                == sum(fibre_parities[fibre] for fibre in chosen_fibres)
            )
            parity = model.new_bool_var(
                f"slack_parity_{shell_eps}_{row_index}"
            )
            model.add_modulo_equality(
                parity,
                intersection + phase,
                2,
            )
            pointwise_slack = model.new_int_var(
                0, (h - 3) // 2, f"pointwise_slack_{shell_eps}_{row_index}"
            )
            normalized_score = shell_eps * sum(
                int(value) * selected[index]
                for index, value in enumerate(feature_row)
            )
            model.add(normalized_score == 3 + 2 * pointwise_slack)
            model.add_modulo_equality(parity, pointwise_slack, 2)
            if p == 7 and c_h == 1 and infinity_value == 1:
                # Here every direction has odd b, phase zero, floor eight,
                # and the four floors saturate each exact type budget.  The
                # complete rank-21 J(7,4) classification at scaled mean
                # eight is unique: parity itself for b=1,5 and
                # (|X cap B|-2)^2 for b=3.
                model.add_allowed_assignments(
                    [record[3], intersection, pointwise_slack],
                    [
                        *([1, t, t & 1] for t in range(0, 2)),
                        *([3, t, (t - 2) ** 2] for t in range(0, 4)),
                        *([5, t, t & 1] for t in range(2, 5)),
                    ],
                )
                p7_infinity_unique_slack_rows += 1
            elif p == 7 and c_h == -1 and infinity_value == 1:
                minimum = model.new_bool_var(
                    f"minimum_phase_one_{shell_eps}_{row_index}"
                )
                model.add(record[5] == 6).only_enforce_if(minimum)
                model.add(record[5] != 6).only_enforce_if(~minimum)
                model.add_allowed_assignments(
                    [record[3], intersection, pointwise_slack],
                    [
                        *([1, t, (t + 1) & 1] for t in range(0, 2)),
                        *([5, t, (t + 1) & 1] for t in range(2, 5)),
                    ],
                ).only_enforce_if(minimum)
                p7_infinity_unique_slack_rows += 1
            pointwise_parity_equalities += 1

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = 15657000 + 10 * p + 2 * infinity_value + int(c_h == 1)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    out = {
        "experiment": "residual_boundary_size_six_affine_cpsat",
        "status": "exact_affine_necessary_condition_not_full_shell",
        "p": p,
        "c_H": c_h,
        "boundary_size": 6,
        "infinity_value": infinity_value,
        "edge_count": h,
        "distinguished_edge": [0, 1],
        "edge_variables": len(selected),
        "plus_affine_rows": len(plus_rows),
        "minus_affine_rows": len(minus_rows),
        "type_budget": type_budget,
        "exact_directional_slack_equalities": len(direction_records),
        "common_type_slack_congruences_modulus": p + 1,
        "pointwise_slack_parity_equalities": pointwise_parity_equalities,
        "p7_infinity_unique_slack_rows": p7_infinity_unique_slack_rows,
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
        "workers": workers,
    }
    if feasible:
        chosen = [
            edge
            for edge, variable in zip(edges, selected)
            if solver.value(variable)
        ]
        chosen_indices = np.fromiter(
            (solver.value(variable) for variable in selected), dtype=np.int8
        )
        degree = [0] * len(C)
        for a, b in chosen:
            degree[a] += 1
            degree[b] += 1
        observed_boundary = [
            vertex for vertex, value in enumerate(degree) if value & 1
        ]
        rows = []
        for (
            direction,
            eps,
            _labels,
            odd_fibres,
            floor,
            slack,
            _fibre_parities,
            phase,
        ) in direction_records:
            rows.append(
                {
                    "direction": list(direction),
                    "eps": eps,
                    "odd_fibres": solver.value(odd_fibres),
                    "floor": solver.value(floor),
                    "slack": solver.value(slack),
                    "phase": phase,
                }
            )
        plus_scores = f_plus @ chosen_indices
        minus_scores = f_minus @ chosen_indices
        product = int(np.prod([int(C[a, b]) for a, b in chosen], dtype=np.int64))
        out.update(
            {
                "boundary": observed_boundary,
                "chosen_edges_H": [list(edge) for edge in chosen],
                "observed_c_H": product,
                "minimum_plus_affine_score": int(plus_scores.min()),
                "maximum_minus_affine_score": int(minus_scores.max()),
                "type_costs": {
                    str(eps): sum(solver.value(value) for value in costs)
                    for eps, costs in costs_by_type.items()
                },
                "type_slacks": {
                    str(eps): sum(solver.value(value) for value in slacks)
                    for eps, slacks in slacks_by_type.items()
                },
                "directions": rows,
                "witness_audit_valid": bool(
                    len(chosen) == h
                    and (0, 1) in chosen
                    and observed_boundary == [
                        index
                        for index, variable in enumerate(boundary)
                        if solver.value(variable)
                    ]
                    and len(observed_boundary) == 6
                    and int(0 in observed_boundary) == infinity_value
                    and product == c_h
                    and int(plus_scores.min()) >= 3
                    and int(minus_scores.max()) <= -3
                ),
            }
        )
        if not out["witness_audit_valid"]:
            raise AssertionError("affine witness failed independent audit")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, default=7)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--infinity", type=int, choices=(0, 1), required=True)
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve(args.p, args.c_h, args.infinity, args.seconds, args.workers)
    print(json.dumps(out, indent=2), flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
