#!/usr/bin/env python3
"""Native-cardinality SAT for one fully fixed p=7 slack-catalog tuple.

For a chosen catalog row in every non-singleton direction, the normalized
affine score at each of the 280 Johnson points is fixed.  Since H has 29
edges, ``epsilon*S=29-2*bad`` turns each score into an exact native
cardinality constraint.  MiniCard also receives the edge count natively;
boundary and Paley-product parities use linear Tseitin XOR chains.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import threading
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
from p7_no_infinity_unsaturated_cpsat import (  # noqa: E402
    atomic_write,
    direction_target_options,
)
from residual_boundary_four_lift_cpsat import geometry, verify_witness  # noqa: E402
from residual_size_four_fixed_pysat import add_xor_chain  # noqa: E402


def solve_case(
    c_h: int,
    fixed_boundary: tuple[int, ...],
    elevated_directions: tuple[int, ...],
    catalog_indices: dict[int, int],
    seconds: float,
    solver_name: str = "minicard",
    full_score_cuts: bool = False,
) -> dict:
    from pysat.solvers import Solver

    if c_h not in (-1, 1):
        raise ValueError("c_h must be +/-1")
    fixed_boundary = tuple(sorted(fixed_boundary))
    elevated_directions = tuple(sorted(elevated_directions))
    if len(fixed_boundary) != 4 or 0 in fixed_boundary:
        raise ValueError("need four finite boundary vertices")
    if len(set(fixed_boundary)) != 4:
        raise ValueError("boundary vertices must be distinct")
    if len(set(elevated_directions)) != len(elevated_directions):
        raise ValueError("elevated directions must be distinct")

    started = time.time()
    data = geometry(7, "affine")
    C = data["C"]
    edges = data["edges"]
    signs = data["edge_signs"]
    direction_data = []
    type_floors = {-1: 0, 1: 0}
    for direction in projective_directions(7):
        eps, labels = field_direction_data(7, direction)
        counts = [0] * 7
        for vertex in fixed_boundary:
            counts[labels[vertex - 1]] += 1
        B = {index for index, count in enumerate(counts) if count & 1}
        phase = int(-eps * c_h == -1)
        floor = scaled_direction_floor(7, len(B), phase)
        type_floors[eps] += floor
        direction_data.append(
            {
                "direction": direction,
                "eps": eps,
                "labels": labels,
                "B": B,
                "phase": phase,
                "floor": floor,
            }
        )
    if any(value not in (24, 32) for value in type_floors.values()):
        raise ValueError(f"boundary is outside the surviving scope: {type_floors}")
    elevated_set = set(elevated_directions)
    for eps in (-1, 1):
        expected = 1 if type_floors[eps] == 24 else 0
        observed = sum(
            index in elevated_set and int(row["eps"]) == eps
            for index, row in enumerate(direction_data)
        )
        if observed != expected:
            raise ValueError(f"type {eps} needs {expected} elevated directions")

    chosen_targets = []
    means_by_type = {-1: [], 1: []}
    for index, row in enumerate(direction_data):
        options = direction_target_options(
            len(row["B"]),
            int(row["phase"]),
            set(row["B"]),
            type_floors[int(row["eps"])],
            index in elevated_set,
        )
        if len(options) == 1:
            if index in catalog_indices and catalog_indices[index] != 0:
                raise ValueError(f"direction {index} has only catalog index zero")
            option_index = 0
        else:
            if index not in catalog_indices:
                raise ValueError(
                    f"direction {index} has {len(options)} rows and needs an index"
                )
            option_index = int(catalog_indices[index])
            if not 0 <= option_index < len(options):
                raise ValueError(f"catalog index outside direction {index}")
        option = options[option_index]
        means_by_type[int(row["eps"])].append(int(option[0]))
        chosen_targets.append(
            {
                "catalog_index": option_index,
                "catalog_total": len(options),
                "scaled_mean": int(option[0]),
                "constant": int(option[1]),
                "pairs": tuple(int(value) for value in option[2:]),
            }
        )
    if any(sum(means_by_type[eps]) != 32 for eps in (-1, 1)):
        raise AssertionError("chosen target tuple violates an exact type mean")

    edge_literals = list(range(1, len(edges) + 1))
    top_id = len(edges)
    solver = Solver(name=solver_name)
    if not solver.supports_atmost():
        raise ValueError(f"solver {solver_name} lacks native cardinality constraints")
    solver.add_atmost(edge_literals, 29)
    solver.add_atmost([-literal for literal in edge_literals], len(edges) - 29)
    solver.add_clause([edge_literals[edges.index((0, 1))]])
    fixed_set = set(fixed_boundary)
    for vertex in range(50):
        incident = [
            edge_literals[index]
            for index, edge in enumerate(edges)
            if vertex in edge
        ]
        top_id = add_xor_chain(solver, incident, int(vertex in fixed_set), top_id)
    negative = [
        edge_literals[index]
        for index, sign in enumerate(signs)
        if int(sign) == -1
    ]
    top_id = add_xor_chain(solver, negative, int(c_h == -1), top_id)

    pair_order = tuple(itertools.combinations(range(7), 2))
    exact_score_constraints = 0
    immediate_score_range_contradictions = 0
    target_rows = []
    for direction_index, (row, target) in enumerate(
        zip(direction_data, chosen_targets)
    ):
        eps = int(row["eps"])
        labels = row["labels"]
        scores = []
        for X in itertools.combinations(range(7), 4):
            X_set = set(X)
            normalized_score = target["constant"] + sum(
                target["pairs"][pair_index]
                * (1 if ((s in X_set) == (t in X_set)) else -1)
                for pair_index, (s, t) in enumerate(pair_order)
            )
            if normalized_score < 3 or normalized_score % 2 == 0:
                raise AssertionError("catalog target has an invalid normalized score")
            if normalized_score > 29:
                solver.add_clause([])
                immediate_score_range_contradictions += 1
                scores.append(normalized_score)
                continue
            bad_count = (29 - normalized_score) // 2
            bad = []
            for edge_index, (a, endpoint) in enumerate(edges):
                y_a = eps if a == 0 else (1 if labels[a - 1] in X_set else -1)
                y_b = 1 if labels[endpoint - 1] in X_set else -1
                if eps * y_a * y_b * int(C[a, endpoint]) < 0:
                    bad.append(edge_literals[edge_index])
            solver.add_atmost(bad, bad_count)
            solver.add_atmost([-literal for literal in bad], len(bad) - bad_count)
            exact_score_constraints += 1
            scores.append(normalized_score)
        target_rows.append(
            {
                "direction_index": direction_index,
                "direction": list(row["direction"]),
                "eps": eps,
                "b": len(row["B"]),
                "phase": int(row["phase"]),
                "floor": int(row["floor"]),
                "catalog_index": target["catalog_index"],
                "catalog_total": target["catalog_total"],
                "scaled_mean": target["scaled_mean"],
                "normalized_score_support": sorted(set(scores)),
            }
        )

    full_score_constraint_count = 0
    if full_score_cuts:
        full_data = geometry(7, "full")
        for eps in (-1, 1):
            for feature in full_data["features"][eps]:
                bad = [
                    edge_literals[index]
                    for index, value in enumerate(feature)
                    if eps * int(value) < 0
                ]
                solver.add_atmost(bad, 13)
                full_score_constraint_count += 1

    timer = threading.Timer(float(seconds), solver.interrupt)
    timer.start()
    try:
        satisfiable = solver.solve_limited(expect_interrupt=True)
    finally:
        timer.cancel()
    status = (
        "SATISFIABLE"
        if satisfiable is True
        else "UNSATISFIABLE"
        if satisfiable is False
        else "UNKNOWN"
    )
    out = {
        "experiment": "p7_no_infinity_fixed_catalog_minicard",
        "status": "exact_fully_fixed_catalog_native_cardinality_sat",
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(fixed_boundary),
        "type_floor_sums": {str(key): value for key, value in type_floors.items()},
        "fixed_elevated_directions": list(elevated_directions),
        "solver_status": status,
        "solver_name": solver_name,
        "feasible": satisfiable is True,
        "finite_infeasibility_certificate": satisfiable is False,
        "exact_score_constraints": exact_score_constraints,
        "full_score_cuts": full_score_cuts,
        "full_score_constraint_count": full_score_constraint_count,
        "immediate_score_range_contradictions": immediate_score_range_contradictions,
        "n_edge_variables": len(edges),
        "n_total_variables": top_id,
        "target_rows": target_rows,
        "accumulated_stats": solver.accum_stats(),
        "elapsed_seconds": time.time() - started,
    }
    if satisfiable is True:
        assignment = {literal for literal in solver.get_model() if literal > 0}
        chosen = [
            list(edge)
            for edge, literal in zip(edges, edge_literals)
            if literal in assignment
        ]
        out["chosen_edges_H"] = chosen
        out["witness_audit"] = verify_witness(
            7,
            c_h,
            chosen,
            0,
            fixed_boundary,
            "full" if full_score_cuts else "affine",
        )
        if not out["witness_audit"]["valid"]:
            raise AssertionError("MiniCard witness failed independent audit")
    solver.delete()
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs=4, required=True)
    parser.add_argument("--elevated-directions", type=int, nargs="+", required=True)
    parser.add_argument(
        "--catalog-index",
        type=int,
        nargs=2,
        action="append",
        metavar=("DIRECTION", "INDEX"),
    )
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument(
        "--solver",
        choices=("minicard", "gluecard3", "gluecard4"),
        default="minicard",
    )
    parser.add_argument(
        "--full-score-cuts",
        action="store_true",
        help="add every cached full-eigenshell score cardinality bound",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    catalog_indices = {}
    for direction, index in args.catalog_index or []:
        if direction in catalog_indices:
            raise ValueError(f"duplicate catalog index for direction {direction}")
        catalog_indices[direction] = index
    out = solve_case(
        args.c_h,
        tuple(args.fixed_boundary),
        tuple(args.elevated_directions),
        catalog_indices,
        args.seconds,
        args.solver,
        args.full_score_cuts,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
