#!/usr/bin/env python3
"""BDD pseudo-Boolean CNF for one fixed p=7 slack-catalog tuple."""
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

from p7_no_infinity_unsaturated_cpsat import (  # noqa: E402
    atomic_write,
    direction_target_options,
)
from p7_unsaturated_gf2_catalog_filter import direction_scope  # noqa: E402
from p7_unsaturated_slack_catalog import _no_linear_interpolation_data  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry, verify_witness  # noqa: E402


def add_xor_chain(solver, literals: list[int], target: int, top_id: int) -> int:
    if not literals:
        if target:
            solver.add_clause([])
        return top_id
    current = literals[0]
    for literal in literals[1:]:
        top_id += 1
        result = top_id
        solver.add_clause([current, literal, -result])
        solver.add_clause([current, -literal, result])
        solver.add_clause([-current, literal, result])
        solver.add_clause([-current, -literal, -result])
        current = result
    solver.add_clause([current if target else -current])
    return top_id


def solve_case(
    c_h: int,
    boundary: tuple[int, ...],
    elevated: tuple[int, ...],
    catalog_indices: dict[int, int],
    all_points: bool,
    solver_name: str,
    seconds: float,
) -> dict:
    from pysat.pb import EncType, PBEnc
    from pysat.solvers import Solver

    boundary = tuple(sorted(boundary))
    elevated = tuple(sorted(elevated))
    if c_h not in (-1, 1) or len(boundary) != 4 or 0 in boundary:
        raise ValueError("need c_h=+/-1 and four finite boundary vertices")
    if len(set(boundary)) != 4 or len(set(elevated)) != len(elevated):
        raise ValueError("boundary and elevated directions must be distinct")
    started = time.time()
    direction_rows, type_floors = direction_scope(c_h, boundary, elevated)
    chosen_targets = []
    means_by_type = {-1: [], 1: []}
    elevated_set = set(elevated)
    for index, row in enumerate(direction_rows):
        options = direction_target_options(
            len(row["B"]),
            int(row["phase"]),
            set(row["B"]),
            type_floors[int(row["eps"])],
            index in elevated_set,
        )
        if len(options) == 1:
            option_index = int(catalog_indices.get(index, 0))
            if option_index != 0:
                raise ValueError(f"direction {index} has only catalog index zero")
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

    data = geometry(7, "affine")
    C = data["C"]
    edges = data["edges"]
    edge_literals = list(range(1, len(edges) + 1))
    solver = Solver(name=solver_name)
    top_id = len(edges)
    clause_count = 0

    def add_exact(literals: list[int], bound: int) -> None:
        nonlocal top_id, clause_count
        encoded = PBEnc.equals(
            lits=literals,
            bound=int(bound),
            top_id=top_id,
            encoding=EncType.bdd,
        )
        solver.append_formula(encoded.clauses)
        top_id = max(top_id, int(encoded.nv))
        clause_count += len(encoded.clauses)

    add_exact(edge_literals, 29)
    solver.add_clause([edge_literals[edges.index((0, 1))]])
    clause_count += 1
    fixed_set = set(boundary)
    for vertex in range(50):
        incident = [
            edge_literals[index]
            for index, edge in enumerate(edges)
            if vertex in edge
        ]
        old_top = top_id
        top_id = add_xor_chain(solver, incident, int(vertex in fixed_set), top_id)
        clause_count += 4 * (top_id - old_top) + 1
    negative = [
        edge_literals[index]
        for index, sign in enumerate(data["edge_signs"])
        if int(sign) == -1
    ]
    old_top = top_id
    top_id = add_xor_chain(solver, negative, int(c_h == -1), top_id)
    clause_count += 4 * (top_id - old_top) + 1

    points = tuple(itertools.combinations(range(7), 4))
    point_indices = (
        tuple(range(35))
        if all_points
        else tuple(int(value) for value in _no_linear_interpolation_data()["pivot_rows"])
    )
    pair_order = tuple(itertools.combinations(range(7), 2))
    exact_score_constraints = 0
    target_rows = []
    encoding_started = time.time()
    for direction_index, (row, target) in enumerate(
        zip(direction_rows, chosen_targets)
    ):
        eps = int(row["eps"])
        labels = row["labels"]
        scores = []
        for point_index in point_indices:
            X = points[point_index]
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
                clause_count += 1
                continue
            bad_count = (29 - normalized_score) // 2
            bad = []
            for edge_index, (a, endpoint) in enumerate(edges):
                y_a = eps if a == 0 else (1 if labels[a - 1] in X_set else -1)
                y_b = 1 if labels[endpoint - 1] in X_set else -1
                if eps * y_a * y_b * int(C[a, endpoint]) < 0:
                    bad.append(edge_literals[edge_index])
            add_exact(bad, bad_count)
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
                "encoded_normalized_score_support": sorted(set(scores)),
            }
        )
    encoding_seconds = time.time() - encoding_started

    timer = threading.Timer(float(seconds), solver.interrupt)
    timer.start()
    solve_started = time.time()
    try:
        satisfiable = solver.solve_limited(expect_interrupt=True)
    finally:
        timer.cancel()
    solve_seconds = time.time() - solve_started
    status = (
        "SATISFIABLE"
        if satisfiable is True
        else "UNSATISFIABLE"
        if satisfiable is False
        else "UNKNOWN"
    )
    out = {
        "experiment": "p7_no_infinity_fixed_catalog_pbcnf",
        "status": "exact_fixed_catalog_bdd_pseudo_boolean_cnf",
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(boundary),
        "type_floor_sums": {str(key): value for key, value in type_floors.items()},
        "fixed_elevated_directions": list(elevated),
        "solver_name": solver_name,
        "solver_status": status,
        "feasible": satisfiable is True,
        "finite_infeasibility_certificate": satisfiable is False,
        "point_set": "all_35" if all_points else "rank_21_interpolation_pivots",
        "exact_score_constraints": exact_score_constraints,
        "n_edge_variables": len(edges),
        "n_total_variables": top_id,
        "n_clauses": clause_count,
        "target_rows": target_rows,
        "encoding_seconds": encoding_seconds,
        "solve_seconds": solve_seconds,
        "elapsed_seconds": time.time() - started,
        "accumulated_stats": solver.accum_stats(),
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
            7, c_h, chosen, 0, boundary, "affine"
        )
        out["encoded_subset_witness_only"] = not out["witness_audit"]["valid"]
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
    parser.add_argument("--all-points", action="store_true")
    parser.add_argument(
        "--solver",
        choices=("kissat404", "cadical195", "cadical300"),
        default="kissat404",
    )
    parser.add_argument("--seconds", type=float, default=300.0)
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
        args.all_points,
        args.solver,
        args.seconds,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
