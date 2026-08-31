#!/usr/bin/env python3
"""SCIP model for one fully fixed p=7 slack-catalog tuple."""
from __future__ import annotations

import argparse
import itertools
import json
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
from p7_no_infinity_unsaturated_cpsat import (  # noqa: E402
    atomic_write,
    direction_target_options,
)
from residual_boundary_four_lift_cpsat import geometry, verify_witness  # noqa: E402


def solve_case(
    c_h: int,
    fixed_boundary: tuple[int, ...],
    elevated_directions: tuple[int, ...],
    catalog_indices: dict[int, int],
    seconds: float,
    workers: int,
) -> dict:
    from pyscipopt import Model, SCIP_PARAMSETTING, quicksum

    fixed_boundary = tuple(sorted(fixed_boundary))
    elevated_directions = tuple(sorted(elevated_directions))
    if c_h not in (-1, 1) or len(fixed_boundary) != 4 or 0 in fixed_boundary:
        raise ValueError("need c_h=+/-1 and four finite boundary vertices")
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

    model = Model("p7_fixed_catalog")
    model.hideOutput(True)
    model.setRealParam("limits/time", float(seconds))
    model.setIntParam("parallel/maxnthreads", int(max(1, workers)))
    model.setPresolve(SCIP_PARAMSETTING.AGGRESSIVE)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    selected = [model.addVar(vtype="B", name=f"edge_{a}_{b}") for a, b in edges]
    model.addCons(quicksum(selected) == 29)
    model.addCons(selected[edges.index((0, 1))] == 1)
    fixed_set = set(fixed_boundary)
    for vertex in range(50):
        incident = [
            selected[index]
            for index, edge in enumerate(edges)
            if vertex in edge
        ]
        quotient = model.addVar(vtype="I", lb=0, ub=24, name=f"degree_half_{vertex}")
        model.addCons(quicksum(incident) - 2 * quotient == int(vertex in fixed_set))
    negative = [
        selected[index]
        for index, sign in enumerate(signs)
        if int(sign) == -1
    ]
    sign_quotient = model.addVar(vtype="I", lb=0, ub=len(negative) // 2, name="sign_half")
    model.addCons(quicksum(negative) - 2 * sign_quotient == int(c_h == -1))

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
                model.addCons(quicksum(selected) <= -1)
                immediate_score_range_contradictions += 1
                scores.append(normalized_score)
                continue
            bad_count = (29 - normalized_score) // 2
            bad = []
            for edge_index, (a, endpoint) in enumerate(edges):
                y_a = eps if a == 0 else (1 if labels[a - 1] in X_set else -1)
                y_b = 1 if labels[endpoint - 1] in X_set else -1
                if eps * y_a * y_b * int(C[a, endpoint]) < 0:
                    bad.append(selected[edge_index])
            model.addCons(quicksum(bad) == bad_count)
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

    model.optimize()
    status = str(model.getStatus())
    n_solutions = int(model.getNSols())
    feasible = n_solutions > 0
    out = {
        "experiment": "p7_no_infinity_fixed_catalog_scip",
        "status": "exact_fully_fixed_catalog_score_equality_mip",
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(fixed_boundary),
        "type_floor_sums": {str(key): value for key, value in type_floors.items()},
        "fixed_elevated_directions": list(elevated_directions),
        "solver_status": status,
        "feasible": feasible,
        "finite_infeasibility_certificate": status == "infeasible",
        "exact_score_constraints": exact_score_constraints,
        "immediate_score_range_contradictions": immediate_score_range_contradictions,
        "n_edge_variables": len(edges),
        "target_rows": target_rows,
        "nodes": int(model.getNNodes()),
        "n_solutions": n_solutions,
        "elapsed_seconds": time.time() - started,
        "workers": workers,
    }
    if feasible:
        solution = model.getBestSol()
        chosen = [
            list(edge)
            for edge, variable in zip(edges, selected)
            if model.getSolVal(solution, variable) > 0.5
        ]
        out["chosen_edges_H"] = chosen
        out["witness_audit"] = verify_witness(
            7, c_h, chosen, 0, fixed_boundary, "affine"
        )
        if not out["witness_audit"]["valid"]:
            raise AssertionError("SCIP witness failed independent audit")
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
    parser.add_argument("--workers", type=int, default=8)
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
        args.workers,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
