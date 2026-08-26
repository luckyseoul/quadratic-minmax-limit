#!/usr/bin/env python3
"""Exact finite-edge solve for one rigid p=7 positive-profile star.

In the last p=7,k0=0 profile, four directions of one quadratic type have
``kd=2`` and the opposite four have ``kd=0``.  All 24 finite edges have the
populated type.  Exact l1 equality turns every signed inter-fibre equation
into an unsigned edge-count equality.  Fixing the five-point infinity star
therefore leaves a compact 0/1 finite-edge incidence model.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from positive_two_point_additive_cpsat import build_geometry  # noqa: E402


def verify_witness(populated_type: int, star: tuple[int, ...], edges: list[list[int]]) -> bool:
    p = 7
    geometry = build_geometry(p)
    chosen = {tuple(edge) for edge in edges}
    if len(chosen) != 24:
        return False
    finite_edges = geometry["finite_edges"]
    signs = geometry["finite_signs"]
    data = geometry["direction_data"]
    direction_of = geometry["edge_directions"]
    indices = [e for e, edge in enumerate(finite_edges) if edge in chosen]
    if len(indices) != 24 or any(signs[e] != populated_type for e in indices):
        return False
    degree = [0] * (p * p)
    for e in indices:
        u, v = finite_edges[e]
        degree[u] += 1
        degree[v] += 1
    if {u for u, value in enumerate(degree) if value & 1} != set(star) ^ {0}:
        return False
    for d, (eps, labels) in enumerate(data):
        kd = 2 if eps == populated_type else 0
        if sum(direction_of[e] == d for e in indices) != 3 * kd:
            return False
        counts = [sum(labels[u] == s for u in star) for s in range(p)]
        special = labels[0]
        for s, t in itertools.combinations(range(p), 2):
            actual = sum(
                {labels[finite_edges[e][0]], labels[finite_edges[e][1]]} == {s, t}
                for e in indices
            )
            expected = (
                2 + int(s == special) + int(t == special) - counts[s] - counts[t]
                if kd == 2
                else counts[s] + counts[t] - int(s == special) - int(t == special)
            )
            if actual != expected:
                return False
    return True


def solve_case(
    populated_type: int,
    star: tuple[int, ...],
    seconds: float,
    workers: int,
    seed: int,
) -> dict:
    from ortools.sat.python import cp_model

    if populated_type not in (-1, 1):
        raise ValueError("populated_type must be +/-1")
    if len(star) != 5 or len(set(star)) != 5 or not all(0 <= u < 49 for u in star):
        raise ValueError("star must contain five distinct points in F_49")
    started = time.time()
    geometry = build_geometry(7)
    finite_edges = geometry["finite_edges"]
    signs = geometry["finite_signs"]
    data = geometry["direction_data"]
    edge_directions = geometry["edge_directions"]
    allowed_indices = [e for e, sign in enumerate(signs) if sign == populated_type]

    model = cp_model.CpModel()
    selected = {
        e: model.new_bool_var(f"edge_{finite_edges[e][0]}_{finite_edges[e][1]}")
        for e in allowed_indices
    }
    model.add(sum(selected.values()) == 24)

    incident_by_vertex = [[] for _ in range(49)]
    for e, variable in selected.items():
        u, v = finite_edges[e]
        incident_by_vertex[u].append(variable)
        incident_by_vertex[v].append(variable)
    finite_boundary = set(star) ^ {0}
    for u, incident in enumerate(incident_by_vertex):
        if u in finite_boundary:
            model.add_bool_xor(incident)
        else:
            model.add_bool_xor([~incident[0], *incident[1:]])

    populated_directions = []
    for d, (eps, labels) in enumerate(data):
        kd = 2 if eps == populated_type else 0
        if kd:
            populated_directions.append(d)
        parallel = [selected[e] for e in allowed_indices if edge_directions[e] == d]
        model.add(sum(parallel) == 3 * kd)
        counts = [sum(labels[u] == s for u in star) for s in range(7)]
        special = labels[0]
        cross_groups = {pair: [] for pair in itertools.combinations(range(7), 2)}
        for e in allowed_indices:
            u, v = finite_edges[e]
            s, t = labels[u], labels[v]
            if s != t:
                cross_groups[tuple(sorted((s, t)))].append(selected[e])
        for (s, t), variables in cross_groups.items():
            required = (
                2 + int(s == special) + int(t == special) - counts[s] - counts[t]
                if kd == 2
                else counts[s] + counts[t] - int(s == special) - int(t == special)
            )
            if required < 0:
                raise ValueError("star violates the proved exact-l1 sign profile")
            model.add(sum(variables) == required)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = seed
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    out = {
        "experiment": "p7_positive_fixed_star_cpsat",
        "status": "exact_finite_edge_model",
        "p": 7,
        "populated_type": populated_type,
        "populated_directions": populated_directions,
        "star": list(star),
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "workers": workers,
        "seed": seed,
    }
    if feasible:
        chosen = [
            list(finite_edges[e]) for e, variable in selected.items() if solver.value(variable)
        ]
        out["finite_edges"] = chosen
        out["witness_verified"] = verify_witness(populated_type, star, chosen)
        if not out["witness_verified"]:
            raise AssertionError("fixed-star witness failed independent verification")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--populated-type", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--star", type=int, nargs=5, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--seed", type=int, default=15655001)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve_case(
        args.populated_type, tuple(sorted(args.star)), args.seconds, args.workers, args.seed
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
