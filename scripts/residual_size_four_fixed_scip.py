#!/usr/bin/env python3
"""SCIP exact edge model for one fixed residual size-four boundary."""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from residual_boundary_four_lift_cpsat import geometry, verify_witness  # noqa: E402


def solve_case(
    p: int,
    c_h: int,
    fixed_boundary: tuple[int, ...],
    shell_mode: str,
    seconds: float,
    workers: int,
) -> dict:
    from pyscipopt import Model, SCIP_PARAMSETTING, quicksum

    if p not in (5, 7) or c_h not in (-1, 1):
        raise ValueError("need p in {5,7} and c_h in {+-1}")
    started = time.time()
    fixed_boundary = tuple(sorted(fixed_boundary))
    data = geometry(p, shell_mode)
    C = data["C"]
    n = int(data["n"])
    edges = data["edges"]
    signs = data["edge_signs"]
    if len(fixed_boundary) != 4 or len(set(fixed_boundary)) != 4:
        raise ValueError("fixed_boundary must have four distinct vertices")
    if not all(0 <= vertex < n for vertex in fixed_boundary):
        raise ValueError("fixed boundary vertex outside graph")

    model = Model(f"residual_size_four_p{p}")
    model.hideOutput(True)
    model.setRealParam("limits/time", float(seconds))
    model.setIntParam("parallel/maxnthreads", int(max(1, workers)))
    model.setHeuristics(SCIP_PARAMSETTING.AGGRESSIVE)
    model.setPresolve(SCIP_PARAMSETTING.AGGRESSIVE)
    selected = [model.addVar(vtype="B", name=f"edge_{a}_{b}") for a, b in edges]
    model.addCons(quicksum(selected) == 4 * p + 1)
    model.addCons(selected[edges.index((0, 1))] == 1)

    fixed_set = set(fixed_boundary)
    for vertex in range(n):
        incident = [
            selected[j] for j, (a, b) in enumerate(edges) if vertex in (a, b)
        ]
        quotient = model.addVar(
            vtype="I", lb=0, ub=(n - 1) // 2, name=f"degree_half_{vertex}"
        )
        model.addCons(quicksum(incident) - 2 * quotient == int(vertex in fixed_set))

    negative = [selected[j] for j, sign in enumerate(signs) if sign == -1]
    sign_quotient = model.addVar(
        vtype="I", lb=0, ub=len(negative) // 2, name="negative_edge_half"
    )
    model.addCons(
        quicksum(negative) - 2 * sign_quotient == int(c_h == -1)
    )

    budget = (p + 1) ** 2 // 2
    half_means_by_type = {-1: [], 1: []}
    infinity_value = int(0 in fixed_boundary)
    for d, direction in enumerate(projective_directions(p)):
        eps, labels = field_direction_data(p, direction)
        counts = [0] * p
        for vertex in fixed_boundary:
            if vertex != 0:
                counts[labels[vertex - 1]] += 1
        b_value = sum(value & 1 for value in counts)
        sign = -eps * c_h
        if infinity_value:
            sign *= eps
        if b_value & 1:
            sign *= -1
        phase = int(sign == -1)
        floor = scaled_direction_floor(p, b_value, phase)
        coefficients = []
        for a, b in edges:
            if a == 0:
                coefficient = 1
            else:
                la, lb = labels[a - 1], labels[b - 1]
                coefficient = p if la == lb else -eps * int(C[a, b])
            coefficients.append(coefficient)
        half_mean = model.addVar(
            vtype="I", lb=floor // 2, ub=p * (4 * p + 1), name=f"half_a_{d}"
        )
        model.addCons(
            2 * half_mean
            == quicksum(
                coefficient * selected[j]
                for j, coefficient in enumerate(coefficients)
            )
            - 3 * p
        )
        half_means_by_type[eps].append(half_mean)
    for eps in (-1, 1):
        model.addCons(2 * quicksum(half_means_by_type[eps]) == budget)

    score_limit = 2 * p - 1
    n_score_constraints = 0
    for eps in (-1, 1):
        for feature in data["features"][eps]:
            bad = np.flatnonzero(eps * feature < 0).tolist()
            model.addCons(quicksum(selected[j] for j in bad) <= score_limit)
            n_score_constraints += 1

    model.optimize()
    status = str(model.getStatus())
    n_solutions = int(model.getNSols())
    feasible = n_solutions > 0
    out = {
        "experiment": "residual_size_four_fixed_scip",
        "status": "exact_fixed_boundary_edge_model",
        "p": p,
        "c_H": c_h,
        "fixed_boundary": list(fixed_boundary),
        "shell_mode": shell_mode,
        "solver_status": status,
        "feasible": feasible,
        "finite_infeasibility_certificate": status == "infeasible",
        "n_edge_variables": len(selected),
        "n_score_constraints": n_score_constraints,
        "n_solutions": n_solutions,
        "nodes": int(model.getNNodes()),
        "gap": float(model.getGap()) if feasible else None,
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
            p,
            c_h,
            chosen,
            infinity_value,
            fixed_boundary,
            shell_mode,
        )
        if not out["witness_audit"]["valid"]:
            raise AssertionError("SCIP witness failed independent audit")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, choices=(5, 7), required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs=4, required=True)
    parser.add_argument("--shell-mode", choices=("affine", "full"), default="affine")
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve_case(
        args.p,
        args.c_h,
        tuple(args.fixed_boundary),
        args.shell_mode,
        args.seconds,
        args.workers,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
