#!/usr/bin/env python3
"""Exact full-eigenshell SCIP model for one prescribed even boundary.

The model selects the ``4p+1`` edges of a residual candidate, fixes the
distinguished edge ``(0,1)``, and imposes the prescribed odd-degree boundary,
Paley edge-product sign, every complete Boolean eigenshell score inequality,
and Proposition 15.632's exact directional means and parity floors.

``infeasible`` is a finite fixed-boundary exclusion.  A feasible edge set is
audited directly against the cached complete shells before it is reported.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
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
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_witness(
    data: dict,
    c_h: int,
    boundary: tuple[int, ...],
    chosen_edges: list[list[int]],
) -> dict:
    edges = data["edges"]
    C = data["C"]
    chosen = {tuple(int(value) for value in edge) for edge in chosen_edges}
    selected = np.asarray([int(edge in chosen) for edge in edges], dtype=np.int16)
    degrees = [0] * int(data["n"])
    for a, b in chosen:
        degrees[a] += 1
        degrees[b] += 1
    observed_boundary = tuple(index for index, degree in enumerate(degrees) if degree & 1)
    product = math.prod(int(C[a, b]) for a, b in chosen)
    score_supports = {}
    for eps in (-1, 1):
        scores = data["features"][eps].astype(np.int16) @ selected
        score_supports[str(eps)] = sorted(int(value) for value in np.unique(scores))
    valid = bool(
        len(chosen) == 4 * int(round(math.sqrt(int(data["n"]) - 1))) + 1
        and (0, 1) in chosen
        and observed_boundary == boundary
        and product == c_h
        and min(score_supports["1"]) >= 3
        and max(score_supports["-1"]) <= -3
    )
    return {
        "valid": valid,
        "edge_count": len(chosen),
        "distinguished_edge_present": (0, 1) in chosen,
        "boundary": list(observed_boundary),
        "c_H": product,
        "plus_score_support": score_supports["1"],
        "minus_score_support": score_supports["-1"],
    }


def solve_case(
    p: int,
    c_h: int,
    fixed_boundary: tuple[int, ...],
    seconds: float,
    workers: int,
) -> dict:
    from pyscipopt import Model, SCIP_PARAMSETTING, quicksum

    if p not in (5, 7) or c_h not in (-1, 1):
        raise ValueError("need p in {5,7} and c_h in {+-1}")
    boundary = tuple(sorted(int(value) for value in fixed_boundary))
    n = p * p + 1
    if not boundary or len(boundary) % 2 or len(set(boundary)) != len(boundary):
        raise ValueError("fixed boundary must be a nonempty even set")
    if not all(0 <= vertex < n for vertex in boundary):
        raise ValueError("fixed boundary vertex is outside the graph")

    started = time.time()
    data = geometry(p, "full")
    C = data["C"]
    edges = data["edges"]
    edge_signs = data["edge_signs"]
    edge_index = {edge: index for index, edge in enumerate(edges)}
    boundary_set = set(boundary)
    h = 4 * p + 1

    model = Model(f"residual_fixed_boundary_full_p{p}")
    model.hideOutput(True)
    model.setRealParam("limits/time", float(seconds))
    model.setIntParam("parallel/maxnthreads", max(1, int(workers)))
    model.setPresolve(SCIP_PARAMSETTING.AGGRESSIVE)
    model.setSeparating(SCIP_PARAMSETTING.AGGRESSIVE)
    model.setHeuristics(SCIP_PARAMSETTING.AGGRESSIVE)
    selected = [model.addVar(vtype="B", name=f"edge_{a}_{b}") for a, b in edges]
    model.addCons(quicksum(selected) == h)
    model.addCons(selected[edge_index[(0, 1)]] == 1)

    incident_indices = [[] for _ in range(n)]
    for index, (a, b) in enumerate(edges):
        incident_indices[a].append(index)
        incident_indices[b].append(index)
    for vertex, indices in enumerate(incident_indices):
        degree_half = model.addVar(vtype="I", lb=0, ub=(n - 1) // 2, name=f"degree_half_{vertex}")
        model.addCons(
            quicksum(selected[index] for index in indices) - 2 * degree_half
            == int(vertex in boundary_set)
        )

    negative_indices = [
        index for index, sign in enumerate(edge_signs) if int(sign) == -1
    ]
    negative_half = model.addVar(
        vtype="I", lb=0, ub=h // 2, name="negative_edge_half"
    )
    model.addCons(
        quicksum(selected[index] for index in negative_indices) - 2 * negative_half
        == int(c_h == -1)
    )

    budget = (p + 1) ** 2 // 2
    half_means = {-1: [], 1: []}
    direction_rows = []
    direction_models = []
    infinity_value = int(0 in boundary_set)
    for direction_index, direction in enumerate(projective_directions(p)):
        eps, labels = field_direction_data(p, direction)
        counts = [0] * p
        for vertex in boundary:
            if vertex:
                counts[labels[vertex - 1]] += 1
        odd_fibres = sum(value & 1 for value in counts)
        odd_fibre_set = {index for index, value in enumerate(counts) if value & 1}
        parity_sign = -eps * c_h
        if infinity_value:
            parity_sign *= eps
        if odd_fibres & 1:
            parity_sign *= -1
        phase = int(parity_sign == -1)
        floor = scaled_direction_floor(p, odd_fibres, phase)
        coefficients = []
        for a, b in edges:
            if a == 0:
                coefficient = 1
            elif labels[a - 1] == labels[b - 1]:
                coefficient = p
            else:
                coefficient = -eps * int(C[a, b])
            coefficients.append(coefficient)
        half_mean = model.addVar(
            vtype="I", lb=floor // 2, ub=budget // 2, name=f"half_mean_{direction_index}"
        )
        model.addCons(
            2 * half_mean
            == quicksum(
                coefficient * selected[index]
                for index, coefficient in enumerate(coefficients)
            )
            - 3 * p
        )
        half_means[int(eps)].append(half_mean)
        direction_rows.append(
            {
                "direction": list(direction),
                "eps": int(eps),
                "odd_fibres": odd_fibres,
                "phase": phase,
                "floor": floor,
            }
        )
        direction_models.append(
            {
                "eps": int(eps),
                "labels": labels,
                "B": odd_fibre_set,
                "phase": phase,
                "floor": floor,
            }
        )
    for eps in (-1, 1):
        model.addCons(2 * quicksum(half_means[eps]) == budget)
        residue_half = model.addVar(
            vtype="I", lb=0, ub=(p - 1) // 2, name=f"common_residue_half_{eps}"
        )
        for direction_index, half_mean in enumerate(half_means[eps]):
            quotient = model.addVar(
                vtype="I", lb=0, ub=budget, name=f"common_residue_quotient_{eps}_{direction_index}"
            )
            model.addCons(
                2 * half_mean - 2 * residue_half == (p + 1) * quotient
            )

    # When one type's parity floors consume its entire exact budget, every
    # directional mean attains its floor.  For p=7 and b in {0,2}, the
    # minimizing quadratic equals its parity lower bound at every affine
    # slice point, so every one of the 35 pointwise slacks is forced.  These
    # equalities are redundant consequences, but expose the rigidity to SCIP.
    saturated_affine_equalities = 0
    if p == 7:
        left = np.asarray([a for a, _b in edges], dtype=np.int16)
        right = np.asarray([b for _a, b in edges], dtype=np.int16)
        for eps in (-1, 1):
            records = [row for row in direction_models if row["eps"] == eps]
            if sum(int(row["floor"]) for row in records) != budget:
                continue
            if any(len(row["B"]) not in (0, 2) for row in records):
                continue
            for record in records:
                labels = record["labels"]
                B = record["B"]
                phase = int(record["phase"])
                for chosen_fibres in itertools.combinations(range(p), (p + 1) // 2):
                    chosen_set = set(chosen_fibres)
                    y = np.empty(n, dtype=np.int8)
                    y[0] = eps
                    y[1:] = np.fromiter(
                        (1 if labels[value] in chosen_set else -1 for value in range(p * p)),
                        dtype=np.int8,
                        count=p * p,
                    )
                    t = len(B & chosen_set)
                    if len(B) == 0:
                        slack = phase
                    elif phase == 0:
                        slack = t * (2 - t)
                    else:
                        slack = (t - 1) ** 2
                    coefficients = (
                        eps
                        * y[left].astype(np.int16)
                        * y[right].astype(np.int16)
                        * C[left, right].astype(np.int16)
                    )
                    model.addCons(
                        quicksum(
                            int(coefficient) * selected[index]
                            for index, coefficient in enumerate(coefficients)
                        )
                        == 3 + 2 * slack
                    )
                    saturated_affine_equalities += 1

    shell_constraints = {}
    for eps in (-1, 1):
        normalized = eps * data["features"][eps].astype(np.int8)
        for row in normalized:
            bad = np.flatnonzero(row < 0).tolist()
            model.addCons(quicksum(selected[index] for index in bad) <= (h - 3) // 2)
        shell_constraints[str(eps)] = int(len(normalized))

    model.optimize()
    solver_status = str(model.getStatus())
    feasible = int(model.getNSols()) > 0
    shell_paths = {
        str(eps): Path(f"/tmp/max{'plus' if eps == 1 else 'minus'}_p{p}.npy")
        for eps in (-1, 1)
    }
    out = {
        "experiment": "residual_fixed_boundary_full_scip",
        "status": "exact_fixed_boundary_complete_eigenshell_model",
        "p": p,
        "c_H": c_h,
        "fixed_boundary": list(boundary),
        "boundary_size": len(boundary),
        "edge_count": h,
        "edge_variables": len(selected),
        "direction_rows": direction_rows,
        "type_floor_sums": {
            str(eps): sum(row["floor"] for row in direction_rows if row["eps"] == eps)
            for eps in (-1, 1)
        },
        "complete_shell_constraints": shell_constraints,
        "saturated_affine_slack_equalities": saturated_affine_equalities,
        "shell_file_sha256": {
            eps: file_sha256(path) for eps, path in shell_paths.items()
        },
        "solver_status": solver_status,
        "feasible": feasible,
        "finite_infeasibility_certificate": solver_status == "infeasible",
        "n_solutions": int(model.getNSols()),
        "nodes": int(model.getNNodes()),
        "gap": float(model.getGap()) if feasible else None,
        "workers": workers,
        "elapsed_seconds": time.time() - started,
    }
    if feasible:
        solution = model.getBestSol()
        chosen = [
            list(edge)
            for edge, variable in zip(edges, selected)
            if model.getSolVal(solution, variable) > 0.5
        ]
        out["chosen_edges_H"] = chosen
        out["witness_audit"] = audit_witness(data, c_h, boundary, chosen)
        if not out["witness_audit"]["valid"]:
            raise AssertionError("SCIP witness failed the direct complete-shell audit")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, choices=(5, 7), required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs="+", required=True)
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve_case(
        args.p,
        args.c_h,
        tuple(args.fixed_boundary),
        args.seconds,
        args.workers,
    )
    if args.output is not None:
        atomic_write(args.output, out)
    compact = {key: value for key, value in out.items() if key != "chosen_edges_H"}
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
