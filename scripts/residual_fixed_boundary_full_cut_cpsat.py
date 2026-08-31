#!/usr/bin/env python3
"""Low-memory full-shell cutting-plane model for one p=7 boundary.

The base pseudo-Boolean model contains the exact edge count, distinguished
edge, boundary and product XORs, all affine shell rows, directional means,
common residues, and any forced minimum-slack equalities.  Complete-shell
rows are separated against each candidate and only violated cardinality cuts
are added.  Thus every infeasibility result and every fully audited witness
has the same mathematical force as the all-at-once model without its memory
cost.
"""
from __future__ import annotations

import argparse
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
from p7_no_infinity_unsaturated_cpsat import atomic_write  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402
from residual_fixed_boundary_full_scip import audit_witness, file_sha256  # noqa: E402


def solve_case(
    c_h: int,
    fixed_boundary: tuple[int, ...],
    seconds: float,
    workers: int,
    round_seconds: float,
    cut_batch: int,
    seed: int,
    fixed_means: dict[int, int] | None = None,
) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    boundary = tuple(sorted(int(value) for value in fixed_boundary))
    if c_h not in (-1, 1):
        raise ValueError("c_H must be +/-1")
    if not boundary or len(boundary) % 2 or len(set(boundary)) != len(boundary):
        raise ValueError("fixed boundary must be a nonempty even set")
    if not all(0 <= value < 50 for value in boundary):
        raise ValueError("fixed boundary vertex is outside the p=7 graph")
    if seconds <= 0 or round_seconds <= 0 or cut_batch <= 0:
        raise ValueError("time limits and cut batch must be positive")

    full_data = geometry(7, "full")
    affine_data = geometry(7, "affine")
    C = full_data["C"]
    edges = full_data["edges"]
    signs = full_data["edge_signs"]
    edge_index = {edge: index for index, edge in enumerate(edges)}
    boundary_set = set(boundary)

    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"edge_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == 29)
    model.add(selected[edge_index[(0, 1)]] == 1)

    incident = [[] for _ in range(50)]
    for index, (a, b) in enumerate(edges):
        incident[a].append(selected[index])
        incident[b].append(selected[index])
    for vertex in range(50):
        if vertex in boundary_set:
            model.add_bool_xor(incident[vertex])
        else:
            model.add_bool_xor([~incident[vertex][0], *incident[vertex][1:]])
    negative = [selected[index] for index, sign in enumerate(signs) if int(sign) == -1]
    if c_h == -1:
        model.add_bool_xor(negative)
    else:
        model.add_bool_xor([~negative[0], *negative[1:]])

    means_by_type = {-1: [], 1: []}
    direction_means = []
    direction_rows = []
    direction_models = []
    infinity_value = int(0 in boundary_set)
    for direction_index, direction in enumerate(projective_directions(7)):
        eps, labels = field_direction_data(7, direction)
        counts = [0] * 7
        for vertex in boundary:
            if vertex:
                counts[labels[vertex - 1]] += 1
        B = {index for index, value in enumerate(counts) if value & 1}
        parity_sign = -eps * c_h
        if infinity_value:
            parity_sign *= eps
        if len(B) & 1:
            parity_sign *= -1
        phase = int(parity_sign == -1)
        floor = int(scaled_direction_floor(7, len(B), phase))
        coefficients = []
        for a, b in edges:
            if a == 0:
                coefficient = 1
            elif labels[a - 1] == labels[b - 1]:
                coefficient = 7
            else:
                coefficient = -eps * int(C[a, b])
            coefficients.append(coefficient)
        mean = model.new_int_var(floor, 32, f"scaled_mean_{direction_index}")
        model.add(
            mean
            == sum(
                coefficient * selected[index]
                for index, coefficient in enumerate(coefficients)
            )
            - 21
        )
        model.add_modulo_equality(0, mean, 2)
        means_by_type[int(eps)].append(mean)
        direction_means.append(mean)
        direction_rows.append(
            {
                "direction": list(direction),
                "eps": int(eps),
                "B": sorted(B),
                "b": len(B),
                "phase": phase,
                "floor": floor,
            }
        )
        direction_models.append(
            {"eps": int(eps), "labels": labels, "B": B, "phase": phase, "floor": floor}
        )
    for eps in (-1, 1):
        model.add(sum(means_by_type[eps]) == 32)
        residue = model.new_int_var(0, 7, f"common_residue_{eps}")
        for mean in means_by_type[eps]:
            model.add_modulo_equality(residue, mean, 8)
    fixed_means = dict(fixed_means or {})
    if len(direction_means) != 8:
        raise AssertionError("failed to restore directional mean order")
    for direction_index, value in fixed_means.items():
        if not 0 <= direction_index < 8:
            raise ValueError("fixed mean direction must lie in 0..7")
        if value < direction_rows[direction_index]["floor"] or value > 32:
            raise ValueError("fixed mean is outside its directional bounds")
        model.add(direction_means[direction_index] == int(value))

    affine_constraints = 0
    for eps in (-1, 1):
        normalized = eps * affine_data["features"][eps].astype(np.int8)
        for row in normalized:
            bad = np.flatnonzero(row < 0).tolist()
            model.add(sum(selected[index] for index in bad) <= 13)
            affine_constraints += 1

    saturated_affine_coefficient_equalities = 0
    for eps in (-1, 1):
        records = [row for row in direction_models if row["eps"] == eps]
        if sum(int(row["floor"]) for row in records) != 32:
            continue
        if any(len(row["B"]) not in (0, 2) for row in records):
            continue
        for saturated_index, record in enumerate(records):
            labels = record["labels"]
            B = record["B"]
            phase = int(record["phase"])
            star_counts = []
            for fibre in range(7):
                star_counts.append(
                    sum(
                        selected[edge_index[(0, value + 1)]]
                        for value, label in enumerate(labels)
                        if label == fibre
                    )
                )
            parallel = sum(
                selected[index]
                for index, (a, b) in enumerate(edges)
                if a != 0 and labels[a - 1] == labels[b - 1]
            )
            gauge = model.new_int_var(
                -50,
                50,
                f"saturated_gauge_{eps}_{saturated_index}",
            )
            target_constant = 3 + 2 * phase if len(B) == 0 else 4
            model.add(
                parallel
                == target_constant + 3 * gauge - sum(star_counts)
            )
            saturated_affine_coefficient_equalities += 1
            for s in range(7):
                for t in range(s + 1, 7):
                    signed_cross = sum(
                        eps * int(C[a, b]) * selected[index]
                        for index, (a, b) in enumerate(edges)
                        if a != 0
                        and {labels[a - 1], labels[b - 1]} == {s, t}
                    )
                    target_pair = 0
                    if len(B) == 2 and {s, t} == B:
                        target_pair = 1 if phase else -1
                    model.add(
                        signed_cross
                        == target_pair
                        + gauge
                        - star_counts[s]
                        - star_counts[t]
                    )
                    saturated_affine_coefficient_equalities += 1

    normalized_full = {
        eps: eps * full_data["features"][eps].astype(np.int8) for eps in (-1, 1)
    }
    added: set[tuple[int, int]] = set()
    rounds = []
    final_status = "UNKNOWN"
    chosen_edges = None
    previous_values = None
    while time.time() - started < seconds:
        remaining = seconds - (time.time() - started)
        if previous_values is not None:
            model.clear_hints()
            for variable, value in zip(selected, previous_values):
                model.add_hint(variable, int(value))
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = min(float(round_seconds), remaining)
        solver.parameters.num_search_workers = int(workers)
        solver.parameters.random_seed = int(seed + len(rounds))
        solver.parameters.cp_model_presolve = True
        solver.parameters.symmetry_level = 3
        status = solver.solve(model)
        status_name = solver.status_name(status)
        record = {
            "round": len(rounds) + 1,
            "solver_status": status_name,
            "cuts_before": len(added),
            "conflicts": int(solver.num_conflicts),
            "branches": int(solver.num_branches),
            "wall_time_seconds": float(solver.wall_time),
        }
        rounds.append(record)
        if status == cp_model.INFEASIBLE:
            final_status = "INFEASIBLE"
            break
        if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
            final_status = status_name
            break

        values = np.asarray(
            [solver.value(variable) for variable in selected], dtype=np.int16
        )
        previous_values = values
        violations = []
        minimum_scores = {}
        for eps in (-1, 1):
            scores = normalized_full[eps].astype(np.int16) @ values
            minimum_scores[str(eps)] = int(scores.min())
            for row_index in np.flatnonzero(scores < 3):
                key = (eps, int(row_index))
                if key not in added:
                    violations.append((3 - int(scores[row_index]), key))
        record["minimum_normalized_scores"] = minimum_scores
        record["violated_unadded_rows"] = len(violations)
        if not violations:
            final_status = "FEASIBLE"
            chosen_edges = [list(edge) for edge, value in zip(edges, values) if value]
            break
        violations.sort(reverse=True)
        cuts_added = 0
        for _excess, (eps, row_index) in violations[:cut_batch]:
            row = normalized_full[eps][row_index]
            bad = np.flatnonzero(row < 0).tolist()
            model.add(sum(selected[index] for index in bad) <= 13)
            added.add((eps, row_index))
            cuts_added += 1
        record["cuts_added"] = cuts_added
        if cuts_added == 0:
            final_status = "SEPARATION_STALLED"
            break

    shell_paths = {
        str(eps): Path(f"/tmp/max{'plus' if eps == 1 else 'minus'}_p7.npy")
        for eps in (-1, 1)
    }
    out = {
        "experiment": "residual_fixed_boundary_full_cut_cpsat",
        "status": "exact_low_memory_complete_eigenshell_cutting_plane",
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(boundary),
        "boundary_size": len(boundary),
        "edge_variables": len(selected),
        "direction_rows": direction_rows,
        "type_floor_sums": {
            str(eps): sum(row["floor"] for row in direction_rows if row["eps"] == eps)
            for eps in (-1, 1)
        },
        "affine_score_constraints": affine_constraints,
        "saturated_affine_coefficient_equalities": saturated_affine_coefficient_equalities,
        "complete_shell_rows_available": {
            str(eps): int(len(normalized_full[eps])) for eps in (-1, 1)
        },
        "separated_full_shell_cuts": len(added),
        "rounds": rounds,
        "solver_status": final_status,
        "feasible": final_status == "FEASIBLE",
        "finite_infeasibility_certificate": final_status == "INFEASIBLE",
        "workers": workers,
        "seed": seed,
        "round_seconds": round_seconds,
        "cut_batch": cut_batch,
        "fixed_scaled_means": {
            str(key): value for key, value in sorted(fixed_means.items())
        },
        "shell_file_sha256": {
            eps: file_sha256(path) for eps, path in shell_paths.items()
        },
        "elapsed_seconds": time.time() - started,
    }
    if chosen_edges is not None:
        out["chosen_edges_H"] = chosen_edges
        out["witness_audit"] = audit_witness(full_data, c_h, boundary, chosen_edges)
        if not out["witness_audit"]["valid"]:
            raise AssertionError("cutting-plane witness failed direct full-shell audit")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs="+", required=True)
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--round-seconds", type=float, default=30.0)
    parser.add_argument("--cut-batch", type=int, default=128)
    parser.add_argument("--seed", type=int, default=15708001)
    parser.add_argument(
        "--fixed-mean",
        type=int,
        nargs=2,
        action="append",
        metavar=("DIRECTION", "SCALED_MEAN"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    fixed_means = {}
    for direction, value in args.fixed_mean or []:
        if direction in fixed_means:
            raise ValueError("duplicate fixed mean direction")
        fixed_means[direction] = value
    out = solve_case(
        args.c_h,
        tuple(args.fixed_boundary),
        args.seconds,
        args.workers,
        args.round_seconds,
        args.cut_batch,
        args.seed,
        fixed_means,
    )
    if args.output is not None:
        atomic_write(args.output, out)
    compact = {key: value for key, value in out.items() if key != "chosen_edges_H"}
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
