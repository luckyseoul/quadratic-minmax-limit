#!/usr/bin/env python3
"""Exact affine edge-lift model for residual boundaries of size four.

The selected variables are the edges of ``H=G union {e}``, with
``|H|=4p+1`` and the distinguished edge ``e=(0,1)`` fixed in ``H``.  The
model imposes every affine Max+ and Max- score inequality, the exact
four-vertex odd-degree boundary, a chosen Paley edge-product ``c_H``, and
Proposition 15.632's exact directional means and parity floors.

``INFEASIBLE`` is a rigorous finite exclusion of the requested branch.
``FEASIBLE`` is only an affine witness and does not satisfy the full
non-affine Max shells automatically.
"""
from __future__ import annotations

import argparse
import functools
import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(
        f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp"
    )
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


@functools.lru_cache(maxsize=None)
def affine_shell(p: int, eigen_sign: int) -> np.ndarray:
    """All one-dimensional affine halfspaces in one Paley eigenshell."""
    if eigen_sign not in (-1, 1):
        raise ValueError("eigen_sign must be +/-1")
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    q = p * p
    m = (p + 1) // 2
    _q, _mul, _add, chi, _frob, _norm, _ia, _ib = field_ctx(p)
    rows = []
    for r, s in projective_directions(p):
        kernel_generator = (s % p) + ((-r) % p) * p
        if chi(kernel_generator) != eigen_sign:
            continue
        labels = np.fromiter(
            ((r * (u % p) + s * (u // p)) % p for u in range(q)),
            dtype=np.int16,
            count=q,
        )
        for chosen in itertools.combinations(range(p), m):
            selected = np.zeros(p, dtype=np.int8)
            selected[list(chosen)] = 1
            y = np.empty(q + 1, dtype=np.int8)
            y[0] = eigen_sign
            y[1:] = 2 * selected[labels] - 1
            rows.append(y)
    Y = np.unique(np.stack(rows), axis=0)
    if not np.array_equal(
        Y.astype(np.int16) @ C.astype(np.int16),
        eigen_sign * p * Y.astype(np.int16),
    ):
        raise AssertionError("affine rows failed their exact eigenshell check")
    return Y


@functools.lru_cache(maxsize=None)
def geometry(p: int, shell_mode: str = "affine") -> dict:
    if shell_mode not in ("affine", "full"):
        raise ValueError("shell_mode must be affine or full")
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    n = p * p + 1
    edges = [(a, b) for a in range(n) for b in range(a + 1, n)]
    left = np.fromiter((a for a, _b in edges), dtype=np.int32)
    right = np.fromiter((b for _a, b in edges), dtype=np.int32)
    signs = C[left, right].astype(np.int8)
    shells = {}
    features = {}
    for eps in (-1, 1):
        if shell_mode == "affine":
            Y = affine_shell(p, eps)
        else:
            path = Path(f"/tmp/max{'plus' if eps == 1 else 'minus'}_p{p}.npy")
            if not path.exists():
                raise FileNotFoundError(path)
            Y = np.rint(np.load(path)).astype(np.int8)
            if Y.shape[1] != n:
                raise ValueError(f"unexpected shell shape in {path}: {Y.shape}")
            if not np.array_equal(
                Y.astype(np.int16) @ C.astype(np.int16),
                eps * p * Y.astype(np.int16),
            ):
                raise AssertionError(f"cached shell {path} failed eigenshell audit")
        F = (Y[:, left] * Y[:, right] * signs).astype(np.int8)
        shells[eps] = Y
        features[eps] = np.unique(np.ascontiguousarray(F), axis=0)
    return {
        "C": C,
        "n": n,
        "edges": edges,
        "edge_signs": signs,
        "shells": shells,
        "features": features,
        "shell_mode": shell_mode,
    }


def verify_witness(
    p: int,
    c_h: int,
    chosen_edges: list[list[int]],
    requested_infinity: int | None,
    fixed_boundary: tuple[int, ...] | None,
    shell_mode: str,
) -> dict:
    data = geometry(p, shell_mode)
    C = data["C"]
    edges = data["edges"]
    chosen = {tuple(edge) for edge in chosen_edges}
    degree = [0] * data["n"]
    for a, b in chosen:
        degree[a] += 1
        degree[b] += 1
    boundary = tuple(v for v, value in enumerate(degree) if value & 1)
    product = math.prod(int(C[a, b]) for a, b in chosen)
    scores = {
        str(eps): [
            int(sum(int(row[j]) for j, edge in enumerate(edges) if edge in chosen))
            for row in data["features"][eps]
        ]
        for eps in (-1, 1)
    }
    valid = bool(
        len(chosen) == 4 * p + 1
        and (0, 1) in chosen
        and len(boundary) == 4
        and product == c_h
        and (requested_infinity is None or int(0 in boundary) == requested_infinity)
        and (fixed_boundary is None or boundary == fixed_boundary)
        and min(scores["1"]) >= 3
        and max(scores["-1"]) <= -3
    )
    return {
        "valid": valid,
        "boundary": list(boundary),
        "c_H": product,
        "plus_score_support": sorted(set(scores["1"])),
        "minus_score_support": sorted(set(scores["-1"])),
    }


def add_p7_saturated_fixed_boundary_equalities(
    model,
    selected,
    edges: list[tuple[int, int]],
    C: np.ndarray,
    fixed_boundary: tuple[int, ...],
) -> int:
    """Materialize the unique minimum-mean p=7 slack by coefficients.

    On the slice ``sum z_s=1``, a quadratic vanishes exactly through the
    degree-two multiples of that linear relation.  Introducing one integer
    parameter ``k_d`` per direction gives the sparse coefficient equations
    below.  They are equivalent to all 35 pointwise score equalities but
    propagate much better than dense shell rows.
    """
    if 0 not in fixed_boundary or len(fixed_boundary) != 4:
        raise ValueError("need infinity plus three finite boundary points")
    finite_points = tuple(v - 1 for v in fixed_boundary if v != 0)
    q = 3
    infinity_edges = [selected[edges.index((0, u + 1))] for u in range(49)]
    infinity_count = model.new_int_var(0, 29, "p7_saturated_infinity_count")
    model.add(infinity_count == sum(infinity_edges))
    model.add_allowed_assignments(
        [infinity_count], [[value] for value in (5, 11, 17, 23, 29)]
    )
    count = 1
    for direction in projective_directions(7):
        eps, labels = field_direction_data(7, direction)
        fibre_counts = [0] * 7
        for u in finite_points:
            fibre_counts[labels[u]] += 1
        B = {s for s, value in enumerate(fibre_counts) if value & 1}
        if len(B) not in (1, 3):
            raise AssertionError("three finite points must give one or three odd fibres")
        star_counts = [
            sum(infinity_edges[u] for u, label in enumerate(labels) if label == s)
            for s in range(7)
        ]
        parallel = sum(
            selected[j]
            for j, (a, b) in enumerate(edges)
            if a != 0 and labels[a - 1] == labels[b - 1]
        )
        k_d = model.new_int_var(-20, 30, f"p7_saturated_k_{direction[0]}_{direction[1]}")
        base = 5 if len(B) == 1 else 2
        model.add(parallel == base + q * k_d - infinity_count)
        count += 1
        for s, t in itertools.combinations(range(7), 2):
            signed_cross = sum(
                eps * int(C[a, b]) * selected[j]
                for j, (a, b) in enumerate(edges)
                if a != 0 and {labels[a - 1], labels[b - 1]} == {s, t}
            )
            if len(B) == 1:
                linear_s = int(s in B)
                linear_t = int(t in B)
                target_pair = 0
            else:
                linear_s = -int(s in B)
                linear_t = -int(t in B)
                target_pair = int(s in B and t in B)
            model.add(
                signed_cross
                == target_pair
                + k_d
                - star_counts[s]
                - star_counts[t]
                + linear_s
                + linear_t
            )
            count += 1
    return count


def solve_case(
    p: int,
    c_h: int,
    seconds: float,
    workers: int,
    infinity_value: int | None = None,
    fixed_boundary: tuple[int, ...] | None = None,
    seed: int = 15652001,
    shell_mode: str = "affine",
) -> dict:
    from ortools.sat.python import cp_model

    if p not in (5, 7):
        raise ValueError("the finite size-four solver is scoped to p=5,7")
    if c_h not in (-1, 1):
        raise ValueError("c_h must be +/-1")
    if infinity_value not in (None, 0, 1):
        raise ValueError("infinity_value must be 0, 1, or None")
    started = time.time()
    data = geometry(p, shell_mode)
    C = data["C"]
    n = data["n"]
    edges = data["edges"]
    signs = data["edge_signs"]
    edge_index = {edge: j for j, edge in enumerate(edges)}
    distinguished = edge_index[(0, 1)]
    if fixed_boundary is not None:
        fixed_boundary = tuple(sorted(fixed_boundary))
        if len(fixed_boundary) != 4 or len(set(fixed_boundary)) != 4:
            raise ValueError("fixed_boundary must have four distinct vertices")
        if not all(0 <= v < n for v in fixed_boundary):
            raise ValueError("fixed boundary vertex outside the graph")
        implied_infinity = int(0 in fixed_boundary)
        if infinity_value is not None and infinity_value != implied_infinity:
            raise ValueError("fixed_boundary contradicts infinity_value")
        infinity_value = implied_infinity

    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"edge_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == 4 * p + 1)
    model.add(selected[distinguished] == 1)

    incident_indices = [[] for _ in range(n)]
    for j, (a, b) in enumerate(edges):
        incident_indices[a].append(j)
        incident_indices[b].append(j)
    if fixed_boundary is not None:
        fixed_set = set(fixed_boundary)
        boundary = [int(v in fixed_set) for v in range(n)]
        for v in range(n):
            incident = [selected[j] for j in incident_indices[v]]
            if boundary[v]:
                model.add_bool_xor(incident)
            else:
                model.add_bool_xor([~incident[0], *incident[1:]])
    else:
        boundary = [model.new_bool_var(f"boundary_{v}") for v in range(n)]
        for v in range(n):
            model.add_modulo_equality(
                boundary[v], sum(selected[j] for j in incident_indices[v]), 2
            )
        model.add(sum(boundary) == 4)
        if infinity_value is not None:
            model.add(boundary[0] == infinity_value)

    negative_edges = [selected[j] for j, sign in enumerate(signs) if sign == -1]
    if c_h == -1:
        model.add_bool_xor(negative_edges)
    else:
        model.add_bool_xor([~negative_edges[0], *negative_edges[1:]])

    # Exact directional means and parity-floor lower bounds.  Their type
    # sums are identities, not relaxations, and substantially strengthen
    # propagation in the saturated p=7 branch.
    budget = (p + 1) ** 2 // 2
    direction_metadata = []
    half_means_by_type = {-1: [], 1: []}
    p7_saturated_profile = bool(
        p == 7 and c_h == 1 and fixed_boundary is not None and infinity_value == 1
    )
    for d, direction in enumerate(projective_directions(p)):
        if p7_saturated_profile:
            continue
        eps, labels = field_direction_data(p, direction)
        if fixed_boundary is not None:
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
        else:
            fibre_parities = []
            for s in range(p):
                parity = model.new_bool_var(f"fibre_parity_{d}_{s}")
                vertices = [
                    boundary[1 + u] for u, label in enumerate(labels) if label == s
                ]
                model.add_modulo_equality(parity, sum(vertices), 2)
                fibre_parities.append(parity)
            b_value = model.new_int_var(0, 4, f"odd_fibres_{d}")
            model.add(b_value == sum(fibre_parities))
            floor = model.new_int_var(0, 2 * p, f"parity_floor_{d}")
            table = []
            for infinity_bit in (0, 1):
                for b in range(5):
                    sign = -eps * c_h
                    if infinity_bit:
                        sign *= eps
                    if b & 1:
                        sign *= -1
                    phase = int(sign == -1)
                    table.append(
                        [infinity_bit, b, scaled_direction_floor(p, b, phase)]
                    )
            model.add_allowed_assignments([boundary[0], b_value, floor], table)

        coefficients = []
        for j, (a, b) in enumerate(edges):
            if a == 0:
                coefficient = 1
            else:
                la, lb = labels[a - 1], labels[b - 1]
                coefficient = p if la == lb else -eps * int(C[a, b])
            coefficients.append(coefficient)
        half_mean = model.new_int_var(0, p * (4 * p + 1), f"half_a_{d}")
        model.add(
            2 * half_mean
            == sum(coefficient * selected[j] for j, coefficient in enumerate(coefficients))
            - 3 * p
        )
        model.add(2 * half_mean >= floor)
        half_means_by_type[eps].append(half_mean)
        direction_metadata.append((direction, eps, b_value, floor, half_mean))
    if not p7_saturated_profile:
        for eps in (-1, 1):
            model.add(2 * sum(half_means_by_type[eps]) == budget)

    saturated_equalities = 0
    if p7_saturated_profile:
        saturated_equalities = add_p7_saturated_fixed_boundary_equalities(
            model, selected, edges, C, fixed_boundary
        )

    n_score_constraints = 0
    if saturated_equalities == 0:
        score_limit = 2 * p - 1
        for eps in (-1, 1):
            for row in data["features"][eps]:
                if eps == 1:
                    bad = np.flatnonzero(row < 0).tolist()
                else:
                    bad = np.flatnonzero(row > 0).tolist()
                model.add(sum(selected[j] for j in bad) <= score_limit)
                n_score_constraints += 1

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    out = {
        "experiment": "residual_boundary_four_lift_cpsat",
        "status": (
            "exact_affine_edge_model_not_full_shell"
            if shell_mode == "affine"
            else "exact_full_shell_edge_model"
        ),
        "shell_mode": shell_mode,
        "p": p,
        "c_H": c_h,
        "boundary_size": 4,
        "infinity_value": infinity_value,
        "fixed_boundary": list(fixed_boundary) if fixed_boundary is not None else None,
        "distinguished_edge": [0, 1],
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "n_edge_variables": len(selected),
        "n_score_constraints": n_score_constraints,
        "type_budget": budget,
        "p7_saturated_slack_equalities": saturated_equalities,
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
        out["direction_rows"] = [
            {
                "direction": list(direction),
                "eps": eps,
                "b": b_value if isinstance(b_value, int) else solver.value(b_value),
                "floor": floor if isinstance(floor, int) else solver.value(floor),
                "a": 2 * solver.value(half_mean),
            }
            for direction, eps, b_value, floor, half_mean in direction_metadata
        ]
        out["witness_audit"] = verify_witness(
            p, c_h, chosen, infinity_value, fixed_boundary, shell_mode
        )
        if not out["witness_audit"]["valid"]:
            raise AssertionError("affine witness failed independent verification")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, choices=(5, 7), required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--infinity", type=int, choices=(0, 1))
    parser.add_argument("--fixed-boundary", type=int, nargs=4)
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=15652001)
    parser.add_argument("--shell-mode", choices=("affine", "full"), default="affine")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve_case(
        args.p,
        args.c_h,
        args.seconds,
        args.workers,
        args.infinity,
        tuple(args.fixed_boundary) if args.fixed_boundary is not None else None,
        args.seed,
        args.shell_mode,
    )
    print(json.dumps(out, indent=2), flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
