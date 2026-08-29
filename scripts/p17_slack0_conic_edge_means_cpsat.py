#!/usr/bin/env python3
"""Exact edge/mean test for the two p=17 slack-zero conic profiles.

Every affine 16-arc in PG(2,17) is conic-minus-two.  Exact comparison with
the phase-labelled arithmetic ledger leaves two tangent-at-infinity
profiles.  This model fixes canonical representatives and imposes the full
69-edge count, odd-degree boundary, Paley product sign, and every exact
same-type directional mean allocation.  INFEASIBLE is an unconditional
exclusion at this necessary-condition level; FEASIBLE is only a witness to
these constraints and UNKNOWN has no mathematical force.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)


P = 17
Q = P * P
N = Q + 1
H_SIZE = 4 * P + 1
CONIC = tuple((t * t % P, t, 1) for t in range(P)) + ((1, 0, 0),)


def atomic_write(path: Path, payload: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(a * b for a, b in zip(left, right)) % P


def determinant(rows: tuple[tuple[int, int, int], ...]) -> int:
    a, b, c = rows
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    ) % P


def chart(line: tuple[int, int, int]) -> tuple[tuple[int, int, int], ...]:
    basis = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    for first in basis:
        for second in basis:
            rows = (first, second, line)
            if first != second and determinant(rows):
                return rows
    raise ArithmeticError("failed to construct an affine chart")


def affine_point(
    matrix: tuple[tuple[int, int, int], ...],
    point: tuple[int, int, int],
) -> tuple[int, int]:
    image = tuple(dot(row, point) for row in matrix)
    inverse = pow(image[2], -1, P)
    return image[0] * inverse % P, image[1] * inverse % P


def canonical_boundary(case: int) -> tuple[int, ...]:
    if case not in (0, 1):
        raise ValueError("case must be zero or one")
    line = (1, 0, 0)
    removed = {0, 1 if case == 0 else 2}
    if {index for index, point in enumerate(CONIC) if dot(line, point) == 0} != {0}:
        raise ArithmeticError("canonical tangent line changed")
    matrix = chart(line)
    points = [
        affine_point(matrix, point)
        for index, point in enumerate(CONIC)
        if index not in removed
    ]
    boundary = tuple(sorted(y * P + x for x, y in points))
    if len(boundary) != 16 or len(set(boundary)) != 16:
        raise ArithmeticError("canonical conic boundary changed")
    return boundary


def conference_data() -> tuple[list[tuple[int, int]], list[int]]:
    _q, _mul, _add, chi, _frob, _norm, _ia, _ib = field_ctx(P)
    edges: list[tuple[int, int]] = []
    signs: list[int] = []
    for left in range(N):
        for right in range(left + 1, N):
            edges.append((left, right))
            if left == 0:
                signs.append(1)
            else:
                a, b = left - 1, right - 1
                difference = (
                    ((a % P - b % P) % P)
                    + ((a // P - b // P) % P) * P
                )
                signs.append(int(chi(difference)))
    return edges, signs


def boundary_directions(boundary: tuple[int, ...], c_h: int) -> list[dict[str, object]]:
    rows = []
    for direction in projective_directions(P):
        eps, labels = field_direction_data(P, direction)
        counts = [0] * P
        for vertex in boundary:
            counts[labels[vertex]] += 1
        b = sum(value & 1 for value in counts)
        phase = int(-eps * c_h == -1)
        rows.append(
            {
                "direction": tuple(direction),
                "eps": int(eps),
                "labels": labels,
                "b": b,
                "phase": phase,
            }
        )
    return rows


def solve(
    case: int,
    seconds: float,
    workers: int,
    seed: int,
    branch: str | None = None,
    boundary_override: tuple[int, ...] | None = None,
    c_h_override: int | None = None,
    linearization_level: int = 0,
    fixed_infinity_degree: int | None = None,
    fixed_infinity_neighbors: tuple[int, ...] | None = None,
) -> dict[str, object]:
    from ortools.sat.python import cp_model

    started = time.time()
    # The defaults are deterministic examples of the two labelled profiles.
    # Overrides let an external complete affine-orbit census test every
    # inequivalent embedding instead of assuming example transitivity.
    c_h = (1 if case == 0 else -1) if c_h_override is None else c_h_override
    boundary = canonical_boundary(case) if boundary_override is None else boundary_override
    if c_h not in (-1, 1):
        raise ValueError("c_h must be -1 or 1")
    if len(boundary) != 16 or len(set(boundary)) != 16 or not all(0 <= x < Q for x in boundary):
        raise ValueError("boundary must contain sixteen distinct affine points")
    boundary_set = set(boundary)
    direction_rows = boundary_directions(boundary, c_h)
    observed = {
        phase: dict(
            sorted(Counter(int(row["b"]) for row in direction_rows if row["phase"] == phase).items())
        )
        for phase in (0, 1)
    }
    expected = (
        {0: {0: 1, 2: 7, 16: 1}, 1: {2: 9}}
        if case == 0
        else {0: {0: 1, 2: 8}, 1: {2: 8, 16: 1}}
    )
    if observed != expected:
        raise ArithmeticError("canonical conic profile changed")

    edges, signs = conference_data()
    edge_index = {edge: index for index, edge in enumerate(edges)}
    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"e_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == H_SIZE)

    incident: list[list[object]] = [[] for _ in range(N)]
    for variable, (left, right) in zip(selected, edges):
        incident[left].append(variable)
        incident[right].append(variable)
    model.add_bool_xor([*incident[0], model.new_constant(1)])
    for vertex in range(Q):
        literals = incident[vertex + 1]
        if vertex in boundary_set:
            model.add_bool_xor(literals)
        else:
            model.add_bool_xor([*literals, model.new_constant(1)])

    negative = [variable for variable, sign in zip(selected, signs) if sign == -1]
    if c_h == -1:
        model.add_bool_xor(negative)
    else:
        model.add_bool_xor([*negative, model.new_constant(1)])

    infinity_degree = sum(
        selected[edge_index[(0, vertex + 1)]] for vertex in range(Q)
    )
    if fixed_infinity_degree is not None:
        model.add(infinity_degree == fixed_infinity_degree)
    if fixed_infinity_neighbors is not None:
        neighbor_set = set(fixed_infinity_neighbors)
        if not all(0 <= vertex < Q for vertex in neighbor_set):
            raise ValueError("infinity neighbors must lie in the affine plane")
        for vertex in range(Q):
            model.add(
                selected[edge_index[(0, vertex + 1)]]
                == int(vertex in neighbor_set)
            )

    allocation_branch = model.new_bool_var("case_one_second_allocation") if case == 1 else None
    phase_one_elevated = []
    phase_zero_elevated = []
    phase_one_elevated_by_direction = {}
    phase_zero_elevated_by_direction = {}
    means_by_phase: dict[int, list[object]] = {0: [], 1: []}
    parallel_counts = []
    rendered_rows = []
    rigid_coefficient_identities = 0
    for index, row in enumerate(direction_rows):
        eps = int(row["eps"])
        labels = row["labels"]
        b = int(row["b"])
        phase = int(row["phase"])
        fibres = [[] for _ in range(P)]
        for vertex, label in enumerate(labels):
            fibres[label].append(vertex)
        odd_fibres = {
            fibre
            for fibre, points in enumerate(fibres)
            if sum(vertex in boundary_set for vertex in points) & 1
        }
        star = [
            sum(selected[edge_index[(0, vertex + 1)]] for vertex in points)
            for points in fibres
        ]
        cross_terms: dict[tuple[int, int], list[object]] = {
            (left, right): []
            for left in range(P)
            for right in range(left + 1, P)
        }
        cross_selected: dict[tuple[int, int], list[object]] = {
            pair: [] for pair in cross_terms
        }
        terms = []
        parallel_selected = []
        for edge_variable, (left, right), sign in zip(selected, edges, signs):
            if left == 0:
                coefficient = 1
            elif labels[left - 1] == labels[right - 1]:
                coefficient = P
                parallel_selected.append(edge_variable)
            else:
                coefficient = -eps * sign
                first, second = labels[left - 1], labels[right - 1]
                pair = (first, second) if first < second else (second, first)
                cross_terms[pair].append(eps * sign * edge_variable)
                cross_selected[pair].append(edge_variable)
            terms.append(coefficient * edge_variable)
        parallel = model.new_int_var(0, H_SIZE, f"parallel_{index}")
        model.add(parallel == sum(parallel_selected))
        parallel_counts.append(parallel)
        mean = model.new_int_var(0, 4 * P, f"mean_{index}")
        model.add(mean == sum(terms) - 3 * P)
        means_by_phase[phase].append(mean)

        role = "fixed"
        rigid_target_sign: int | None = None
        rigid_guard = None
        if case == 0 and phase == 0:
            model.add(mean == {0: 0, 2: 18, 16: 36}[b])
            if b == 0:
                rigid_target_sign = 0
            elif b == 2:
                rigid_target_sign = -1
        elif case == 0 and phase == 1:
            elevated = model.new_bool_var(f"phase_one_elevated_{index}")
            model.add(mean == 16 + 18 * elevated)
            phase_one_elevated.append(elevated)
            phase_one_elevated_by_direction[index] = elevated
            role = "one_of_nine_elevated"
            rigid_target_sign = 1
            rigid_guard = ~elevated
        elif case == 1 and phase == 1:
            model.add(mean == {2: 16, 16: 34}[b])
            if b == 2:
                rigid_target_sign = 1
        elif case == 1 and b == 0:
            assert allocation_branch is not None
            model.add(mean == 18 * allocation_branch)
            role = "allocation_branch_anchor"
            rigid_target_sign = 0
            rigid_guard = ~allocation_branch
        else:
            assert allocation_branch is not None and phase == 0 and b == 2
            elevated = model.new_bool_var(f"phase_zero_elevated_{index}")
            model.add(elevated + allocation_branch <= 1)
            model.add(mean == 18 + 18 * elevated)
            phase_zero_elevated.append(elevated)
            phase_zero_elevated_by_direction[index] = elevated
            role = "one_elevated_unless_anchor_branch"
            rigid_target_sign = -1
            rigid_guard = ~elevated

        if rigid_target_sign is not None:
            gauge = model.new_int_var(-H_SIZE, H_SIZE, f"gauge_{index}")
            target_sum = rigid_target_sign * (b * (b - 1) // 2)
            # Sum the 136 cell identities explicitly.  This exposes the
            # parallel-count/gauge arithmetic that otherwise remains hidden
            # behind thousands of edge expressions during presolve.
            aggregate_identity = model.add(
                infinity_degree + P * parallel - 3 * P - mean
                == target_sum
                + (P * (P - 1) // 2) * gauge
                - (P - 1) * infinity_degree
            )
            if rigid_guard is not None:
                aggregate_identity.only_enforce_if(rigid_guard)
            absolute_rhs = []
            for first in range(P):
                for second in range(first + 1, P):
                    target = (
                        rigid_target_sign
                        if first in odd_fibres and second in odd_fibres
                        else 0
                    )
                    rhs = target + gauge - star[first] - star[second]
                    constraint = model.add(
                        sum(cross_terms[first, second])
                        == rhs
                    )
                    positive_cut = model.add(
                        sum(cross_selected[first, second]) >= rhs
                    )
                    negative_cut = model.add(
                        sum(cross_selected[first, second]) >= -rhs
                    )
                    absolute = model.new_int_var(
                        0, H_SIZE, f"absolute_rhs_{index}_{first}_{second}"
                    )
                    absolute_positive = model.add(absolute >= rhs)
                    absolute_negative = model.add(absolute >= -rhs)
                    absolute_rhs.append(absolute)
                    rigid_coefficient_identities += 1
                    if rigid_guard is not None:
                        constraint.only_enforce_if(rigid_guard)
                        positive_cut.only_enforce_if(rigid_guard)
                        negative_cut.only_enforce_if(rigid_guard)
                        absolute_positive.only_enforce_if(rigid_guard)
                        absolute_negative.only_enforce_if(rigid_guard)
            aggregate_l1 = model.add(
                sum(absolute_rhs)
                <= H_SIZE - infinity_degree - sum(parallel_selected)
            )
            if rigid_guard is not None:
                aggregate_l1.only_enforce_if(rigid_guard)
        rendered_rows.append(
            {
                "direction": list(row["direction"]),
                "eps": eps,
                "phase": phase,
                "b": b,
                "role": role,
            }
        )

    if case == 0:
        model.add(sum(phase_one_elevated) == 1)
        if branch is not None:
            direction = int(branch)
            if direction not in phase_one_elevated_by_direction:
                raise ValueError("case-zero branch is not a phase-one b=2 direction")
            model.add(phase_one_elevated_by_direction[direction] == 1)
    else:
        assert allocation_branch is not None
        model.add(sum(phase_zero_elevated) + allocation_branch == 1)
        if branch == "anchor":
            model.add(allocation_branch == 1)
        elif branch is not None:
            direction = int(branch)
            if direction not in phase_zero_elevated_by_direction:
                raise ValueError("case-one branch is not a phase-zero b=2 direction")
            model.add(allocation_branch == 0)
            model.add(phase_zero_elevated_by_direction[direction] == 1)
    for phase in (0, 1):
        model.add(sum(means_by_phase[phase]) == 162)
    # Every finite affine edge is parallel to exactly one projective
    # direction; infinity-star edges are parallel to none.
    model.add(sum(parallel_counts) == H_SIZE - infinity_degree)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    solver.parameters.linearization_level = int(linearization_level)
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    result: dict[str, object] = {
        "experiment": "p17_slack0_conic_edge_means_cpsat",
        "status": (
            "exact_fixed_boundary_edge_parity_product_mean_and_"
            "rigid_coefficient_model"
        ),
        "p": P,
        "c_H": c_h,
        "case": case,
        "fixed_allocation_branch": branch,
        "fixed_infinity_degree": fixed_infinity_degree,
        "fixed_infinity_neighbors": (
            list(fixed_infinity_neighbors)
            if fixed_infinity_neighbors is not None
            else None
        ),
        "fixed_boundary": list(boundary),
        "phase_profiles_b": observed,
        "direction_rows": rendered_rows,
        "edge_variables": len(selected),
        "rigid_coefficient_identities": rigid_coefficient_identities,
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "workers": workers,
        "linearization_level": linearization_level,
        "seed": seed,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
    }
    if feasible:
        result["allocation_branch"] = (
            solver.value(allocation_branch) if allocation_branch is not None else None
        )
        result["chosen_edges_H"] = [
            list(edge)
            for edge, variable in zip(edges, selected)
            if solver.value(variable)
        ]
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=int, choices=(0, 1), required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15700001)
    parser.add_argument(
        "--branch",
        help="fix the nonrigid allocation: a direction index, or 'anchor' in case one",
    )
    parser.add_argument(
        "--boundary",
        help="comma-separated sixteen-point affine boundary override",
    )
    parser.add_argument("--c-h", type=int, choices=(-1, 1))
    parser.add_argument("--linearization-level", type=int, choices=(0, 1, 2), default=0)
    parser.add_argument("--infinity-degree", type=int)
    parser.add_argument("--infinity-neighbors")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    boundary = (
        tuple(int(value) for value in args.boundary.split(","))
        if args.boundary is not None
        else None
    )
    result = solve(
        args.case,
        args.seconds,
        args.workers,
        args.seed,
        args.branch,
        boundary,
        args.c_h,
        args.linearization_level,
        args.infinity_degree,
        (
            tuple(int(value) for value in args.infinity_neighbors.split(","))
            if args.infinity_neighbors is not None
            else None
        ),
    )
    if args.output is not None:
        atomic_write(args.output, result)
    rendered = dict(result)
    rendered.pop("chosen_edges_H", None)
    print(json.dumps(rendered, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
