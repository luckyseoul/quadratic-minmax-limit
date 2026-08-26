#!/usr/bin/env python3
"""Exact sparse model for the remaining p=7 four-finite boundaries.

The saturated model fixes each directional slack at its minimum catalog.
Here a type whose floor sum is 24 has exactly eight units of scaled mean
left.  Complete Johnson-slice catalogs supply every allowed minimum or
next-mean target, and the exact type identities select targets whose means
sum to 32.  Sparse coefficient equations then couple those targets to the
29 selected Paley edges.

``INFEASIBLE`` is a rigorous exclusion for the fixed boundary.  A feasible
row is independently checked against all 280 affine halfspace scores.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
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
from p7_unsaturated_slack_catalog import (  # noqa: E402
    exact_slack_catalog_values,
    mapped_target_catalog_rows,
)
from residual_boundary_four_lift_cpsat import geometry, verify_witness  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(
        f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp"
    )
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def _value(solver, expression) -> int:
    return int(expression) if isinstance(expression, int) else int(solver.value(expression))


def add_universal_excess_target(
    model,
    direction_index: int,
    B: set[int],
    phase: int,
    use_value_table: bool = False,
    catalog_range: tuple[int, int] | None = None,
):
    """Encode the 1,764-vector mass-ten catalog without a flat table.

    Write the elevated slack as ``A_min+2Q``.  Proposition 15.649 proves
    that every admissible ``Q`` has mass ten, values in ``{0,1,2}``, one of
    four exact histograms, and satisfies the pair-incidence reconstruction
    identity below.  Conversely that identity is an exact description of
    degree-two values on ``J(7,4)``.
    """
    points = tuple(itertools.combinations(range(7), 4))
    pairs = tuple(itertools.combinations(range(7), 2))
    values = [
        model.new_int_var(0, 2, f"excess_{direction_index}_{index}")
        for index in range(35)
    ]
    model.add(sum(values) == 10)
    vertex_masses = []
    for vertex in range(7):
        mass = model.new_int_var(0, 10, f"excess_U_{direction_index}_{vertex}")
        model.add(
            mass == sum(values[index] for index, X in enumerate(points) if vertex in X)
        )
        vertex_masses.append(mass)
    pair_masses = {}
    for s, t in pairs:
        mass = model.new_int_var(0, 10, f"excess_T_{direction_index}_{s}_{t}")
        model.add(
            mass
            == sum(
                values[index]
                for index, X in enumerate(points)
                if s in X and t in X
            )
        )
        pair_masses[s, t] = mass
    for index, X in enumerate(points):
        model.add(
            6 * values[index]
            == 2 * sum(pair_masses[pair] for pair in itertools.combinations(X, 2))
            - 3 * sum(vertex_masses[vertex] for vertex in X)
            + 36
        )

    ones = []
    twos = []
    for index, value in enumerate(values):
        is_one = model.new_bool_var(f"excess_one_{direction_index}_{index}")
        is_two = model.new_bool_var(f"excess_two_{direction_index}_{index}")
        model.add(value == 1).only_enforce_if(is_one)
        model.add(value != 1).only_enforce_if(~is_one)
        model.add(value == 2).only_enforce_if(is_two)
        model.add(value != 2).only_enforce_if(~is_two)
        ones.append(is_one)
        twos.append(is_two)
    n_ones = model.new_int_var(0, 10, f"excess_n_ones_{direction_index}")
    n_twos = model.new_int_var(0, 5, f"excess_n_twos_{direction_index}")
    model.add(n_ones == sum(ones))
    model.add(n_twos == sum(twos))
    model.add_allowed_assignments(
        [n_ones, n_twos],
        [(0, 5), (6, 2), (8, 1), (10, 0)],
    )
    catalog_total = 1764
    catalog_count = catalog_total
    if use_value_table or catalog_range is not None:
        universal_rows = [
            [int(value) // 2 for value in row]
            for row in exact_slack_catalog_values(0, 0, 8)
        ]
        if len(universal_rows) != catalog_total:
            raise AssertionError("universal excess value table is incomplete")
        if catalog_range is not None:
            start, stop = catalog_range
            if not 0 <= start < stop <= catalog_total:
                raise ValueError(
                    f"invalid catalog range [{start},{stop}) for {catalog_total} rows"
                )
            universal_rows = universal_rows[start:stop]
        catalog_count = len(universal_rows)
        model.add_allowed_assignments(values, universal_rows)

    target_constant = model.new_int_var(-16, 24, f"target_constant_{direction_index}")
    target_pairs = tuple(
        model.new_int_var(-16, 24, f"target_pair_{direction_index}_{s}_{t}")
        for s, t in pairs
    )
    for index, X in enumerate(points):
        minimum = (sum(vertex in X for vertex in B) + phase) & 1
        model.add(
            target_constant
            + sum(
                target_pairs[pair_index]
                * (1 if ((s in X) == (t in X)) else -1)
                for pair_index, (s, t) in enumerate(pairs)
            )
            == 3 + 2 * minimum + 4 * values[index]
        )
    return target_constant, target_pairs, catalog_count, catalog_total


def direction_target_options(
    b: int,
    phase: int,
    B: set[int],
    type_floor_sum: int,
    elevated: bool | None = None,
) -> tuple[tuple[int, ...], ...]:
    """Rows ``(scaled_mean, constant, 21 pair coefficients)``."""
    floor = scaled_direction_floor(7, b, phase)
    if type_floor_sum == 32:
        if elevated:
            raise ValueError("a saturated type cannot elevate a direction")
        means = [floor]
    elif elevated is None:
        means = [floor, floor + 8]
    elif elevated:
        means = [floor + 8]
    else:
        means = [floor]
    rows = []
    for scaled_mean in means:
        try:
            catalog = mapped_target_catalog_rows(
                b, phase, scaled_mean, B
            )
        except (AssertionError, ValueError) as exc:
            raise ValueError(
                "missing complete target catalog for "
                f"b={b}, phase={phase}, mean={scaled_mean}"
            ) from exc
        rows.extend((scaled_mean, *target) for target in catalog)
    unique = tuple(sorted(set(rows)))
    if len(unique) != len(rows):
        raise AssertionError("direction target table contains a duplicate")
    return unique


def solve_case(
    c_h: int,
    fixed_boundary: tuple[int, ...],
    seconds: float,
    workers: int,
    seed: int,
    elevated_directions: tuple[int, ...] | None = None,
    universal_value_table: bool = False,
    catalog_ranges: dict[int, tuple[int, int]] | None = None,
    direct_score_cuts: bool = False,
    pointwise_score_equalities: bool = False,
    pointwise_only: bool = False,
    full_score_cuts: bool = False,
) -> dict:
    from ortools.sat.python import cp_model

    if c_h not in (-1, 1):
        raise ValueError("c_h must be +/-1")
    fixed_boundary = tuple(sorted(fixed_boundary))
    if len(fixed_boundary) != 4 or len(set(fixed_boundary)) != 4:
        raise ValueError("need four distinct boundary vertices")
    if 0 in fixed_boundary:
        raise ValueError("this solver handles four finite boundary vertices")

    started = time.time()
    data = geometry(7, "affine")
    C = data["C"]
    edges = data["edges"]
    edge_index = {edge: index for index, edge in enumerate(edges)}
    signs = data["edge_signs"]

    # Boundary geometry determines every parity phase and floor before the
    # edge model is built.
    direction_data = []
    type_floors = {-1: 0, 1: 0}
    for direction in projective_directions(7):
        eps, labels = field_direction_data(7, direction)
        boundary_counts = [0] * 7
        for vertex in fixed_boundary:
            boundary_counts[labels[vertex - 1]] += 1
        B = {s for s, count in enumerate(boundary_counts) if count & 1}
        b = len(B)
        phase = int(-eps * c_h == -1)
        floor = scaled_direction_floor(7, b, phase)
        type_floors[eps] += floor
        direction_data.append(
            {
                "direction": direction,
                "eps": eps,
                "labels": labels,
                "B": B,
                "b": b,
                "phase": phase,
                "floor": floor,
            }
        )
    if any(value not in (24, 32) for value in type_floors.values()):
        raise ValueError(f"boundary is outside the surviving p=7 scope: {type_floors}")
    if type_floors == {-1: 32, 1: 32}:
        raise ValueError("boundary is doubly saturated; use Proposition 15.654")
    catalog_ranges = dict(catalog_ranges or {})
    if pointwise_only:
        pointwise_score_equalities = True
    if elevated_directions is not None:
        elevated_set = set(elevated_directions)
        if len(elevated_set) != len(elevated_directions):
            raise ValueError("elevated direction indices must be distinct")
        if not all(0 <= index < 8 for index in elevated_set):
            raise ValueError("elevated direction index outside 0..7")
        for eps in (-1, 1):
            chosen = [
                index
                for index in elevated_set
                if int(direction_data[index]["eps"]) == eps
            ]
            expected = 1 if type_floors[eps] == 24 else 0
            if len(chosen) != expected:
                raise ValueError(
                    f"type {eps} needs exactly {expected} elevated directions"
                )
    else:
        elevated_set = set()
    if set(catalog_ranges) - elevated_set:
        raise ValueError("catalog ranges may restrict only elevated directions")

    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"edge_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == 29)
    model.add(selected[edge_index[(0, 1)]] == 1)

    fixed_set = set(fixed_boundary)
    incident_indices = [[] for _ in range(50)]
    for index, (a, b) in enumerate(edges):
        incident_indices[a].append(index)
        incident_indices[b].append(index)
    for vertex in range(50):
        incident = [selected[index] for index in incident_indices[vertex]]
        if vertex in fixed_set:
            model.add_bool_xor(incident)
        else:
            model.add_bool_xor([~incident[0], *incident[1:]])

    negative = [selected[index] for index, sign in enumerate(signs) if sign == -1]
    if c_h == -1:
        model.add_bool_xor(negative)
    else:
        model.add_bool_xor([~negative[0], *negative[1:]])

    infinity_edges = [selected[edge_index[(0, u + 1)]] for u in range(49)]
    infinity_count = model.new_int_var(0, 29, "infinity_count")
    model.add(infinity_count == sum(infinity_edges))
    pair_order = tuple(itertools.combinations(range(7), 2))
    means_by_type = {-1: [], 1: []}
    coefficient_constraints = 0
    pointwise_score_constraint_count = 0

    for d, row in enumerate(direction_data):
        eps = int(row["eps"])
        labels = row["labels"]
        b = int(row["b"])
        phase = int(row["phase"])
        fixed_elevated = elevated_directions is not None and d in elevated_set
        universal_excess = bool(
            fixed_elevated and (b, phase) in {(0, 0), (2, 0), (2, 1)}
        )
        if universal_excess and d not in catalog_ranges:
            scaled_mean = int(row["floor"]) + 8
            target_constant, target_pairs, option_count, option_total_count = (
                add_universal_excess_target(
                model,
                d,
                set(row["B"]),
                phase,
                universal_value_table,
                catalog_ranges.get(d),
                )
            )
            target_encoding = (
                "mass_ten_pair_incidence_reconstruction_with_value_partition"
                if d in catalog_ranges
                else "mass_ten_pair_incidence_reconstruction"
            )
        else:
            all_options = direction_target_options(
                b,
                phase,
                set(row["B"]),
                type_floors[eps],
                None if elevated_directions is None else fixed_elevated,
            )
            option_total_count = len(all_options)
            if d in catalog_ranges:
                start, stop = catalog_ranges[d]
                if not 0 <= start < stop <= option_total_count:
                    raise ValueError(
                        f"invalid catalog range [{start},{stop}) for "
                        f"direction {d} with {option_total_count} rows"
                    )
                options = all_options[start:stop]
            else:
                options = all_options
            if len(options) == 1:
                scaled_mean = int(options[0][0])
                target_constant = int(options[0][1])
                target_pairs = tuple(int(value) for value in options[0][2:])
            else:
                columns = tuple(zip(*options))
                scaled_mean = model.new_int_var(
                    min(columns[0]), max(columns[0]), f"scaled_mean_{d}"
                )
                target_constant = model.new_int_var(
                    min(columns[1]), max(columns[1]), f"target_constant_{d}"
                )
                target_pairs = tuple(
                    model.new_int_var(
                        min(columns[2 + index]),
                        max(columns[2 + index]),
                        f"target_pair_{d}_{s}_{t}",
                    )
                    for index, (s, t) in enumerate(pair_order)
                )
                model.add_allowed_assignments(
                    [scaled_mean, target_constant, *target_pairs], options
                )
            option_count = len(options)
            target_encoding = (
                "complete_target_table_partition"
                if d in catalog_ranges
                else "complete_target_table"
            )
        means_by_type[eps].append(scaled_mean)

        if not pointwise_only:
            star_counts = [
                sum(infinity_edges[u] for u, label in enumerate(labels) if label == s)
                for s in range(7)
            ]
            parallel = sum(
                selected[index]
                for index, (a, b) in enumerate(edges)
                if a != 0 and labels[a - 1] == labels[b - 1]
            )
            k_d = model.new_int_var(-30, 40, f"coefficient_kernel_{d}")
            model.add(
                parallel == target_constant + 3 * k_d - infinity_count
            )
            coefficient_constraints += 1
            for pair_index, (s, t) in enumerate(pair_order):
                signed_cross = sum(
                    eps * int(C[a, b]) * selected[index]
                    for index, (a, b) in enumerate(edges)
                    if a != 0 and {labels[a - 1], labels[b - 1]} == {s, t}
                )
                model.add(
                    signed_cross
                    == target_pairs[pair_index]
                    + k_d
                    - star_counts[s]
                    - star_counts[t]
                )
                coefficient_constraints += 1
        if pointwise_score_equalities:
            for X in itertools.combinations(range(7), 4):
                target_value = target_constant + sum(
                    target_pairs[pair_index]
                    * (1 if ((s in X) == (t in X)) else -1)
                    for pair_index, (s, t) in enumerate(pair_order)
                )
                bad_edges = []
                X_set = set(X)
                for index, (a, endpoint) in enumerate(edges):
                    y_a = eps if a == 0 else (1 if labels[a - 1] in X_set else -1)
                    y_b = 1 if labels[endpoint - 1] in X_set else -1
                    if eps * y_a * y_b * int(C[a, endpoint]) < 0:
                        bad_edges.append(selected[index])
                model.add(2 * sum(bad_edges) + target_value == 29)
                pointwise_score_constraint_count += 1
        row["target_option_count"] = option_count
        row["target_option_total_count"] = option_total_count
        row["target_encoding"] = target_encoding
        row["scaled_mean_variable"] = scaled_mean

    for eps in (-1, 1):
        model.add(sum(means_by_type[eps]) == 32)

    if direct_score_cuts and full_score_cuts:
        raise ValueError("affine and full score cuts are mutually exclusive")
    direct_score_constraint_count = 0
    full_score_constraint_count = 0
    if direct_score_cuts or full_score_cuts:
        score_data = geometry(7, "full" if full_score_cuts else "affine")
        score_limit = 13
        for eps in (-1, 1):
            for feature in score_data["features"][eps]:
                bad = [
                    selected[index]
                    for index, value in enumerate(feature)
                    if eps * int(value) < 0
                ]
                model.add(sum(bad) <= score_limit)
                if full_score_cuts:
                    full_score_constraint_count += 1
                else:
                    direct_score_constraint_count += 1

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    out = {
        "experiment": "p7_no_infinity_unsaturated_cpsat",
        "status": "exact_complete_catalog_coefficient_edge_model",
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(fixed_boundary),
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "finite_infeasibility_certificate": status == cp_model.INFEASIBLE,
        "type_floor_sums": {str(key): value for key, value in type_floors.items()},
        "coefficient_constraints": coefficient_constraints,
        "fixed_elevated_directions": (
            sorted(elevated_set) if elevated_directions is not None else None
        ),
        "universal_value_table": universal_value_table,
        "catalog_ranges": {
            str(direction): list(bounds)
            for direction, bounds in sorted(catalog_ranges.items())
        },
        "catalog_partition_basis": (
            "mapped_target_catalog_rows_lexicographic_v1"
            if catalog_ranges
            else None
        ),
        "direct_score_cuts": direct_score_cuts,
        "direct_score_constraint_count": direct_score_constraint_count,
        "full_score_cuts": full_score_cuts,
        "full_score_constraint_count": full_score_constraint_count,
        "pointwise_score_equalities": pointwise_score_equalities,
        "pointwise_score_constraint_count": pointwise_score_constraint_count,
        "pointwise_only": pointwise_only,
        "direction_rows": [
            {
                "direction": list(row["direction"]),
                "eps": row["eps"],
                "b": row["b"],
                "phase": row["phase"],
                "floor": row["floor"],
                "target_option_count": row["target_option_count"],
                "target_option_total_count": row["target_option_total_count"],
                "target_encoding": row["target_encoding"],
                "scaled_mean": _value(solver, row["scaled_mean_variable"])
                if feasible
                else None,
            }
            for row in direction_data
        ],
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
        out["witness_audit"] = verify_witness(
            7,
            c_h,
            chosen,
            0,
            fixed_boundary,
            "full" if full_score_cuts else "affine",
        )
        if not out["witness_audit"]["valid"]:
            raise AssertionError("catalog coefficient witness failed direct audit")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs=4, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15655001)
    parser.add_argument("--elevated-directions", type=int, nargs="*")
    parser.add_argument("--universal-value-table", action="store_true")
    parser.add_argument(
        "--catalog-range",
        type=int,
        nargs=3,
        action="append",
        metavar=("DIRECTION", "START", "STOP"),
        help="restrict one elevated direction to a half-open catalog interval",
    )
    parser.add_argument(
        "--direct-score-cuts",
        action="store_true",
        help="add all 280 redundant affine score cardinality bounds",
    )
    parser.add_argument(
        "--full-score-cuts",
        action="store_true",
        help="add every cached full-eigenshell score cardinality bound",
    )
    parser.add_argument(
        "--pointwise-score-equalities",
        action="store_true",
        help="add all 280 exact bad-edge/slack cardinality equalities",
    )
    parser.add_argument(
        "--pointwise-only",
        action="store_true",
        help="use exact pointwise score equalities without coefficient equations",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    catalog_ranges = {}
    for direction, start, stop in args.catalog_range or []:
        if direction in catalog_ranges:
            raise ValueError(f"duplicate catalog range for direction {direction}")
        catalog_ranges[direction] = (start, stop)
    out = solve_case(
        args.c_h,
        tuple(args.fixed_boundary),
        args.seconds,
        args.workers,
        args.seed,
        tuple(args.elevated_directions)
        if args.elevated_directions is not None
        else None,
        args.universal_value_table,
        catalog_ranges,
        args.direct_score_cuts,
        args.pointwise_score_equalities,
        args.pointwise_only,
        args.full_score_cuts,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
