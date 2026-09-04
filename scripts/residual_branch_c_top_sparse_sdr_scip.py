#!/usr/bin/env python3
"""SCIP backend for the restricted branch-C pair-collision Möbius gate.

The mathematical model is documented in
``residual_branch_c_top_sparse_sdr_cpsat.py``.  SCIP is used here because
the compact-triangle condition is a large collection of linear l1 bounds;
its LP relaxation is materially stronger than the corresponding CP-SAT
absolute-value encoding.  Every returned solution is rebuilt as a physical
graph and replayed independently.

As in the CP-SAT implementation, the sole cancellation is required to be one
isolated opposite-sign collision between two selected paired-SDR options.
This backend does not encode clean three-half ``2:1`` overlaps or any other
higher-multiplicity cancellation geometry.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from e1_gmin_m4_inversion_antisymmetric_radon import (  # noqa: E402
    _negative_edge,
    edge_radon_image,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (  # noqa: E402
    paley_direction_sign,
    paley_edge_sign,
)
from io_atomic import write_json_atomic  # noqa: E402
from scripts.residual_branch_c_aux_sdr_cpsat import _profile, build_chart  # noqa: E402
from scripts.residual_branch_c_top_sparse_sdr_cpsat import (  # noqa: E402
    PhysicalOption,
    _fixed_edge,
    _fixed_word_capacity_replay,
    _graph_from_orbit_values,
    _literal_centres,
    _option_graph,
    _projected_key,
    generate_physical_options,
)


def _option_projection(
    p: int,
    directions: tuple[tuple[int, int], ...],
    option: PhysicalOption,
) -> Counter[tuple[int, tuple[object, ...]]]:
    out: Counter[tuple[int, tuple[object, ...]]] = Counter()
    for edge in _option_graph(p, option):
        tau = paley_edge_sign(p, edge)
        for direction_index, direction in enumerate(directions):
            coefficient = paley_direction_sign(p, direction) * tau
            out[(direction_index, _projected_key(p, direction, edge))] += coefficient
    return out


def _collision_correction(
    p: int,
    directions: tuple[tuple[int, int], ...],
    orbit,
) -> Counter[tuple[int, tuple[object, ...]]]:
    out: Counter[tuple[int, tuple[object, ...]]] = Counter()
    for edge in (orbit, _negative_edge(p, orbit)):
        tau = paley_edge_sign(p, edge)
        for direction_index, direction in enumerate(directions):
            coefficient = -paley_direction_sign(p, direction) * tau
            out[(direction_index, _projected_key(p, direction, edge))] += coefficient
    return out


def solve_scale(
    p: int,
    fixed_direction_index: int | None,
    alphas: tuple[int, ...],
    c: int,
    seconds: float,
    workers: int,
    relax_sparse_budget: bool,
) -> dict[str, object]:
    from pyscipopt import Model, quicksum

    chart = build_chart(p, fixed_direction_index)
    m = (p + 1) // 2
    r = (p - 3) // 4
    options = generate_physical_options(chart, alphas, c)
    model = Model(f"branch_c_top_sparse_p{p}_c{c}")
    model.hideOutput()
    model.setParam("limits/time", seconds)
    model.setParam("parallel/maxnthreads", workers)
    model.setParam("randomization/randomseedshift", c)
    selected = tuple(
        model.addVar(vtype="B", name=f"option_{index}")
        for index in range(len(options))
    )

    target_rows: list[list[int]] = [[] for _ in range(m)]
    auxiliary_rows: list[list[int]] = [[] for _ in range(p)]
    for index, option in enumerate(options):
        target_rows[option.first_target].append(index)
        target_rows[option.second_target].append(index)
        auxiliary_rows[option.first_auxiliary].append(index)
        auxiliary_rows[option.second_auxiliary].append(index)
    if any(not row for row in target_rows):
        return {
            "c": c,
            "status": "infeasible_by_empty_target_row",
            "option_count": len(options),
        }
    for row in target_rows:
        model.addCons(quicksum(selected[index] for index in row) == 1)
    for row in auxiliary_rows:
        if row:
            model.addCons(quicksum(selected[index] for index in row) <= 1)
    model.addCons(
        quicksum(
            (
                int(chart.type_by_coordinate[option.first_auxiliary] == 1)
                + int(chart.type_by_coordinate[option.second_auxiliary] == 1)
            )
            * selected[index]
            for index, option in enumerate(options)
        )
        == m - 2
    )

    option_orbits = tuple(dict(option.orbit_coefficients) for option in options)
    # Restricted collision ansatz: exactly one isolated opposite collision
    # between two selected options.  A clean three-half 2:1 overlap is a valid
    # cancellation-one pattern but is intentionally outside this model.
    collision_pairs = []
    physical_conflicts = 0
    for first_index in range(len(options)):
        first = options[first_index]
        first_targets = {first.first_target, first.second_target}
        first_aux = {first.first_auxiliary, first.second_auxiliary}
        for second_index in range(first_index + 1, len(options)):
            second = options[second_index]
            if first_targets & {second.first_target, second.second_target}:
                continue
            if first_aux & {second.first_auxiliary, second.second_auxiliary}:
                continue
            common = set(option_orbits[first_index]) & set(option_orbits[second_index])
            same = tuple(
                orbit for orbit in common
                if option_orbits[first_index][orbit] == option_orbits[second_index][orbit]
            )
            opposite = tuple(sorted(common - set(same)))
            if same or len(opposite) > 1:
                model.addCons(selected[first_index] + selected[second_index] <= 1)
                physical_conflicts += 1
            elif len(opposite) == 1:
                both = model.addVar(
                    vtype="B", name=f"collision_{first_index}_{second_index}"
                )
                model.addCons(both <= selected[first_index])
                model.addCons(both <= selected[second_index])
                model.addCons(
                    both >= selected[first_index] + selected[second_index] - 1
                )
                collision_pairs.append(
                    (first_index, second_index, opposite[0], both)
                )
    model.addCons(quicksum(row[3] for row in collision_pairs) == 1)

    directions = projective_functionals(p)
    literal_centres = _literal_centres(chart, alphas)
    option_projections = tuple(
        _option_projection(p, directions, option) for option in options
    )
    collision_projections = tuple(
        _collision_correction(p, directions, row[2])
        for row in collision_pairs
    )
    cell_terms: defaultdict[
        tuple[int, tuple[object, ...]], list[tuple[int, object]]
    ] = defaultdict(list)
    for option_index, projection in enumerate(option_projections):
        for cell, coefficient in projection.items():
            cell_terms[cell].append((coefficient, selected[option_index]))
    for pair_index, projection in enumerate(collision_projections):
        variable = collision_pairs[pair_index][3]
        for cell, coefficient in projection.items():
            cell_terms[cell].append((coefficient, variable))

    fixed = _fixed_edge(chart, c)
    fixed_tau = paley_edge_sign(p, fixed)
    fixed_constants: defaultdict[tuple[int, tuple[object, ...]], int] = defaultdict(int)
    for direction_index, direction in enumerate(directions):
        coefficient = paley_direction_sign(p, direction) * fixed_tau
        fixed_constants[
            (direction_index, _projected_key(p, direction, fixed))
        ] += coefficient

    def expression(cell):
        return (
            quicksum(coefficient * variable for coefficient, variable in cell_terms[cell])
            + fixed_constants[cell]
        )

    row_l1_variables = []
    parallel_expressions = []
    residual_expressions = []
    for direction_index, direction in enumerate(directions):
        direction_sign = paley_direction_sign(p, direction)
        parallel = expression((direction_index, ("P",)))
        parallel_expressions.append(parallel)
        if direction_sign == 1:
            model.addCons(parallel >= 3)
            literal = literal_centres[direction_index]
        else:
            model.addCons(parallel >= r + 2)
            literal = None
        cells = []
        for left in range(p):
            for right in range(left + 1, p):
                cell = (direction_index, ("K", left, right))
                value = expression(cell)
                if literal is not None and literal in (left, right):
                    value += 1
                cells.append(value)
        residual_expressions.append(tuple(cells))
        if relax_sparse_budget:
            absolutes = []
            for cell_index, value in enumerate(cells):
                absolute = model.addVar(
                    lb=0.0,
                    ub=2.0 * m,
                    vtype="C",
                    name=f"abs_{direction_index}_{cell_index}",
                )
                model.addCons(absolute >= value)
                model.addCons(absolute >= -value)
                absolutes.append(absolute)
            row_l1 = model.addVar(
                lb=0.0,
                ub=3.0 * m * p,
                vtype="C",
                name=f"row_l1_{direction_index}",
            )
            model.addCons(row_l1 == quicksum(absolutes))
            row_l1_variables.append(row_l1)
    if relax_sparse_budget:
        model.setObjective(quicksum(row_l1_variables), "minimize")
    else:
        model.setObjective(0.0, "minimize")

    started = time.monotonic()
    cut_iterations = 0
    l1_cuts_added = 0
    maximum_seen_violation = 0
    accepted_solution = None
    while True:
        remaining = seconds - (time.monotonic() - started)
        if remaining <= 0:
            break
        model.setParam("limits/time", remaining)
        model.optimize()
        solution = model.getBestSol()
        if solution is None or relax_sparse_budget:
            accepted_solution = solution
            break

        violated = []
        for direction_index, cells in enumerate(residual_expressions):
            values = tuple(
                int(round(model.getSolVal(solution, value))) for value in cells
            )
            l1 = sum(abs(value) for value in values)
            parallel = int(
                round(model.getSolVal(solution, parallel_expressions[direction_index]))
            )
            budget = 3 * (parallel - 3)
            if l1 > budget:
                violated.append((direction_index, values, l1 - budget))
                maximum_seen_violation = max(maximum_seen_violation, l1 - budget)
        if not violated:
            accepted_solution = solution
            break

        # l1(v)=max_{sigma in {+1,-1}^n} sigma.v.  The sign vector of the
        # current violating row supplies a globally valid separating cut.
        model.freeTransform()
        for direction_index, values, _violation in violated:
            signs = tuple(1 if value >= 0 else -1 for value in values)
            model.addCons(
                quicksum(
                    sign * value
                    for sign, value in zip(
                        signs,
                        residual_expressions[direction_index],
                        strict=True,
                    )
                )
                <= 3 * (parallel_expressions[direction_index] - 3),
                name=f"l1_cut_{cut_iterations}_{direction_index}",
            )
            l1_cuts_added += 1
        cut_iterations += 1

    wall = time.monotonic() - started
    status = str(model.getStatus())
    result: dict[str, object] = {
        "c": c,
        "status": status,
        "option_count": len(options),
        "candidate_orbit_count": len(set().union(*(set(row) for row in option_orbits))),
        "physical_conflict_pair_count": physical_conflicts,
        "single_opposite_collision_pair_count": len(collision_pairs),
        "collision_model_scope": (
            "one isolated opposite collision between two selected paired-SDR "
            "options; no orbit shared by three or more constituent halves"
        ),
        "higher_multiplicity_overlap_models_included": False,
        "compact_triangle_l1_budget_enforced": not relax_sparse_budget,
        "wall_seconds": wall,
        "nodes": model.getNNodes(),
        "primal_bound": model.getPrimalbound(),
        "dual_bound": model.getDualbound(),
        "gap": model.getGap(),
        "l1_cut_iterations": cut_iterations,
        "l1_cuts_added": l1_cuts_added,
        "maximum_seen_l1_violation": maximum_seen_violation,
    }
    solution = accepted_solution
    if solution is None:
        return result

    chosen_indices = tuple(
        index
        for index, variable in enumerate(selected)
        if model.getSolVal(solution, variable) > 0.5
    )
    chosen = tuple(options[index] for index in chosen_indices)
    replay_counter = Counter()
    for option in chosen:
        replay_counter.update(dict(option.orbit_coefficients))
    if any(abs(value) > 1 for value in replay_counter.values()):
        raise ArithmeticError("the SCIP witness failed ternary replay")
    replayed_orbits = {
        orbit: value for orbit, value in replay_counter.items() if value
    }
    cancellations = (
        m * (p - 1) - sum(abs(value) for value in replay_counter.values())
    ) // 2
    if cancellations != 1:
        raise ArithmeticError("the SCIP witness failed the cancellation replay")
    graph = _graph_from_orbit_values(p, replayed_orbits, fixed)
    image = edge_radon_image(
        p, {edge: paley_edge_sign(p, edge) for edge in graph}
    )
    replay_rows = []
    for direction_index, direction in enumerate(directions):
        direction_sign = paley_direction_sign(p, direction)
        parallel = direction_sign * image.get(("P", direction_index), 0)
        literal = literal_centres.get(direction_index)
        residual = []
        for left in range(p):
            for right in range(left + 1, p):
                value = direction_sign * image.get(
                    ("K", direction_index, left, right), 0
                )
                if literal is not None and literal in (left, right):
                    value += 1
                residual.append(value)
        l1 = sum(abs(value) for value in residual)
        budget = 3 * (parallel - 3)
        replay_rows.append(
            {
                "direction_index": direction_index,
                "direction_sign": direction_sign,
                "parallel_count": parallel,
                "literal_centre": literal,
                "residual_coefficient_l1": l1,
                "compact_triangle_l1_budget": budget,
                "within_budget": l1 <= budget,
            }
        )
    fixed_word = _fixed_word_capacity_replay(
        chart, alphas, c, graph, replayed_orbits
    )
    graph_bytes = json.dumps(graph, separators=(",", ":")).encode()
    if len(graph) != m * (p - 1) - 1:
        raise ArithmeticError("the SCIP graph has the wrong top-end size")
    if not relax_sparse_budget and not all(row["within_budget"] for row in replay_rows):
        raise ArithmeticError("the SCIP graph failed its l1 constraints")
    result["witness"] = {
        "selected_option_indices": list(chosen_indices),
        "selected_options": [
            {
                "targets": [option.first_target, option.second_target],
                "signs": [option.first_sign, option.second_sign],
                "auxiliaries": [option.first_auxiliary, option.second_auxiliary],
            }
            for option in chosen
        ],
        "used_orbit_count": len(replayed_orbits),
        "graph_edge_count": len(graph),
        "graph_sha256": hashlib.sha256(graph_bytes).hexdigest(),
        "graph_edges": [[list(edge[0]), list(edge[1])] for edge in graph],
        "row_replay": replay_rows,
        "total_residual_coefficient_l1": sum(
            row["residual_coefficient_l1"] for row in replay_rows
        ),
        "total_compact_triangle_l1_budget": sum(
            row["compact_triangle_l1_budget"] for row in replay_rows
        ),
        "all_rows_within_compact_triangle_l1_budget": all(
            row["within_budget"] for row in replay_rows
        ),
        "fixed_word_atom_capacity_replay": fixed_word,
        "exact_physical_replay": True,
    }
    return result


def run(args: argparse.Namespace) -> dict[str, object]:
    chart = build_chart(args.p, args.fixed_direction_index)
    alphas = _profile(
        args.alpha_mode,
        args.p,
        len(chart.hard_coordinates),
        args.profile_seed,
    )
    skipped = {
        int(value) for value in args.skip_scales.split(",") if value.strip()
    }
    scales = tuple(
        c
        for c in range(1, args.p)
        if (c - 1) % args.shard_count == args.shard_index
        and c not in skipped
    )
    records = []
    started = time.monotonic()
    for c in scales:
        record = solve_scale(
            args.p,
            args.fixed_direction_index,
            alphas,
            c,
            args.seconds_per_scale,
            args.workers,
            args.relax_sparse_budget,
        )
        records.append(record)
        print(
            f"c={c} status={record['status']} witness={'witness' in record} "
            f"cuts={record.get('l1_cuts_added', 0)} wall={record.get('wall_seconds', 0):.3f}",
            flush=True,
        )
        if args.stop_on_solution and "witness" in record:
            break
    alpha_bytes = b"".join(
        int(value).to_bytes(4, byteorder="little", signed=False)
        for value in alphas
    )
    return {
        "schema": "residual_branch_c_top_sparse_sdr_scip_v1",
        "scope": (
            "one centre profile and fixed direction in the paired-SDR isolated "
            "two-option collision subfamily; necessary compact-triangle l1 gate"
        ),
        "collision_model_scope": (
            "one isolated opposite collision between two selected paired-SDR "
            "options; clean three-half 2:1 overlaps are outside scope"
        ),
        "higher_multiplicity_overlap_models_included": False,
        "p": args.p,
        "fixed_direction_index": chart.fixed_direction_index,
        "alpha_mode": args.alpha_mode,
        "profile_seed": args.profile_seed,
        "alphas": list(alphas),
        "alpha_sha256": hashlib.sha256(alpha_bytes).hexdigest(),
        "shard": [args.shard_index, args.shard_count],
        "scale_count": len(records),
        "solution_count": sum("witness" in row for row in records),
        "infeasible_count": sum(row["status"] == "infeasible" for row in records),
        "records": records,
        "elapsed_seconds": time.monotonic() - started,
        "interpretation": (
            "A solution passes only the necessary compact-triangle l1 gate; exact atom "
            "decomposition remains. Infeasibility excludes only the recorded "
            "profile/scale/fixed direction inside the isolated pair-collision "
            "paired-SDR subfamily; higher-multiplicity overlaps remain open."
        ),
    }


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, default=31)
    parser.add_argument("--fixed-direction-index", type=int)
    parser.add_argument(
        "--alpha-mode", choices=("constant", "ramp", "alternating", "random"), default="constant"
    )
    parser.add_argument("--profile-seed", type=int, default=15766)
    parser.add_argument("--seconds-per-scale", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--skip-scales", default="")
    parser.add_argument("--stop-on-solution", action="store_true")
    parser.add_argument("--relax-sparse-budget", action="store_true")
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    result = run(args)
    if args.output is None:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        write_json_atomic(args.output, result)
        print(
            f"wrote {args.output}: scales={result['scale_count']} "
            f"solutions={result['solution_count']} infeasible={result['infeasible_count']}"
        )


if __name__ == "__main__":
    main()
