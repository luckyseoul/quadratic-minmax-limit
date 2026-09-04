#!/usr/bin/env python3
"""Exact restricted pair-collision gate for the branch-C top endpoint.

This is the next layer after :mod:`residual_branch_c_aux_sdr_cpsat`.
For one prime, centre profile, fixed direction, and fixed-edge scale it
selects complementary target pairs and their forced auxiliary directions,
but also constructs the *physical* localized-Mobius source.  The model
requires

* one option through every hard target;
* globally distinct auxiliaries with the endpoint Paley-type quota;
* a ternary sum of the Mobius anti-chains;
* the exact requested number of opposite-sign orbit cancellations; and
* the sharp coefficient-support budget forced by the compact triangle
  atoms in every projected row.

At the top branch-C endpoint the cancellation count is one.  This model
represents it by one isolated opposite-sign orbit collision between two
selected paired-SDR options.  That is a restricted subfamily, not an exact
encoding of every one-cancellation ternary sum: in particular, clean
three-half ``2:1`` overlaps exist and lie outside this model.

Within that restricted subfamily, the surviving ``m(p-1)-2`` orbits together
with the prescribed antipodal edge give a simple graph of the required size.
In a hard row its coefficient vector, after subtracting the
complement-literal star, must have l1 norm at most ``3(P_L-3)``.  In an
opposite row the same bound applies without a literal star.  These are
necessary because every remaining atom has exactly three signed triangle
edges.  They are not sufficient for an atom decomposition, so SAT here is
still not a residual-(ii) proof; UNSAT excludes only the recorded centre
profile, scale, fixed direction, and isolated pair-collision subfamily.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from e1_gmin_m4_adaptive_mobius_pairing import (  # noqa: E402
    forced_affine_auxiliary_pair,
)
from e1_gmin_m4_inversion_antisymmetric_radon import (  # noqa: E402
    Edge,
    Point,
    _functional_value,
    _negative_edge,
    edge_radon_image,
    localized_star_trade,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (  # noqa: E402
    paley_direction_sign,
    paley_edge_sign,
)
from e1_gmin_m4_symmetric_fixed_edge_elimination import (  # noqa: E402
    orbit_fixed_word,
)
from io_atomic import write_json_atomic  # noqa: E402
from scripts.residual_branch_c_aux_sdr_cpsat import (  # noqa: E402
    Chart,
    _profile,
    build_chart,
)


Orbit = Edge
TargetKey = tuple[object, ...]


@dataclass(frozen=True)
class PhysicalOption:
    first_target: int
    second_target: int
    first_sign: int
    second_sign: int
    first_auxiliary: int
    second_auxiliary: int
    orbit_coefficients: tuple[tuple[Orbit, int], ...]


def _scaled(p: int, scalar: int, value: Point) -> Point:
    return scalar * value[0] % p, scalar * value[1] % p


def _canonical_orbit_coefficient(p: int, source: dict[Edge, int]) -> dict[Orbit, int]:
    """Compress an antisymmetric edge chain to one coefficient per orbit."""
    out: dict[Orbit, int] = {}
    for edge, value in source.items():
        negative = _negative_edge(p, edge)
        if edge < negative:
            if source.get(negative) != -value:
                raise ArithmeticError("a Mobius source lost antisymmetry")
            out[edge] = value
    return out


def _half_orbits(
    chart: Chart,
    target_coordinate: int,
    alpha: int,
    auxiliary_coordinate: int,
    sign: int,
    c: int,
) -> dict[Orbit, int]:
    """Build the exactly scaled physical Mobius anti-chain for one target."""
    p = chart.p
    # X=L/j has evaluation alpha at x0.  The singleton condition on
    # x=c*x0 is N(x)=2*sign, hence N(x0)=2*sign/c.
    target = _scaled(
        p, alpha, chart.direction_by_coordinate[target_coordinate]
    )
    auxiliary = _scaled(
        p,
        2 * sign * pow(c, -1, p),
        chart.direction_by_coordinate[auxiliary_coordinate],
    )
    source = localized_star_trade(p, target, auxiliary, 1)
    out = _canonical_orbit_coefficient(p, source)
    if len(out) != p - 1 or set(out.values()) - {-1, 1}:
        raise ArithmeticError("a scaled Mobius half changed support")
    return out


def generate_physical_options(
    chart: Chart, alphas: Sequence[int], c: int
) -> tuple[PhysicalOption, ...]:
    """Enumerate all signed complementary pairs without sign deduplication."""
    p = chart.p
    targets = chart.hard_coordinates
    if len(alphas) != len(targets) or any(alpha % p == 0 for alpha in alphas):
        raise ValueError("need one nonzero alpha for every hard target")
    options: list[PhysicalOption] = []
    for first in range(len(targets)):
        for second in range(first + 1, len(targets)):
            for first_sign in (-1, 1):
                for second_sign in (-1, 1):
                    try:
                        forced = forced_affine_auxiliary_pair(
                            p,
                            targets[first],
                            targets[second],
                            int(alphas[first]),
                            int(alphas[second]),
                            c,
                            first_sign,
                            second_sign,
                        )
                    except ValueError:
                        continue
                    first_auxiliary = int(forced["first_auxiliary_coordinate_U"])
                    second_auxiliary = int(forced["second_auxiliary_coordinate_V"])
                    total = Counter(
                        _half_orbits(
                            chart,
                            targets[first],
                            int(alphas[first]),
                            first_auxiliary,
                            first_sign,
                            c,
                        )
                    )
                    total.update(
                        _half_orbits(
                            chart,
                            targets[second],
                            int(alphas[second]),
                            second_auxiliary,
                            second_sign,
                            c,
                        )
                    )
                    total = Counter({edge: value for edge, value in total.items() if value})
                    if any(abs(value) > 1 for value in total.values()):
                        # Same-sign internal collisions cannot occur in a ternary sum.
                        continue
                    options.append(
                        PhysicalOption(
                            first,
                            second,
                            first_sign,
                            second_sign,
                            first_auxiliary,
                            second_auxiliary,
                            tuple(sorted(total.items())),
                        )
                    )
    return tuple(options)


def _projected_key(p: int, direction: Point, edge: Edge) -> TargetKey:
    left = _functional_value(p, direction, edge[0])
    right = _functional_value(p, direction, edge[1])
    if left == right:
        return ("P",)
    first, second = sorted((left, right))
    return ("K", first, second)


def _literal_centres(chart: Chart, alphas: Sequence[int]) -> dict[int, int]:
    """Return the hard literal label in each canonical projective row."""
    p = chart.p
    directions = projective_functionals(p)
    alpha_by_coordinate = dict(zip(chart.hard_coordinates, alphas, strict=True))
    out: dict[int, int] = {}
    for direction_index, direction in enumerate(directions):
        if paley_direction_sign(p, direction) != 1:
            continue
        evaluation = _functional_value(p, direction, chart.x0)
        if evaluation == 0:
            raise ArithmeticError("the opposite fixed direction was classified as hard")
        normalized = _scaled(p, pow(evaluation, -1, p), direction)
        coordinate = chart.direction_by_coordinate.index(normalized)
        alpha = int(alpha_by_coordinate[coordinate]) % p
        out[direction_index] = evaluation * pow(alpha, -1, p) % p
    return out


def _fixed_edge(chart: Chart, c: int) -> Edge:
    p = chart.p
    point = _scaled(p, c, chart.x0)
    edge = tuple(sorted((point, ((-point[0]) % p, (-point[1]) % p))))
    if edge[0] == edge[1]:
        raise ArithmeticError("the forced antipodal edge collapsed")
    return edge  # type: ignore[return-value]


def _graph_from_orbit_values(
    p: int, orbit_values: dict[Orbit, int], fixed_edge: Edge
) -> tuple[Edge, ...]:
    graph: set[Edge] = set()
    for orbit, value in orbit_values.items():
        if value not in (-1, 1):
            raise ArithmeticError("a replayed orbit value is not ternary")
        signed_difference = paley_edge_sign(p, orbit) * value
        graph.add(orbit if signed_difference == 1 else _negative_edge(p, orbit))
    if fixed_edge in graph:
        raise ArithmeticError("the fixed edge was already in a nonfixed orbit")
    graph.add(fixed_edge)
    return tuple(sorted(graph))


def _option_graph(p: int, option: PhysicalOption) -> tuple[Edge, ...]:
    """Return the sixty selected graph edges before cross-option cancellation."""
    graph = []
    for orbit, value in option.orbit_coefficients:
        signed_difference = paley_edge_sign(p, orbit) * value
        graph.append(
            orbit if signed_difference == 1 else _negative_edge(p, orbit)
        )
    if len(set(graph)) != len(graph):
        raise ArithmeticError("one physical option repeated a graph edge")
    return tuple(sorted(graph))


def _fixed_word_capacity_replay(
    chart: Chart,
    alphas: Sequence[int],
    c: int,
    graph: Sequence[Edge],
    orbit_values: dict[Orbit, int],
) -> dict[str, object]:
    """Replay (E.29i) for the actual support and parallel allocation."""
    p = chart.p
    directions = projective_functionals(p)
    literal_centres = _literal_centres(chart, alphas)
    squares = tuple(sorted({value * value % p for value in range(1, p)}))
    square_index = {value: index for index, value in enumerate(squares)}
    block_parity: defaultdict[tuple[int, int], int] = defaultdict(int)
    for orbit in orbit_values:
        record = orbit_fixed_word(p, orbit)
        if not int(record["fixed_word_weight"]):
            continue
        midpoint = tuple(record["midpoint"])
        difference = tuple(record["difference_representative"])
        annihilators = tuple(
            index
            for index, direction in enumerate(directions)
            if _functional_value(p, direction, midpoint) == 0
        )
        if len(annihilators) != 1:
            raise ArithmeticError("a nonzero Phi block lost its direction")
        direction_index = annihilators[0]
        evaluation = _functional_value(
            p, directions[direction_index], difference
        )
        beta = evaluation * evaluation % p
        block_parity[(direction_index, square_index[beta])] ^= 1

    signed_source = {edge: paley_edge_sign(p, edge) for edge in graph}
    image = edge_radon_image(p, signed_source)
    singleton = _scaled(p, c, chart.x0)
    rows = []
    all_feasible = True
    for direction_index, direction in enumerate(directions):
        direction_sign = paley_direction_sign(p, direction)
        parallel = direction_sign * image.get(("P", direction_index), 0)
        bits = [
            block_parity[(direction_index, index)]
            for index in range(len(squares))
        ]
        literal_cell = None
        if direction_sign == 1:
            centre = literal_centres[direction_index]
            literal_cell = square_index[centre * centre % p]
            bits[literal_cell] ^= 1
        singleton_value = _functional_value(p, direction, singleton)
        singleton_cell = None
        if singleton_value:
            singleton_cell = square_index[
                singleton_value * singleton_value % p
            ]
            bits[singleton_cell] ^= 1
        atom_count = parallel - 3
        weight = sum(bits)
        feasible = bool(
            atom_count >= 0
            and weight <= atom_count
            and weight % 2 == atom_count % 2
        )
        all_feasible &= feasible
        rows.append(
            {
                "direction_index": direction_index,
                "direction_sign": direction_sign,
                "parallel_count": parallel,
                "atom_count": atom_count,
                "required_fixed_atom_weight": weight,
                "literal_square_cell": literal_cell,
                "singleton_square_cell": singleton_cell,
                "capacity_and_parity_feasible": feasible,
            }
        )
    return {
        "nonzero_Phi_block_parity_weight": sum(block_parity.values()),
        "rows": rows,
        "all_rows_capacity_and_parity_feasible": all_feasible,
    }


def solve_sparse_scale(
    chart: Chart,
    alphas: Sequence[int],
    c: int,
    cancellations: int,
    seconds: float,
    workers: int,
    random_seed: int,
    enforce_sparse_budget: bool = True,
) -> dict[str, object]:
    """Solve one exact physical-orbit and compact-coefficient model."""
    from ortools.sat.python import cp_model

    p = chart.p
    h = (p - 1) // 2
    m = (p + 1) // 2
    r = (p - 3) // 4
    options = generate_physical_options(chart, alphas, c)
    model = cp_model.CpModel()
    selected = tuple(model.new_bool_var(f"option_{index}") for index in range(len(options)))

    target_rows: list[list[int]] = [[] for _ in range(m)]
    auxiliary_rows: list[list[int]] = [[] for _ in range(p)]
    orbit_rows: defaultdict[Orbit, list[tuple[int, int]]] = defaultdict(list)
    for index, option in enumerate(options):
        target_rows[option.first_target].append(index)
        target_rows[option.second_target].append(index)
        auxiliary_rows[option.first_auxiliary].append(index)
        auxiliary_rows[option.second_auxiliary].append(index)
        for orbit, coefficient in option.orbit_coefficients:
            orbit_rows[orbit].append((index, coefficient))
    if any(not row for row in target_rows):
        return {
            "c": c,
            "status": "INFEASIBLE_BY_EMPTY_TARGET_ROW",
            "option_count": len(options),
        }
    for row in target_rows:
        model.add_exactly_one(selected[index] for index in row)
    for row in auxiliary_rows:
        if row:
            model.add_at_most_one(selected[index] for index in row)
    model.add(
        sum(
            (
                int(chart.type_by_coordinate[option.first_auxiliary] == 1)
                + int(chart.type_by_coordinate[option.second_auxiliary] == 1)
            )
            * selected[index]
            for index, option in enumerate(options)
        )
        == m - 2
    )

    # This restricted top-end model realizes the sole cancellation as one
    # isolated opposite collision shared by two selected options.  This is
    # not equivalent to arbitrary cancellation-one ternary overlap: clean
    # three-half 2:1 overlaps exist and are intentionally outside this model.
    if cancellations != 1:
        raise ValueError("the sparse top-end formulation currently requires one cancellation")
    option_orbits = tuple(dict(option.orbit_coefficients) for option in options)
    collision_pairs: list[tuple[int, int, Orbit, object]] = []
    physical_conflicts = 0
    for first_index in range(len(options)):
        first_option = options[first_index]
        first_aux = {
            first_option.first_auxiliary, first_option.second_auxiliary
        }
        first_targets = {
            first_option.first_target, first_option.second_target
        }
        for second_index in range(first_index + 1, len(options)):
            second_option = options[second_index]
            # These pairs are already mutually exclusive by the exact-cover
            # and auxiliary-SDR constraints, so no collision indicator is needed.
            if first_targets & {
                second_option.first_target, second_option.second_target
            } or first_aux & {
                second_option.first_auxiliary, second_option.second_auxiliary
            }:
                continue
            common = set(option_orbits[first_index]) & set(option_orbits[second_index])
            same = tuple(
                orbit for orbit in common
                if option_orbits[first_index][orbit] == option_orbits[second_index][orbit]
            )
            opposite = tuple(sorted(common - set(same)))
            if same or len(opposite) > 1:
                model.add(selected[first_index] + selected[second_index] <= 1)
                physical_conflicts += 1
            elif len(opposite) == 1:
                both = model.new_bool_var(
                    f"collision_{first_index}_{second_index}"
                )
                model.add(both <= selected[first_index])
                model.add(both <= selected[second_index])
                model.add(both >= selected[first_index] + selected[second_index] - 1)
                collision_pairs.append(
                    (first_index, second_index, opposite[0], both)
                )
    model.add_exactly_one(pair[3] for pair in collision_pairs)
    expected_support = m * (p - 1) - 2

    directions = projective_functionals(p)
    literal_centres = _literal_centres(chart, alphas)
    fixed_edge = _fixed_edge(chart, c)
    fixed_tau = paley_edge_sign(p, fixed_edge)
    cell_terms: defaultdict[tuple[int, TargetKey], list[object]] = defaultdict(list)
    option_graphs = tuple(_option_graph(p, option) for option in options)
    for option_index, graph_edges in enumerate(option_graphs):
        variable = selected[option_index]
        for edge in graph_edges:
            tau = paley_edge_sign(p, edge)
            for direction_index, direction in enumerate(directions):
                coefficient = paley_direction_sign(p, direction) * tau
                key = _projected_key(p, direction, edge)
                cell_terms[(direction_index, key)].append(coefficient * variable)

    # The two options at the unique opposite-sign collision select opposite
    # physical sides of the same inversion orbit.  The top-end graph omits
    # both, so subtract that central pair from every projected cell.
    for _first_index, _second_index, orbit, variable in collision_pairs:
        for edge in (orbit, _negative_edge(p, orbit)):
            tau = paley_edge_sign(p, edge)
            for direction_index, direction in enumerate(directions):
                coefficient = -paley_direction_sign(p, direction) * tau
                key = _projected_key(p, direction, edge)
                cell_terms[(direction_index, key)].append(coefficient * variable)

    fixed_constants: defaultdict[tuple[int, TargetKey], int] = defaultdict(int)
    for direction_index, direction in enumerate(directions):
        coefficient = paley_direction_sign(p, direction) * fixed_tau
        key = _projected_key(p, direction, fixed_edge)
        fixed_constants[(direction_index, key)] += coefficient

    row_l1 = []
    parallel_expressions = []
    for direction_index, direction in enumerate(directions):
        direction_sign = paley_direction_sign(p, direction)
        parallel = (
            sum(cell_terms[(direction_index, ("P",))])
            + fixed_constants[(direction_index, ("P",))]
        )
        parallel_expressions.append(parallel)
        if direction_sign == 1:
            model.add(parallel >= 3)
            literal = literal_centres[direction_index]
        else:
            model.add(parallel >= r + 2)
            literal = None

        absolute_cells = []
        for left in range(p):
            for right in range(left + 1, p):
                key = ("K", left, right)
                expression = (
                    sum(cell_terms[(direction_index, key)])
                    + fixed_constants[(direction_index, key)]
                )
                # Subtracting the hard baseline star, whose coefficient is
                # -1, adds one on every edge incident with the literal label.
                if literal is not None and literal in (left, right):
                    expression += 1
                bound = len(cell_terms[(direction_index, key)]) + 2
                absolute = model.new_int_var(
                    0, max(1, bound), f"abs_{direction_index}_{left}_{right}"
                )
                model.add_abs_equality(absolute, expression)
                absolute_cells.append(absolute)
        total_l1 = model.new_int_var(0, 3 * m * p, f"row_l1_{direction_index}")
        model.add(total_l1 == sum(absolute_cells))
        if enforce_sparse_budget:
            model.add(total_l1 <= 3 * (parallel - 3))
        row_l1.append(total_l1)

    if not enforce_sparse_budget:
        model.minimize(sum(row_l1))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = random_seed
    solver.parameters.log_search_progress = False
    started = time.monotonic()
    status = solver.solve(model)
    wall = time.monotonic() - started
    name = solver.status_name(status)
    result: dict[str, object] = {
        "c": c,
        "status": name,
        "option_count": len(options),
        "candidate_orbit_count": len(orbit_rows),
        "physical_conflict_pair_count": physical_conflicts,
        "single_opposite_collision_pair_count": len(collision_pairs),
        "collision_model_scope": (
            "one isolated opposite collision between two selected paired-SDR "
            "options; no orbit shared by three or more constituent halves"
        ),
        "higher_multiplicity_overlap_models_included": False,
        "required_cancellations": cancellations,
        "required_used_orbits": expected_support,
        "compact_triangle_l1_budget_enforced": enforce_sparse_budget,
        "wall_seconds": wall,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "best_objective_bound": solver.best_objective_bound,
    }
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return result

    chosen_indices = tuple(
        index for index in range(len(options)) if solver.value(selected[index])
    )
    chosen = tuple(options[index] for index in chosen_indices)
    replay_counter: Counter[Orbit] = Counter()
    for option in chosen:
        replay_counter.update(dict(option.orbit_coefficients))
    if any(abs(value) > 1 for value in replay_counter.values()):
        raise ArithmeticError("the pairwise formulation replay lost ternarity")
    replayed_orbits = {
        orbit: value for orbit, value in replay_counter.items() if value
    }
    replay_cancellations = (
        m * (p - 1) - sum(abs(value) for value in replay_counter.values())
    ) // 2
    if replay_cancellations != 1:
        raise ArithmeticError("the pairwise formulation replay lost its unique cancellation")
    graph = _graph_from_orbit_values(p, replayed_orbits, fixed_edge)
    signed_source = {edge: paley_edge_sign(p, edge) for edge in graph}
    image = edge_radon_image(p, signed_source)
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
    if len(graph) != expected_support + 1:
        raise ArithmeticError("the replayed top-end graph has wrong size")
    if enforce_sparse_budget and not all(row["within_budget"] for row in replay_rows):
        raise ArithmeticError("a solver witness failed the exact l1 replay")
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
        "graph_edges": [
            [list(edge[0]), list(edge[1])] for edge in graph
        ],
        "row_replay": replay_rows,
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
    scales = tuple(
        value
        for value in range(1, args.p)
        if (value - 1) % args.shard_count == args.shard_index
    )
    records = []
    started = time.monotonic()
    for c in scales:
        row = solve_sparse_scale(
            chart,
            alphas,
            c,
            args.cancellations,
            args.seconds_per_scale,
            args.workers,
            args.solver_seed + c,
            enforce_sparse_budget=not args.relax_sparse_budget,
        )
        records.append(row)
        if args.stop_on_sat and row["status"] in ("OPTIMAL", "FEASIBLE"):
            break
    alpha_bytes = b"".join(
        int(value).to_bytes(4, byteorder="little", signed=False)
        for value in alphas
    )
    return {
        "schema": "residual_branch_c_top_sparse_sdr_cpsat_v1",
        "scope": (
            "one centre profile and fixed direction in the paired-SDR isolated "
            "two-option collision subfamily; physical Mobius ternarity and "
            "necessary compact-triangle l1 budgets"
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
        "sat_count": sum(
            row["status"] in ("OPTIMAL", "FEASIBLE") for row in records
        ),
        "infeasible_count": sum(
            str(row["status"]).startswith("INFEASIBLE") for row in records
        ),
        "unknown_count": sum(
            row["status"] not in ("OPTIMAL", "FEASIBLE")
            and not str(row["status"]).startswith("INFEASIBLE")
            for row in records
        ),
        "records": records,
        "elapsed_seconds": time.monotonic() - started,
        "interpretation": (
            "SAT passes a necessary nonfixed-cell sparsity gate but still needs exact "
            "triangle-atom decomposition. UNSAT excludes only the recorded "
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
    parser.add_argument("--solver-seed", type=int, default=15766)
    parser.add_argument("--cancellations", type=int, default=1)
    parser.add_argument("--seconds-per-scale", type=float, default=60.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--stop-on-sat", action="store_true")
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
            f"SAT={result['sat_count']} infeasible={result['infeasible_count']} "
            f"unknown={result['unknown_count']}"
        )


if __name__ == "__main__":
    main()
