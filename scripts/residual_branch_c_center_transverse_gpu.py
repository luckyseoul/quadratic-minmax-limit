#!/usr/bin/env python3
"""GPU coordinate search over the exact p=31 physical centre fibre.

This is a bounded construction search, not an exhaustive certificate.  It
keeps the sixteen scaled ``(L,M)`` halves from the direct parallel design,
fixes its proved clean collision between halves 2 and 13, and varies the
remaining fourteen nonzero centres plus the antipodal fixed edge.  Every
candidate is evaluated on all ``32*binom(31,2)`` literal-corrected transverse
Radon cells.  Pairwise orbit intersections other than the fixed cancelling
pair are rejected before scoring.

The tensor identity is exact: each half/centre table contains the selected
physical graph edges and its own hard literal-star correction.  The two
physical edges removed by the cancelling orbit are subtracted with
multiplicity, and the fixed edge is added.  Reported extremizers are replayed
from their actual 479 graph edges with the repository integer Radon map.

Passing the compact-atom l1 budgets is necessary, not sufficient.  This
program therefore searches for a serious next-stage candidate but cannot by
itself prove residual (ii).
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations, product
import json
import os
from pathlib import Path
import random
import socket
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from e1_gmin_m4_inversion_antisymmetric_radon import (  # noqa: E402
    Edge,
    _functional_value,
    _negative_edge,
    edge_radon_image,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (  # noqa: E402
    paley_direction_sign,
    paley_edge_sign,
)
from e1_gmin_m4_p31_direct_mobius_parallel_design import (  # noqa: E402
    CANCELLATION_DIRECTION_INDEX,
    HALVES,
    P,
    PHYSICAL_CENTERS,
    PHYSICAL_FIXED_POINT,
    _canonical_center,
    _oriented_orbit_coefficients,
    _projective_index,
    centered_physical_graph,
    centered_physical_parallel_design_certificate,
    transverse_compact_l1_diagnostic,
)
from io_atomic import write_json_atomic  # noqa: E402


PAIR_FIRST = 2
PAIR_SECOND = 13
PAIR_FIRST_CENTER = 28
PAIR_SECOND_CENTER = 1


def _cell_catalog() -> tuple[tuple[tuple[int, int], ...], dict[tuple[int, int], int]]:
    cells = tuple(combinations(range(P), 2))
    return cells, {cell: index for index, cell in enumerate(cells)}


def _edge_projection(
    edge: Edge,
    directions: tuple[tuple[int, int], ...],
    direction_signs: tuple[int, ...],
    cell_index: dict[tuple[int, int], int],
) -> np.ndarray:
    out = np.zeros((P + 1, len(cell_index)), dtype=np.int16)
    tau = paley_edge_sign(P, edge)
    for direction_index, direction in enumerate(directions):
        first = _functional_value(P, direction, edge[0])
        second = _functional_value(P, direction, edge[1])
        if first == second:
            continue
        cell = tuple(sorted((first, second)))
        out[direction_index, cell_index[cell]] += direction_signs[direction_index] * tau
    return out


def _fixed_edges(
    directions: tuple[tuple[int, int], ...],
) -> tuple[Edge, ...]:
    fixed = set()
    target = directions[CANCELLATION_DIRECTION_INDEX]
    for point in product(range(P), repeat=2):
        if point == (0, 0) or _functional_value(P, target, point) != 0:
            continue
        negative = (-point[0] % P, -point[1] % P)
        fixed.add(tuple(sorted((point, negative))))
    if len(fixed) != (P - 1) // 2:
        raise ArithmeticError("the fixed-edge fibre changed")
    return tuple(sorted(fixed))


def build_exact_tables() -> dict[str, object]:
    directions = projective_functionals(P)
    direction_signs = tuple(
        paley_direction_sign(P, direction) for direction in directions
    )
    cells, cell_index = _cell_catalog()
    orbit_cache = tuple(
        tuple(
            _oriented_orbit_coefficients(direction, auxiliary, center)
            for center in range(1, P)
        )
        for direction, auxiliary in HALVES
    )

    half_tables = np.zeros(
        (len(HALVES), P - 1, P + 1, len(cells)), dtype=np.int16
    )
    for half_index, ((target, _auxiliary), center_rows) in enumerate(
        zip(HALVES, orbit_cache, strict=True)
    ):
        target_index = _projective_index(target)
        for center_index, orbit_map in enumerate(center_rows):
            row = half_tables[half_index, center_index]
            for orbit, value in orbit_map.items():
                edge = orbit if value == 1 else _negative_edge(P, orbit)
                row += _edge_projection(
                    edge, directions, direction_signs, cell_index
                )
            canonical_center = _canonical_center(target, center_index + 1)[1]
            for other in range(P):
                if other == canonical_center:
                    continue
                row[
                    target_index,
                    cell_index[tuple(sorted((canonical_center, other)))],
                ] += 1

    first_map = orbit_cache[PAIR_FIRST][PAIR_FIRST_CENTER - 1]
    second_map = orbit_cache[PAIR_SECOND][PAIR_SECOND_CENTER - 1]
    common = set(first_map) & set(second_map)
    if len(common) != 1:
        raise ArithmeticError("the frozen collision pair changed")
    collision_orbit = next(iter(common))
    if first_map[collision_orbit] != -second_map[collision_orbit]:
        raise ArithmeticError("the frozen common orbit stopped cancelling")
    collision_edges = (
        collision_orbit
        if first_map[collision_orbit] == 1
        else _negative_edge(P, collision_orbit),
        collision_orbit
        if second_map[collision_orbit] == 1
        else _negative_edge(P, collision_orbit),
    )
    collision_table = sum(
        (
            _edge_projection(edge, directions, direction_signs, cell_index)
            for edge in collision_edges
        ),
        start=np.zeros((P + 1, len(cells)), dtype=np.int16),
    )

    fixed_edges = _fixed_edges(directions)
    fixed_tables = np.stack(
        [
            _edge_projection(edge, directions, direction_signs, cell_index)
            for edge in fixed_edges
        ]
    )

    pair_bad = np.zeros(
        (len(HALVES), len(HALVES), P - 1, P - 1), dtype=np.bool_
    )
    for first in range(len(HALVES)):
        for second in range(first + 1, len(HALVES)):
            for first_center in range(P - 1):
                first_support = set(orbit_cache[first][first_center])
                for second_center in range(P - 1):
                    bad = bool(
                        first_support & set(orbit_cache[second][second_center])
                    )
                    pair_bad[first, second, first_center, second_center] = bad
                    pair_bad[second, first, second_center, first_center] = bad

    profile = centered_physical_parallel_design_certificate()[
        "final_parallel_profile"
    ]
    budgets = np.asarray([3 * (int(value) - 3) for value in profile], dtype=np.int16)
    return {
        "directions": directions,
        "direction_signs": direction_signs,
        "cells": cells,
        "cell_index": cell_index,
        "orbit_cache": orbit_cache,
        "half_tables": half_tables,
        "collision_orbit": collision_orbit,
        "collision_edges": collision_edges,
        "collision_table": collision_table,
        "fixed_edges": fixed_edges,
        "fixed_tables": fixed_tables,
        "pair_bad": pair_bad,
        "budgets": budgets,
    }


def _valid_assignment(choices: np.ndarray, pair_bad: np.ndarray) -> bool:
    for first in range(len(HALVES)):
        for second in range(first + 1, len(HALVES)):
            if (first, second) == (PAIR_FIRST, PAIR_SECOND):
                continue
            if pair_bad[first, second, choices[first], choices[second]]:
                return False
    return True


def _random_assignments(
    count: int, pair_bad: np.ndarray, seed: int
) -> np.ndarray:
    rng = random.Random(seed)
    out = []
    frozen = {
        PAIR_FIRST: PAIR_FIRST_CENTER - 1,
        PAIR_SECOND: PAIR_SECOND_CENTER - 1,
    }
    baseline = np.asarray([center - 1 for center in PHYSICAL_CENTERS], dtype=np.int16)
    if not _valid_assignment(baseline, pair_bad):
        raise ArithmeticError("the published center witness is no longer valid")
    out.append(baseline)
    while len(out) < count:
        choices = np.full(len(HALVES), -1, dtype=np.int16)
        for index, value in frozen.items():
            choices[index] = value
        order = [
            index for index in range(len(HALVES)) if index not in frozen
        ]
        rng.shuffle(order)
        valid = True
        for index in order:
            options = list(range(P - 1))
            rng.shuffle(options)
            selected = None
            for option in options:
                if all(
                    other == index
                    or choices[other] < 0
                    or not pair_bad[index, other, option, choices[other]]
                    for other in range(len(HALVES))
                ):
                    selected = option
                    break
            if selected is None:
                valid = False
                break
            choices[index] = selected
        if valid and _valid_assignment(choices, pair_bad):
            out.append(choices)
    return np.stack(out)


def _cpu_tensor(
    choices: np.ndarray, fixed_index: int, tables: dict[str, object]
) -> np.ndarray:
    half_tables = tables["half_tables"]
    assert isinstance(half_tables, np.ndarray)
    fixed_tables = tables["fixed_tables"]
    collision_table = tables["collision_table"]
    assert isinstance(fixed_tables, np.ndarray)
    assert isinstance(collision_table, np.ndarray)
    out = fixed_tables[fixed_index].copy() - collision_table
    for half_index, center_index in enumerate(choices):
        out += half_tables[half_index, int(center_index)]
    return out


def _score_tensor(tensor: np.ndarray, budgets: np.ndarray) -> tuple[int, int, tuple[int, ...]]:
    row_l1 = tuple(int(value) for value in np.abs(tensor).sum(axis=1))
    excess = tuple(max(0, value - int(budget)) for value, budget in zip(row_l1, budgets, strict=True))
    return sum(excess), max(excess), row_l1


def validate_published_witness(tables: dict[str, object]) -> dict[str, object]:
    fixed_edges = tables["fixed_edges"]
    assert isinstance(fixed_edges, tuple)
    published_fixed = tuple(
        sorted(
            (
                PHYSICAL_FIXED_POINT,
                (-PHYSICAL_FIXED_POINT[0] % P, -PHYSICAL_FIXED_POINT[1] % P),
            )
        )
    )
    fixed_index = fixed_edges.index(published_fixed)
    choices = np.asarray([center - 1 for center in PHYSICAL_CENTERS], dtype=np.int16)
    tensor = _cpu_tensor(choices, fixed_index, tables)
    budgets = tables["budgets"]
    assert isinstance(budgets, np.ndarray)
    total, maximum, row_l1 = _score_tensor(tensor, budgets)
    frozen = transverse_compact_l1_diagnostic()
    expected = tuple(
        int(row["transverse_residual_l1"]) for row in frozen["rows"]
    )
    if row_l1 != expected or total != 5068 or maximum != 194:
        raise ArithmeticError("the additive tensor failed the public graph replay")
    return {
        "centers": list(PHYSICAL_CENTERS),
        "fixed_edge_index": fixed_index,
        "fixed_edge": [list(point) for point in published_fixed],
        "total_positive_l1_excess": total,
        "maximum_l1_excess": maximum,
        "row_l1": list(row_l1),
        "exact_public_graph_replay": True,
    }


def _candidate_valid_mask(
    choices: np.ndarray,
    variable: int,
    option: int,
    pair_bad: np.ndarray,
) -> np.ndarray:
    valid = np.ones(len(choices), dtype=np.bool_)
    for other in range(len(HALVES)):
        if other == variable:
            continue
        valid &= ~pair_bad[
            variable, other, option, choices[:, other]
        ]
    return valid


def gpu_coordinate_search(
    tables: dict[str, object],
    chains: int,
    restarts: int,
    sweeps: int,
    seed: int,
    primary: str,
) -> dict[str, object]:
    # This host's CUDA headers are newer than its NVRTC; disabling the optional
    # CUB reducer uses CuPy's compatible built-in reduction kernel.
    os.environ.setdefault("CUPY_ACCELERATORS", "")
    import cupy as cp

    half_tables = tables["half_tables"]
    fixed_tables = tables["fixed_tables"]
    collision_table = tables["collision_table"]
    pair_bad = tables["pair_bad"]
    budgets = tables["budgets"]
    assert isinstance(half_tables, np.ndarray)
    assert isinstance(fixed_tables, np.ndarray)
    assert isinstance(collision_table, np.ndarray)
    assert isinstance(pair_bad, np.ndarray)
    assert isinstance(budgets, np.ndarray)

    half_gpu = cp.asarray(half_tables)
    fixed_gpu = cp.asarray(fixed_tables)
    collision_gpu = cp.asarray(collision_table)
    budgets_gpu = cp.asarray(budgets)[None, :]
    rng = random.Random(seed)
    mutable = [
        index
        for index in range(len(HALVES))
        if index not in (PAIR_FIRST, PAIR_SECOND)
    ]
    evaluated = 0
    best: dict[str, object] | None = None
    started = time.perf_counter()

    for restart in range(restarts):
        choices = _random_assignments(
            chains, pair_bad, rng.randrange(1 << 62)
        )
        fixed_choice = np.asarray(
            [rng.randrange(len(fixed_tables)) for _ in range(chains)],
            dtype=np.int16,
        )
        if restart == 0:
            published = validate_published_witness(tables)
            fixed_choice[0] = int(published["fixed_edge_index"])

        state = fixed_gpu[cp.asarray(fixed_choice)] - collision_gpu[None, :, :]
        for half_index in range(len(HALVES)):
            state += half_gpu[
                half_index, cp.asarray(choices[:, half_index])
            ]

        for _sweep in range(sweeps):
            changed = False
            rng.shuffle(mutable)
            for variable in mutable:
                old = half_gpu[
                    variable, cp.asarray(choices[:, variable])
                ]
                base = state - old
                best_total = np.full(chains, np.iinfo(np.int64).max, dtype=np.int64)
                best_max = np.full(chains, np.iinfo(np.int64).max, dtype=np.int64)
                best_option = choices[:, variable].copy()
                for option in range(P - 1):
                    valid = _candidate_valid_mask(
                        choices, variable, option, pair_bad
                    )
                    if not valid.any():
                        continue
                    candidate = base + half_gpu[variable, option]
                    l1 = cp.abs(candidate).sum(axis=2)
                    excess = cp.maximum(l1 - budgets_gpu, 0)
                    total_score = cp.asnumpy(excess.sum(axis=1))
                    max_score = cp.asnumpy(excess.max(axis=1))
                    total_score[~valid] = np.iinfo(np.int64).max
                    max_score[~valid] = np.iinfo(np.int64).max
                    if primary == "total":
                        improve = (total_score < best_total) | (
                            (total_score == best_total) & (max_score < best_max)
                        )
                    else:
                        improve = (max_score < best_max) | (
                            (max_score == best_max) & (total_score < best_total)
                        )
                    best_total[improve] = total_score[improve]
                    best_max[improve] = max_score[improve]
                    best_option[improve] = option
                    evaluated += int(valid.sum())
                if np.any(best_option != choices[:, variable]):
                    changed = True
                choices[:, variable] = best_option
                state = base + half_gpu[
                    variable, cp.asarray(best_option)
                ]

            base = state - fixed_gpu[cp.asarray(fixed_choice)]
            best_total = np.full(chains, np.iinfo(np.int64).max, dtype=np.int64)
            best_max = np.full(chains, np.iinfo(np.int64).max, dtype=np.int64)
            best_fixed = fixed_choice.copy()
            for option in range(len(fixed_tables)):
                candidate = base + fixed_gpu[option]
                l1 = cp.abs(candidate).sum(axis=2)
                excess = cp.maximum(l1 - budgets_gpu, 0)
                total_score = cp.asnumpy(excess.sum(axis=1))
                max_score = cp.asnumpy(excess.max(axis=1))
                if primary == "total":
                    improve = (total_score < best_total) | (
                        (total_score == best_total) & (max_score < best_max)
                    )
                else:
                    improve = (max_score < best_max) | (
                        (max_score == best_max) & (total_score < best_total)
                    )
                best_total[improve] = total_score[improve]
                best_max[improve] = max_score[improve]
                best_fixed[improve] = option
                evaluated += chains
            if np.any(best_fixed != fixed_choice):
                changed = True
            fixed_choice = best_fixed
            state = base + fixed_gpu[cp.asarray(fixed_choice)]
            if not changed:
                break

        l1 = cp.asnumpy(cp.abs(state).sum(axis=2))
        excess = np.maximum(l1 - budgets[None, :], 0)
        totals = excess.sum(axis=1)
        maxima = excess.max(axis=1)
        winner = min(
            range(chains),
            key=lambda index: (
                (int(totals[index]), int(maxima[index]))
                if primary == "total"
                else (int(maxima[index]), int(totals[index]))
            ),
        )
        record = {
            "centers": [int(value) + 1 for value in choices[winner]],
            "fixed_edge_index": int(fixed_choice[winner]),
            "total_positive_l1_excess": int(totals[winner]),
            "maximum_l1_excess": int(maxima[winner]),
            "row_l1": [int(value) for value in l1[winner]],
            "restart": restart,
        }
        record_key = (
            (int(record["total_positive_l1_excess"]), int(record["maximum_l1_excess"]))
            if primary == "total"
            else (int(record["maximum_l1_excess"]), int(record["total_positive_l1_excess"]))
        )
        best_key = None if best is None else (
            (int(best["total_positive_l1_excess"]), int(best["maximum_l1_excess"]))
            if primary == "total"
            else (int(best["maximum_l1_excess"]), int(best["total_positive_l1_excess"]))
        )
        if best_key is None or record_key < best_key:
            best = record
            print(json.dumps({"new_best": best}), flush=True)
        if best["total_positive_l1_excess"] == 0:
            break

    cp.cuda.Stream.null.synchronize()
    elapsed = time.perf_counter() - started
    assert best is not None
    return {
        "gpu_device": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        "elapsed_seconds": elapsed,
        "objective_states_evaluated": evaluated,
        "objective_states_per_second": evaluated / elapsed,
        "chains": chains,
        "restarts_completed": restart + 1,
        "maximum_sweeps": sweeps,
        "primary_objective": primary,
        "best": best,
    }


def exact_graph_replay(
    centers: list[int], fixed_index: int, tables: dict[str, object]
) -> dict[str, object]:
    orbit_cache = tables["orbit_cache"]
    fixed_edges = tables["fixed_edges"]
    directions = tables["directions"]
    direction_signs = tables["direction_signs"]
    budgets = tables["budgets"]
    assert isinstance(orbit_cache, tuple)
    assert isinstance(fixed_edges, tuple)
    assert isinstance(directions, tuple)
    assert isinstance(direction_signs, tuple)
    assert isinstance(budgets, np.ndarray)
    choices = np.asarray([center - 1 for center in centers], dtype=np.int16)
    pair_bad = tables["pair_bad"]
    assert isinstance(pair_bad, np.ndarray)
    if not _valid_assignment(choices, pair_bad):
        raise ArithmeticError("reported candidate acquired an extra overlap")

    total: Counter[Edge] = Counter()
    for half_index, center in enumerate(centers):
        total.update(orbit_cache[half_index][center - 1])
    vanished = tuple(sorted(orbit for orbit, value in total.items() if value == 0))
    surviving = {orbit: value for orbit, value in total.items() if value}
    if len(vanished) != 1 or len(surviving) != 478 or any(
        abs(value) != 1 for value in surviving.values()
    ):
        raise ArithmeticError("reported candidate is not the one-cancellation ternary support")
    graph = [
        orbit if value == 1 else _negative_edge(P, orbit)
        for orbit, value in surviving.items()
    ]
    fixed_edge = fixed_edges[fixed_index]
    graph.append(fixed_edge)
    graph = sorted(graph)
    if len(graph) != 479 or len(set(graph)) != 479:
        raise ArithmeticError("reported candidate is not a 479-edge simple graph")
    source = {edge: paley_edge_sign(P, edge) for edge in graph}
    image = edge_radon_image(P, source)
    row_l1 = []
    parallel = []
    for direction_index, direction in enumerate(directions):
        sign = direction_signs[direction_index]
        parallel.append(sign * image.get(("P", direction_index), 0))
        coefficients = []
        literal = None
        if sign == 1:
            half_index = next(
                index
                for index, (target, _auxiliary) in enumerate(HALVES)
                if _projective_index(target) == direction_index
            )
            literal = _canonical_center(HALVES[half_index][0], centers[half_index])[1]
        for left, right in combinations(range(P), 2):
            value = sign * image.get(
                ("K", direction_index, left, right), 0
            )
            if literal is not None and literal in (left, right):
                value += 1
            coefficients.append(value)
        row_l1.append(sum(abs(value) for value in coefficients))
    excess = [
        max(0, value - int(budget))
        for value, budget in zip(row_l1, budgets, strict=True)
    ]
    graph_bytes = json.dumps(graph, separators=(",", ":")).encode()
    import hashlib

    return {
        "centers": centers,
        "fixed_edge": [list(point) for point in fixed_edge],
        "cancelled_orbit": [list(point) for point in vanished[0]],
        "graph_edge_count": len(graph),
        "graph_sha256": hashlib.sha256(graph_bytes).hexdigest(),
        "parallel_profile": parallel,
        "row_l1": row_l1,
        "total_positive_l1_excess": sum(excess),
        "maximum_l1_excess": max(excess),
        "passes_all_l1_rows": not any(excess),
        "exact_integer_graph_replay": True,
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    tables = build_exact_tables()
    self_test = validate_published_witness(tables)
    search = gpu_coordinate_search(
        tables,
        chains=args.chains,
        restarts=args.restarts,
        sweeps=args.sweeps,
        seed=args.seed,
        primary=args.primary,
    )
    best = search["best"]
    assert isinstance(best, dict)
    replay = exact_graph_replay(
        list(best["centers"]), int(best["fixed_edge_index"]), tables
    )
    if (
        replay["row_l1"] != best["row_l1"]
        or replay["total_positive_l1_excess"]
        != best["total_positive_l1_excess"]
        or replay["maximum_l1_excess"] != best["maximum_l1_excess"]
    ):
        raise ArithmeticError("GPU extremizer failed exact CPU replay")
    result = {
        "schema": "residual_branch_c_center_transverse_gpu_v1",
        "classification": (
            "bounded exact-objective center search in one fixed p31 half/auxiliary fibre"
        ),
        "host": socket.gethostname(),
        "p": P,
        "fixed_collision_halves": [PAIR_FIRST, PAIR_SECOND],
        "fixed_collision_centers": [PAIR_FIRST_CENTER, PAIR_SECOND_CENTER],
        "collision_orbit": [list(point) for point in tables["collision_orbit"]],
        "published_graph_self_test": self_test,
        "search": search,
        "best_exact_replay": replay,
        "l1_pass_is_sufficient_for_atom_decomposition": False,
        "fibre_exhausted": False,
        "residual_ii_closed": False,
    }
    if args.output:
        write_json_atomic(args.output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chains", type=int, default=64)
    parser.add_argument("--restarts", type=int, default=4)
    parser.add_argument("--sweeps", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--primary", choices=("total", "max"), default="total")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.chains < 1 or args.restarts < 1 or args.sweeps < 1:
        parser.error("chains, restarts, and sweeps must be positive")
    result = run(args)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
