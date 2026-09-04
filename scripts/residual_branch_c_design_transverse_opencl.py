#!/usr/bin/env python3
"""Search exact p=31 moved auxiliary designs with the OpenCL row objective.

The input is a checked list of sixteen-half designs and clean collision
seeds.  This runner rebuilds all physical centre tables from the displayed
functionals, verifies the parallel profile and unique opposite collision,
then calls the independent OpenCL coordinate backend.  It supports parity
sharding by ``design_index`` so NUKA and a CUDA worker can cover disjoint
quotient fibres.

This remains a bounded search in supplied Möbius-design fibres.  Failure to
find a compact-l1 pass is not an exhaustive exclusion of those fibres or a
closure of residual (ii).
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import socket
import sys
import time
from typing import Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from e1_gmin_m4_inversion_antisymmetric_radon import (  # noqa: E402
    Edge,
    Functional,
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
    P,
    _canonical_center,
    _oriented_orbit_coefficients,
    _physical_parallel_profile,
    _projective_index,
    _spatial_direction_index,
)
from io_atomic import write_json_atomic  # noqa: E402
from scripts.residual_branch_c_center_transverse_gpu import (  # noqa: E402
    _cell_catalog,
    _edge_projection,
)
from scripts.residual_branch_c_center_transverse_opencl import (  # noqa: E402
    _valid_assignment,
    opencl_coordinate_search,
)


def _functional(value: Sequence[int]) -> Functional:
    if len(value) != 2:
        raise ValueError("each functional must have two coordinates")
    out = int(value[0]) % P, int(value[1]) % P
    if out == (0, 0):
        raise ValueError("a zero functional is invalid")
    return out


def _halves(value: object) -> tuple[tuple[Functional, Functional], ...]:
    if not isinstance(value, list):
        raise ValueError("design halves must be a list")
    out = []
    for row in value:
        if not isinstance(row, list) or len(row) != 2:
            raise ValueError("each half must display [target, auxiliary]")
        out.append((_functional(row[0]), _functional(row[1])))
    if len(out) != 16:
        raise ValueError("a p31 top design must contain sixteen halves")
    return tuple(out)


def load_designs(path: Path) -> tuple[dict[str, object], ...]:
    payload = json.loads(path.read_text())
    rows = payload if isinstance(payload, list) else payload.get("designs")
    if rows is None and isinstance(payload, dict):
        rows = payload.get("searchable_designs")
    if not isinstance(rows, list):
        raise ValueError("design JSON needs a top-level designs list")
    out = []
    seen = set()
    for offset, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError("each design record must be an object")
        design_index = int(
            row.get(
                "design_index",
                row.get("component_id", row.get("move_id", offset)),
            )
        )
        if design_index in seen:
            raise ValueError("design_index must be unique")
        seen.add(design_index)
        normalized = dict(row)
        normalized["design_index"] = design_index
        normalized["halves"] = _halves(row.get("halves"))
        seeds = row.get("clean_collision_seeds")
        if seeds is None and row.get("canonical_collision_seed") is not None:
            seeds = [row["canonical_collision_seed"]]
        if not isinstance(seeds, list):
            raise ValueError("clean_collision_seeds must be a list")
        normalized["clean_collision_seeds"] = tuple(seeds)
        normalized.setdefault(
            "clean_collision_seed_count",
            30 if row.get("canonical_collision_seed") is not None else len(seeds),
        )
        out.append(normalized)
    return tuple(out)


def _canonical_seed(design: dict[str, object]) -> dict[str, object]:
    seeds = design["clean_collision_seeds"]
    assert isinstance(seeds, tuple)

    def key(seed: object) -> tuple[object, ...]:
        if not isinstance(seed, dict):
            raise ValueError("a collision seed must be an object")
        return (
            tuple(int(value) for value in seed["half_indices"]),
            tuple(int(value) for value in seed["centers"]),
            tuple(tuple(int(x) for x in point) for point in seed["orbit"]),
        )

    seed = min(seeds, key=key)
    assert isinstance(seed, dict)
    return seed


def _fixed_edges(direction_index: int) -> tuple[Edge, ...]:
    directions = projective_functionals(P)
    target = directions[direction_index]
    fixed = set()
    for point in product(range(P), repeat=2):
        if point == (0, 0) or _functional_value(P, target, point) != 0:
            continue
        negative = -point[0] % P, -point[1] % P
        fixed.add(tuple(sorted((point, negative))))
    out = tuple(sorted(fixed))
    if len(out) != (P - 1) // 2 or any(
        _spatial_direction_index(edge) != direction_index for edge in out
    ):
        raise ArithmeticError("the requested antipodal fixed fibre is invalid")
    return out


def build_design_tables(
    design: dict[str, object], seed: dict[str, object]
) -> dict[str, object]:
    halves = design["halves"]
    assert isinstance(halves, tuple)
    half_indices = tuple(int(value) for value in seed["half_indices"])
    centers = tuple(int(value) for value in seed["centers"])
    if len(half_indices) != 2 or len(set(half_indices)) != 2:
        raise ValueError("this backend requires a two-half collision seed")
    if len(centers) != 2 or any(not 1 <= value < P for value in centers):
        raise ValueError("collision centers must be two nonzero field values")
    if any(not 0 <= value < len(halves) for value in half_indices):
        raise ValueError("collision half index is out of range")

    directions = projective_functionals(P)
    direction_signs = tuple(
        paley_direction_sign(P, direction) for direction in directions
    )
    hard = tuple(index for index, sign in enumerate(direction_signs) if sign == 1)
    opposite = tuple(index for index, sign in enumerate(direction_signs) if sign == -1)
    target_indices = tuple(_projective_index(target) for target, _aux in halves)
    auxiliary_indices = tuple(_projective_index(aux) for _target, aux in halves)
    if set(target_indices) != set(hard) or len(set(target_indices)) != 16:
        raise ArithmeticError("the moved design lost its hard-target SDR")
    if len(set(auxiliary_indices)) != 16:
        raise ArithmeticError("the moved design lost its auxiliary SDR")

    cells, cell_index = _cell_catalog()
    orbit_cache = tuple(
        tuple(
            _oriented_orbit_coefficients(direction, auxiliary, center)
            for center in range(1, P)
        )
        for direction, auxiliary in halves
    )
    first_map = orbit_cache[half_indices[0]][centers[0] - 1]
    second_map = orbit_cache[half_indices[1]][centers[1] - 1]
    common = set(first_map) & set(second_map)
    if len(common) != 1:
        raise ArithmeticError("the displayed seed is not a unique pair collision")
    collision_orbit = next(iter(common))
    coefficients = first_map[collision_orbit], second_map[collision_orbit]
    if coefficients[0] != -coefficients[1]:
        raise ArithmeticError("the displayed collision does not cancel")
    stated_orbit = tuple(
        tuple(int(value) for value in point) for point in seed["orbit"]
    )
    stated_coefficients = tuple(int(value) for value in seed["coefficients"])
    if collision_orbit != stated_orbit or coefficients != stated_coefficients:
        raise ArithmeticError("collision seed metadata failed exact replay")
    collision_direction = _spatial_direction_index(collision_orbit)
    fixed_direction = int(
        design.get("fixed_direction_index", seed.get("spatial_direction_index"))
    )
    if collision_direction != fixed_direction:
        raise ArithmeticError("collision and fixed directions differ")

    half_tables = np.zeros(
        (len(halves), P - 1, P + 1, len(cells)), dtype=np.int16
    )
    for half_index, (((target, _auxiliary), center_rows)) in enumerate(
        zip(halves, orbit_cache, strict=True)
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

    collision_edges = (
        collision_orbit
        if coefficients[0] == 1
        else _negative_edge(P, collision_orbit),
        collision_orbit
        if coefficients[1] == 1
        else _negative_edge(P, collision_orbit),
    )
    collision_table = sum(
        (
            _edge_projection(edge, directions, direction_signs, cell_index)
            for edge in collision_edges
        ),
        start=np.zeros((P + 1, len(cells)), dtype=np.int16),
    )
    fixed_edges = _fixed_edges(fixed_direction)
    fixed_tables = np.stack(
        [
            _edge_projection(edge, directions, direction_signs, cell_index)
            for edge in fixed_edges
        ]
    )

    pair_bad = np.zeros(
        (len(halves), len(halves), P - 1, P - 1), dtype=np.bool_
    )
    supports = tuple(
        tuple(set(orbit_map) for orbit_map in center_rows)
        for center_rows in orbit_cache
    )
    for first in range(len(halves)):
        for second in range(first + 1, len(halves)):
            for first_center in range(P - 1):
                for second_center in range(P - 1):
                    bad = bool(
                        supports[first][first_center]
                        & supports[second][second_center]
                    )
                    pair_bad[first, second, first_center, second_center] = bad
                    pair_bad[second, first, second_center, first_center] = bad

    frozen_assignment = np.full(len(halves), -1, dtype=np.int16)
    for index, center in zip(half_indices, centers, strict=True):
        frozen_assignment[index] = center - 1
    # At the displayed pair itself there is exactly one shared orbit.  All
    # other half pairs are screened dynamically as their centers are chosen.
    if not pair_bad[
        half_indices[0], half_indices[1], centers[0] - 1, centers[1] - 1
    ]:
        raise ArithmeticError("the pair-intersection table lost the seed")

    raw = [0] * (P + 1)
    for target, auxiliary in halves:
        profile = _physical_parallel_profile(target, auxiliary)
        for index, value in enumerate(profile):
            raw[index] += value
    stated_raw = tuple(int(value) for value in design["raw_parallel_profile"])
    if tuple(raw) != stated_raw:
        raise ArithmeticError("the moved design raw parallel profile changed")
    final = list(raw)
    final[collision_direction] -= 2
    final[fixed_direction] += 1
    if (
        sorted(final[index] for index in hard) != [14] * 14 + [15] * 2
        or sorted(final[index] for index in opposite) != [15] * 3 + [16] * 13
    ):
        raise ArithmeticError("the moved design does not reach the top quotas")
    budgets = np.asarray([3 * (value - 3) for value in final], dtype=np.int16)
    return {
        "halves": halves,
        "directions": directions,
        "direction_signs": direction_signs,
        "cells": cells,
        "cell_index": cell_index,
        "orbit_cache": orbit_cache,
        "half_tables": half_tables,
        "collision_orbit": collision_orbit,
        "collision_coefficients": coefficients,
        "collision_edges": collision_edges,
        "collision_table": collision_table,
        "fixed_direction_index": fixed_direction,
        "fixed_edges": fixed_edges,
        "fixed_tables": fixed_tables,
        "pair_bad": pair_bad,
        "raw_parallel_profile": tuple(raw),
        "final_parallel_profile": tuple(final),
        "budgets": budgets,
    }


def exact_design_graph_replay(
    centers: Sequence[int], fixed_index: int, tables: dict[str, object]
) -> dict[str, object]:
    halves = tables["halves"]
    orbit_cache = tables["orbit_cache"]
    fixed_edges = tables["fixed_edges"]
    directions = tables["directions"]
    direction_signs = tables["direction_signs"]
    budgets = np.asarray(tables["budgets"], dtype=np.int16)
    assert isinstance(halves, tuple)
    assert isinstance(orbit_cache, tuple)
    assert isinstance(fixed_edges, tuple)
    assert isinstance(directions, tuple)
    assert isinstance(direction_signs, tuple)
    if len(centers) != len(halves) or any(not 1 <= int(value) < P for value in centers):
        raise ValueError("invalid exact-replay center tuple")

    total: Counter[Edge] = Counter()
    for half_index, center in enumerate(centers):
        total.update(orbit_cache[half_index][int(center) - 1])
    vanished = tuple(sorted(orbit for orbit, value in total.items() if value == 0))
    surviving = {orbit: value for orbit, value in total.items() if value}
    if len(vanished) != 1 or len(surviving) != 478 or any(
        abs(value) != 1 for value in surviving.values()
    ):
        raise ArithmeticError("candidate lost one-cancellation ternary support")
    graph = [
        orbit if value == 1 else _negative_edge(P, orbit)
        for orbit, value in surviving.items()
    ]
    fixed_edge = fixed_edges[fixed_index]
    graph.append(fixed_edge)
    graph = sorted(graph)
    if len(graph) != 479 or len(set(graph)) != 479:
        raise ArithmeticError("candidate is not a 479-edge simple graph")

    source = {edge: paley_edge_sign(P, edge) for edge in graph}
    image = edge_radon_image(P, source)
    literal_centers = {
        _projective_index(target): _canonical_center(target, int(center))[1]
        for (target, _auxiliary), center in zip(halves, centers, strict=True)
    }
    parallel = []
    row_l1 = []
    for direction_index, sign in enumerate(direction_signs):
        parallel.append(sign * image.get(("P", direction_index), 0))
        literal = literal_centers.get(direction_index) if sign == 1 else None
        coefficients = []
        for left, right in combinations(range(P), 2):
            value = sign * image.get(("K", direction_index, left, right), 0)
            if literal is not None and literal in (left, right):
                value += 1
            coefficients.append(value)
        row_l1.append(sum(abs(value) for value in coefficients))
    if tuple(parallel) != tables["final_parallel_profile"]:
        raise ArithmeticError("exact graph replay missed the final parallel profile")
    excess = [
        max(0, value - int(budget))
        for value, budget in zip(row_l1, budgets, strict=True)
    ]
    graph_bytes = json.dumps(graph, separators=(",", ":")).encode()
    return {
        "centers": [int(value) for value in centers],
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


def _selected(design_index: int, parity: str) -> bool:
    return parity == "all" or design_index % 2 == (1 if parity == "odd" else 0)


def run(args: argparse.Namespace) -> dict[str, object]:
    design_list_sha256 = hashlib.sha256(args.design_list.read_bytes()).hexdigest()
    requested_indices = (
        None if args.design_index is None else set(args.design_index)
    )
    designs = tuple(
        design
        for design in load_designs(args.design_list)
        if int(design.get("clean_collision_seed_count", 0)) == 30
        and _selected(int(design["design_index"]), args.design_parity)
        and (
            requested_indices is None
            or int(design["design_index"]) in requested_indices
        )
    )
    if not designs:
        raise ValueError("the requested design shard is empty")
    result: dict[str, object] = {
        "schema": "residual_branch_c_design_transverse_opencl_v1",
        "classification": "bounded exact-objective search over supplied moved-design quotient fibres",
        "host": socket.gethostname(),
        "design_list": str(args.design_list),
        "design_list_sha256": design_list_sha256,
        "design_parity": args.design_parity,
        "selected_design_indices": [int(row["design_index"]) for row in designs],
        "seed": args.seed,
        "runs": [],
        "residual_ii_closed": False,
        "design_fibres_exhausted": False,
    }
    started = time.perf_counter()
    runs = result["runs"]
    assert isinstance(runs, list)
    for ordinal, design in enumerate(designs):
        design_index = int(design["design_index"])
        seed = _canonical_seed(design)
        tables = build_design_tables(design, seed)
        half_indices = tuple(int(value) for value in seed["half_indices"])
        collision_centers = tuple(int(value) for value in seed["centers"])
        search = opencl_coordinate_search(
            tables,
            chains=args.chains,
            restarts=args.restarts,
            sweeps=args.sweeps,
            seed=args.seed + 1_000_003 * design_index,
            primary=args.primary,
            time_limit=args.seconds_per_design,
            frozen_pair=half_indices,
            frozen_centers=collision_centers,
        )
        best = search["best"]
        replay = exact_design_graph_replay(
            best["centers"], int(best["fixed_edge_index"]), tables
        )
        for key in ("row_l1", "total_positive_l1_excess", "maximum_l1_excess"):
            if replay[key] != best[key]:
                raise ArithmeticError(
                    f"design {design_index} failed exact CPU replay at {key}"
                )
        record = {
            "design_index": design_index,
            "move_id": design.get("move_id"),
            "design_sha256": design.get("design_sha256"),
            "canonical_collision_seed": seed,
            "search": search,
            "best_exact_replay": replay,
        }
        runs.append(record)
        result["elapsed_seconds"] = time.perf_counter() - started
        result["designs_completed"] = len(runs)
        result["best"] = min(
            runs,
            key=lambda row: (
                row["best_exact_replay"]["total_positive_l1_excess"],
                row["best_exact_replay"]["maximum_l1_excess"],
            ),
        )
        if args.output:
            write_json_atomic(args.output, result)
        print(
            json.dumps(
                {
                    "design_complete": design_index,
                    "ordinal": ordinal,
                    "best": replay,
                }
            ),
            flush=True,
        )
        if replay["passes_all_l1_rows"]:
            break
    result["elapsed_seconds"] = time.perf_counter() - started
    result["designs_completed"] = len(runs)
    result["all_selected_designs_completed"] = len(runs) == len(designs)
    if args.output:
        write_json_atomic(args.output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--design-list", type=Path, required=True)
    parser.add_argument("--design-parity", choices=("odd", "even", "all"), default="odd")
    parser.add_argument("--design-index", type=int, action="append")
    parser.add_argument("--chains", type=int, default=512)
    parser.add_argument("--restarts", type=int, default=100_000)
    parser.add_argument("--sweeps", type=int, default=8)
    parser.add_argument("--seconds-per-design", type=float, default=10.0)
    parser.add_argument("--seed", type=int, default=920_260_904)
    parser.add_argument("--primary", choices=("total", "max"), default="total")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if (
        args.chains < 1
        or args.restarts < 1
        or args.sweeps < 1
        or args.seconds_per_design <= 0
    ):
        parser.error("all search bounds must be positive")
    print(json.dumps(run(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
