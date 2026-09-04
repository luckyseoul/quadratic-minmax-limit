#!/usr/bin/env python3
"""Bounded V100 search across exact p=31 auxiliary-profile fibres.

This is a construction search, not an exhaustive certificate.  It starts
from the published sixteen-half top parallel design and enumerates every
nontrivial two-target replacement which preserves its *entire* raw
32-direction parallel vector.  Thus auxiliary projective directions and
relative scales genuinely change while the corrected top profile remains
exact.  Only moved designs admitting a clean opposite-oriented one-orbit
cancellation in the required fixed direction are sent to the GPU.

There are thirty scalar copies of every clean collision seed.  Simultaneous
dilation of all centers and the fixed edge merely permutes affine labels in
every transverse row, so their row objectives are identical.  The search
stores all physical seeds but optimizes one canonical representative of each
scalar orbit.  Each reported winner is rebuilt from its actual 479 graph
edges and replayed with the integer Radon map.

Passing the compact-atom l1 budgets is necessary, not sufficient, for the
remaining symmetric Boolean fibre.  A bounded miss proves no obstruction.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import hashlib
from itertools import combinations, product
import json
import os
from pathlib import Path
import random
import socket
import sys
import time
from typing import Iterable, Sequence

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
    _parallel_formula,
    paley_direction_sign,
    paley_edge_sign,
)
import e1_gmin_m4_p31_direct_mobius_parallel_design as frozen  # noqa: E402
from io_atomic import write_json_atomic  # noqa: E402


P = 31
DIRECTIONS = tuple(projective_functionals(P))
DIRECTION_SIGNS = tuple(paley_direction_sign(P, row) for row in DIRECTIONS)
HARD = tuple(index for index, sign in enumerate(DIRECTION_SIGNS) if sign == 1)
OPPOSITE = tuple(index for index, sign in enumerate(DIRECTION_SIGNS) if sign == -1)
FIXED_DIRECTION = frozen.PHYSICAL_FIXED_DIRECTION_INDEX
CANCELLED_DIRECTION = frozen.PHYSICAL_CANCELLATION_DIRECTION_INDEX


@dataclass(frozen=True)
class HalfChoice:
    target_index: int
    auxiliary_index: int
    relative_scale: int

    @property
    def target(self) -> Functional:
        return DIRECTIONS[self.target_index]

    @property
    def auxiliary(self) -> Functional:
        row = DIRECTIONS[self.auxiliary_index]
        return (
            self.relative_scale * row[0] % P,
            self.relative_scale * row[1] % P,
        )


@dataclass(frozen=True)
class CollisionSeed:
    first_half: int
    second_half: int
    first_center: int
    second_center: int
    orbit: Edge
    first_coefficient: int
    second_coefficient: int


def _functional_scale(functional: Functional, canonical: Functional) -> int:
    if canonical[0]:
        return functional[0] * pow(canonical[0], -1, P) % P
    return functional[1] * pow(canonical[1], -1, P) % P


def _choice_profile(choice: HalfChoice) -> tuple[int, ...]:
    return tuple(
        _parallel_formula(P, choice.target, choice.auxiliary, row)
        for row in DIRECTIONS
    )


def canonical_frozen_design() -> tuple[HalfChoice, ...]:
    """Normalize the published halves to canonical target representatives."""
    out: list[HalfChoice | None] = [None] * len(HARD)
    for target, auxiliary in frozen.HALVES:
        target_index = frozen._projective_index(target)
        position = HARD.index(target_index)
        target_scale = _functional_scale(target, DIRECTIONS[target_index])
        auxiliary_index = frozen._projective_index(auxiliary)
        auxiliary_scale = _functional_scale(
            auxiliary, DIRECTIONS[auxiliary_index]
        )
        relative_scale = auxiliary_scale * pow(target_scale, -1, P) % P
        out[position] = HalfChoice(
            target_index, auxiliary_index, relative_scale
        )
    if any(row is None for row in out):
        raise ArithmeticError("the frozen design lost a hard target")
    result = tuple(row for row in out if row is not None)
    if len({row.auxiliary_index for row in result}) != len(HARD):
        raise ArithmeticError("the frozen auxiliary SDR changed")
    return result


def canonical_frozen_centers() -> tuple[int, ...]:
    out: list[int | None] = [None] * len(HARD)
    for (target, _auxiliary), center in zip(
        frozen.HALVES, frozen.PHYSICAL_CENTERS, strict=True
    ):
        target_index, canonical_center = frozen._canonical_center(target, center)
        out[HARD.index(target_index)] = canonical_center
    if any(center is None for center in out):
        raise ArithmeticError("the frozen center normalization changed")
    return tuple(int(center) for center in out)


def _all_options() -> tuple[tuple[tuple[HalfChoice, tuple[int, ...]], ...], ...]:
    rows = []
    for target_index in HARD:
        options = []
        for auxiliary_index in range(P + 1):
            if auxiliary_index == target_index:
                continue
            for relative_scale in range(1, P):
                choice = HalfChoice(
                    target_index, auxiliary_index, relative_scale
                )
                options.append((choice, _choice_profile(choice)))
        rows.append(tuple(options))
    return tuple(rows)


def raw_parallel_profile(design: Sequence[HalfChoice]) -> tuple[int, ...]:
    if len(design) != len(HARD):
        raise ValueError("a design needs sixteen halves")
    return tuple(
        sum(_choice_profile(choice)[direction] for choice in design)
        for direction in range(P + 1)
    )


def exact_pair_profile_moves() -> tuple[tuple[HalfChoice, ...], ...]:
    """Enumerate every two-target move preserving the frozen raw profile."""
    base = canonical_frozen_design()
    options = _all_options()
    used = {choice.auxiliary_index for choice in base}
    base_key = tuple(
        (choice.auxiliary_index, choice.relative_scale) for choice in base
    )
    found: dict[tuple[tuple[int, int], ...], tuple[HalfChoice, ...]] = {}
    base_profiles = tuple(_choice_profile(choice) for choice in base)
    for first, second in combinations(range(len(HARD)), 2):
        old_sum = tuple(
            a + b
            for a, b in zip(
                base_profiles[first], base_profiles[second], strict=True
            )
        )
        second_by_profile: defaultdict[
            tuple[int, ...], list[HalfChoice]
        ] = defaultdict(list)
        for choice, profile in options[second]:
            second_by_profile[profile].append(choice)
        unavailable = used - {
            base[first].auxiliary_index,
            base[second].auxiliary_index,
        }
        for first_choice, first_profile in options[first]:
            if first_choice.auxiliary_index in unavailable:
                continue
            needed = tuple(
                total - value
                for total, value in zip(old_sum, first_profile, strict=True)
            )
            for second_choice in second_by_profile.get(needed, ()):
                if (
                    second_choice.auxiliary_index in unavailable
                    or second_choice.auxiliary_index
                    == first_choice.auxiliary_index
                ):
                    continue
                moved = list(base)
                moved[first] = first_choice
                moved[second] = second_choice
                key = tuple(
                    (choice.auxiliary_index, choice.relative_scale)
                    for choice in moved
                )
                if key != base_key:
                    found[key] = tuple(moved)
    result = tuple(found[key] for key in sorted(found))
    frozen_raw = raw_parallel_profile(base)
    if (
        len(result) != 32
        or any(raw_parallel_profile(row) != frozen_raw for row in result)
        or any(
            len({choice.auxiliary_index for choice in row}) != len(HARD)
            for row in result
        )
        or any(
            sum(a != b for a, b in zip(row, base, strict=True)) != 2
            for row in result
        )
    ):
        raise ArithmeticError("the exact two-target move catalog changed")
    return result


def design_sha256(design: Sequence[HalfChoice]) -> str:
    payload = [
        [choice.target_index, choice.auxiliary_index, choice.relative_scale]
        for choice in design
    ]
    return hashlib.sha256(
        json.dumps(payload, separators=(",", ":")).encode()
    ).hexdigest()


def _orbit_cache(
    design: Sequence[HalfChoice],
) -> tuple[tuple[dict[Edge, int], ...], ...]:
    return tuple(
        tuple(
            frozen._oriented_orbit_coefficients(
                choice.target, choice.auxiliary, center
            )
            for center in range(1, P)
        )
        for choice in design
    )


def clean_collision_seeds(
    design: Sequence[HalfChoice],
    cache: tuple[tuple[dict[Edge, int], ...], ...] | None = None,
) -> tuple[CollisionSeed, ...]:
    """Return all clean pair collisions in the required spatial direction."""
    if cache is None:
        cache = _orbit_cache(design)
    seeds = []
    for first, second in combinations(range(len(HARD)), 2):
        for first_center in range(1, P):
            first_map = cache[first][first_center - 1]
            for second_center in range(1, P):
                second_map = cache[second][second_center - 1]
                common = set(first_map) & set(second_map)
                if len(common) != 1:
                    continue
                orbit = next(iter(common))
                if (
                    first_map[orbit] == -second_map[orbit]
                    and frozen._spatial_direction_index(orbit)
                    == CANCELLED_DIRECTION
                ):
                    seeds.append(
                        CollisionSeed(
                            first,
                            second,
                            first_center,
                            second_center,
                            orbit,
                            first_map[orbit],
                            second_map[orbit],
                        )
                    )
    return tuple(seeds)


def _cell_catalog() -> tuple[
    tuple[tuple[int, int], ...], dict[tuple[int, int], int]
]:
    cells = tuple(combinations(range(P), 2))
    return cells, {cell: index for index, cell in enumerate(cells)}


def _accumulate_edge_projection(
    out: np.ndarray,
    edge: Edge,
    cell_index: dict[tuple[int, int], int],
) -> None:
    """Add one sparse edge projection without allocating a dense tensor."""
    tau = paley_edge_sign(P, edge)
    for direction_index, direction in enumerate(DIRECTIONS):
        first = _functional_value(P, direction, edge[0])
        second = _functional_value(P, direction, edge[1])
        if first != second:
            out[
                direction_index,
                cell_index[tuple(sorted((first, second)))],
            ] += DIRECTION_SIGNS[direction_index] * tau


def _edge_projection(
    edge: Edge,
    cell_index: dict[tuple[int, int], int],
) -> np.ndarray:
    out = np.zeros((P + 1, len(cell_index)), dtype=np.int16)
    _accumulate_edge_projection(out, edge, cell_index)
    return out


def _fixed_edges() -> tuple[Edge, ...]:
    direction = DIRECTIONS[FIXED_DIRECTION]
    fixed = set()
    for point in product(range(P), repeat=2):
        if point == (0, 0) or _functional_value(P, direction, point):
            continue
        negative = (-point[0] % P, -point[1] % P)
        fixed.add(tuple(sorted((point, negative))))
    result = tuple(sorted(fixed))
    if len(result) != (P - 1) // 2:
        raise ArithmeticError("the fixed edge fibre changed")
    return result


def build_design_tables(
    design: Sequence[HalfChoice],
    seed: CollisionSeed,
    cache: tuple[tuple[dict[Edge, int], ...], ...] | None = None,
) -> dict[str, object]:
    if cache is None:
        cache = _orbit_cache(design)
    cells, cell_index = _cell_catalog()
    half_tables = np.zeros(
        (len(HARD), P - 1, P + 1, len(cells)), dtype=np.int16
    )
    for half_index, (choice, center_rows) in enumerate(
        zip(design, cache, strict=True)
    ):
        for center_index, orbit_map in enumerate(center_rows):
            row = half_tables[half_index, center_index]
            for orbit, value in orbit_map.items():
                edge = orbit if value == 1 else _negative_edge(P, orbit)
                _accumulate_edge_projection(row, edge, cell_index)
            center = center_index + 1
            for other in range(P):
                if other != center:
                    row[
                        choice.target_index,
                        cell_index[tuple(sorted((center, other)))],
                    ] += 1

    first_map = cache[seed.first_half][seed.first_center - 1]
    second_map = cache[seed.second_half][seed.second_center - 1]
    collision_edges = (
        seed.orbit
        if first_map[seed.orbit] == 1
        else _negative_edge(P, seed.orbit),
        seed.orbit
        if second_map[seed.orbit] == 1
        else _negative_edge(P, seed.orbit),
    )
    collision_table = np.zeros((P + 1, len(cells)), dtype=np.int16)
    for edge in collision_edges:
        _accumulate_edge_projection(collision_table, edge, cell_index)
    fixed_edges = _fixed_edges()
    fixed_tables = np.zeros(
        (len(fixed_edges), P + 1, len(cells)), dtype=np.int16
    )
    for index, edge in enumerate(fixed_edges):
        _accumulate_edge_projection(fixed_tables[index], edge, cell_index)

    pair_bad = np.zeros(
        (len(HARD), len(HARD), P - 1, P - 1), dtype=np.bool_
    )
    supports = tuple(
        tuple(set(row) for row in center_rows) for center_rows in cache
    )
    for first, second in combinations(range(len(HARD)), 2):
        for first_center in range(P - 1):
            for second_center in range(P - 1):
                bad = bool(
                    supports[first][first_center]
                    & supports[second][second_center]
                )
                pair_bad[
                    first, second, first_center, second_center
                ] = bad
                pair_bad[
                    second, first, second_center, first_center
                ] = bad

    final_profile = frozen.centered_physical_parallel_design_certificate()[
        "final_parallel_profile"
    ]
    budgets = np.asarray(
        [3 * (int(value) - 3) for value in final_profile], dtype=np.int16
    )
    return {
        "orbit_cache": cache,
        "half_tables": half_tables,
        "collision_table": collision_table,
        "fixed_edges": fixed_edges,
        "fixed_tables": fixed_tables,
        "pair_bad": pair_bad,
        "budgets": budgets,
        "cells": cells,
        "cell_index": cell_index,
    }


def _valid_assignment(
    choices: np.ndarray,
    pair_bad: np.ndarray,
    collision_pair: tuple[int, int],
) -> bool:
    for first, second in combinations(range(len(HARD)), 2):
        if (first, second) == collision_pair:
            continue
        if pair_bad[first, second, choices[first], choices[second]]:
            return False
    return True


def _random_assignments(
    count: int,
    pair_bad: np.ndarray,
    collision: CollisionSeed,
    rng: random.Random,
) -> np.ndarray:
    pair = (collision.first_half, collision.second_half)
    out = []
    attempts = 0
    while len(out) < count and attempts < 1000 * count:
        attempts += 1
        choices = np.full(len(HARD), -1, dtype=np.int16)
        choices[collision.first_half] = collision.first_center - 1
        choices[collision.second_half] = collision.second_center - 1
        order = [index for index in range(len(HARD)) if index not in pair]
        rng.shuffle(order)
        valid = True
        for index in order:
            options = list(range(P - 1))
            rng.shuffle(options)
            selected = next(
                (
                    option
                    for option in options
                    if all(
                        choices[other] < 0
                        or not pair_bad[index, other, option, choices[other]]
                        for other in range(len(HARD))
                        if other != index
                    )
                ),
                None,
            )
            if selected is None:
                valid = False
                break
            choices[index] = selected
        if valid and _valid_assignment(choices, pair_bad, pair):
            out.append(choices)
    if len(out) != count:
        raise RuntimeError(
            f"only generated {len(out)} of {count} disjoint assignments"
        )
    return np.stack(out)


def _candidate_valid_mask(
    choices: np.ndarray,
    variable: int,
    option: int,
    pair_bad: np.ndarray,
) -> np.ndarray:
    valid = np.ones(len(choices), dtype=np.bool_)
    for other in range(len(HARD)):
        if other != variable:
            valid &= ~pair_bad[variable, other, option, choices[:, other]]
    return valid


def gpu_center_search(
    tables: dict[str, object],
    collision: CollisionSeed,
    *,
    chains: int,
    restarts: int,
    sweeps: int,
    seed: int,
    primary: str,
) -> dict[str, object]:
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
    pair = (collision.first_half, collision.second_half)
    mutable = [index for index in range(len(HARD)) if index not in pair]
    evaluated = 0
    best: dict[str, object] | None = None
    started = time.perf_counter()

    for restart in range(restarts):
        choices = _random_assignments(chains, pair_bad, collision, rng)
        fixed_choice = np.asarray(
            [rng.randrange(len(fixed_tables)) for _ in range(chains)],
            dtype=np.int16,
        )
        state = fixed_gpu[cp.asarray(fixed_choice)] - collision_gpu[None, :, :]
        for half_index in range(len(HARD)):
            state += half_gpu[
                half_index, cp.asarray(choices[:, half_index])
            ]

        for _sweep in range(sweeps):
            changed = False
            rng.shuffle(mutable)
            for variable in mutable:
                old = half_gpu[variable, cp.asarray(choices[:, variable])]
                base = state - old
                best_total = np.full(
                    chains, np.iinfo(np.int64).max, dtype=np.int64
                )
                best_max = np.full(
                    chains, np.iinfo(np.int64).max, dtype=np.int64
                )
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
                    totals = cp.asnumpy(excess.sum(axis=1))
                    maxima = cp.asnumpy(excess.max(axis=1))
                    totals[~valid] = np.iinfo(np.int64).max
                    maxima[~valid] = np.iinfo(np.int64).max
                    if primary == "total":
                        improve = (totals < best_total) | (
                            (totals == best_total) & (maxima < best_max)
                        )
                    else:
                        improve = (maxima < best_max) | (
                            (maxima == best_max) & (totals < best_total)
                        )
                    best_total[improve] = totals[improve]
                    best_max[improve] = maxima[improve]
                    best_option[improve] = option
                    evaluated += int(valid.sum())
                changed |= bool(np.any(best_option != choices[:, variable]))
                choices[:, variable] = best_option
                state = base + half_gpu[variable, cp.asarray(best_option)]

            base = state - fixed_gpu[cp.asarray(fixed_choice)]
            best_total = np.full(
                chains, np.iinfo(np.int64).max, dtype=np.int64
            )
            best_max = np.full(
                chains, np.iinfo(np.int64).max, dtype=np.int64
            )
            best_fixed = fixed_choice.copy()
            for option in range(len(fixed_tables)):
                candidate = base + fixed_gpu[option]
                l1 = cp.abs(candidate).sum(axis=2)
                excess = cp.maximum(l1 - budgets_gpu, 0)
                totals = cp.asnumpy(excess.sum(axis=1))
                maxima = cp.asnumpy(excess.max(axis=1))
                if primary == "total":
                    improve = (totals < best_total) | (
                        (totals == best_total) & (maxima < best_max)
                    )
                else:
                    improve = (maxima < best_max) | (
                        (maxima == best_max) & (totals < best_total)
                    )
                best_total[improve] = totals[improve]
                best_max[improve] = maxima[improve]
                best_fixed[improve] = option
                evaluated += chains
            changed |= bool(np.any(best_fixed != fixed_choice))
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
        key = (
            (record["total_positive_l1_excess"], record["maximum_l1_excess"])
            if primary == "total"
            else (record["maximum_l1_excess"], record["total_positive_l1_excess"])
        )
        best_key = None if best is None else (
            (best["total_positive_l1_excess"], best["maximum_l1_excess"])
            if primary == "total"
            else (best["maximum_l1_excess"], best["total_positive_l1_excess"])
        )
        if best_key is None or key < best_key:
            best = record
        if record["total_positive_l1_excess"] == 0:
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
    design: Sequence[HalfChoice],
    collision: CollisionSeed,
    centers: Sequence[int],
    fixed_index: int,
    tables: dict[str, object],
) -> dict[str, object]:
    cache = tables["orbit_cache"]
    fixed_edges = tables["fixed_edges"]
    budgets = tables["budgets"]
    pair_bad = tables["pair_bad"]
    assert isinstance(cache, tuple)
    assert isinstance(fixed_edges, tuple)
    assert isinstance(budgets, np.ndarray)
    assert isinstance(pair_bad, np.ndarray)
    choices = np.asarray([int(center) - 1 for center in centers], dtype=np.int16)
    pair = (collision.first_half, collision.second_half)
    if not _valid_assignment(choices, pair_bad, pair):
        raise ArithmeticError("GPU winner acquired an extra overlap")

    total: Counter[Edge] = Counter()
    occurrences: defaultdict[Edge, list[tuple[int, int]]] = defaultdict(list)
    for half_index, center in enumerate(centers):
        orbit_map = cache[half_index][int(center) - 1]
        total.update(orbit_map)
        for orbit, coefficient in orbit_map.items():
            occurrences[orbit].append((half_index, coefficient))
    vanished = tuple(sorted(orbit for orbit, value in total.items() if value == 0))
    surviving = {orbit: value for orbit, value in total.items() if value}
    if (
        vanished != (collision.orbit,)
        or len(surviving) != 478
        or any(abs(value) != 1 for value in surviving.values())
    ):
        raise ArithmeticError("GPU winner lost clean one-orbit ternarity")
    graph = [
        orbit if value == 1 else _negative_edge(P, orbit)
        for orbit, value in surviving.items()
    ]
    fixed_edge = fixed_edges[fixed_index]
    graph.append(fixed_edge)
    graph = tuple(sorted(graph))
    if len(graph) != 479 or len(set(graph)) != 479:
        raise ArithmeticError("GPU winner is not a 479-edge simple graph")

    image = edge_radon_image(
        P, {edge: paley_edge_sign(P, edge) for edge in graph}
    )
    parallel = tuple(
        DIRECTION_SIGNS[index] * image.get(("P", index), 0)
        for index in range(P + 1)
    )
    expected_parallel = tuple(
        frozen.centered_physical_parallel_design_certificate()[
            "final_parallel_profile"
        ]
    )
    if parallel != expected_parallel:
        raise ArithmeticError("GPU winner failed its exact parallel replay")

    rows = []
    for direction_index, sign in enumerate(DIRECTION_SIGNS):
        coefficients: dict[tuple[int, int], int] = {}
        for left, right in combinations(range(P), 2):
            value = sign * image.get(
                ("K", direction_index, left, right), 0
            )
            if sign == 1:
                center = centers[HARD.index(direction_index)]
                if center in (left, right):
                    value += 1
            if value:
                coefficients[(left, right)] = value
        l1 = sum(abs(value) for value in coefficients.values())
        budget = int(budgets[direction_index])
        rows.append(
            {
                "direction_index": direction_index,
                "l1": l1,
                "budget": budget,
                "positive_excess": max(0, l1 - budget),
                "nonzero_cell_count": len(coefficients),
            }
        )
    graph_bytes = json.dumps(graph, separators=(",", ":")).encode()
    return {
        "centers": tuple(int(center) for center in centers),
        "fixed_edge": fixed_edge,
        "cancelled_orbit": collision.orbit,
        "collision_half_indices": pair,
        "graph_edge_count": len(graph),
        "graph_sha256": hashlib.sha256(graph_bytes).hexdigest(),
        "parallel_profile": parallel,
        "row_l1": tuple(row["l1"] for row in rows),
        "total_positive_l1_excess": sum(
            row["positive_excess"] for row in rows
        ),
        "maximum_l1_excess": max(row["positive_excess"] for row in rows),
        "violating_row_count": sum(row["positive_excess"] > 0 for row in rows),
        "passes_all_l1_rows": all(row["positive_excess"] == 0 for row in rows),
        "rows": tuple(rows),
        "exact_integer_graph_replay": True,
    }


def _collision_record(seed: CollisionSeed) -> dict[str, object]:
    return {
        "half_indices": (seed.first_half, seed.second_half),
        "centers": (seed.first_center, seed.second_center),
        "orbit": seed.orbit,
        "coefficients": (seed.first_coefficient, seed.second_coefficient),
        "spatial_direction_index": CANCELLED_DIRECTION,
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    base = canonical_frozen_design()
    moved = exact_pair_profile_moves()
    catalog = []
    searchable = []
    for design_index, design in enumerate(moved):
        cache = _orbit_cache(design)
        collisions = clean_collision_seeds(design, cache)
        record = {
            "design_index": design_index,
            "move_id": design_index,
            "design_sha256": design_sha256(design),
            "halves": tuple(
                (choice.target, choice.auxiliary) for choice in design
            ),
            "choices": tuple(
                {
                    "target_direction_index": row.target_index,
                    "auxiliary_direction_index": row.auxiliary_index,
                    "relative_scale": row.relative_scale,
                }
                for row in design
            ),
            "raw_parallel_profile": raw_parallel_profile(design),
            "clean_collision_seed_count": len(collisions),
            "clean_collision_seeds": tuple(
                _collision_record(seed) for seed in collisions
            ),
            "fixed_direction_index": FIXED_DIRECTION,
        }
        catalog.append(record)
        if collisions:
            # The thirty listed seeds form one simultaneous-dilation orbit.
            representative = min(
                collisions,
                key=lambda row: (
                    row.first_half,
                    row.second_half,
                    row.first_center,
                    row.second_center,
                    row.orbit,
                ),
            )
            searchable.append((design_index, design, cache, representative))

    if len(moved) != 32 or sum(row["clean_collision_seed_count"] for row in catalog) != 630:
        raise ArithmeticError("the moved-design collision catalog changed")
    catalog_payload = {
        "schema": "residual_branch_c_auxiliary_pair_move_catalog_v1",
        "p": P,
        "hard_target_order": HARD,
        "fixed_direction_index": FIXED_DIRECTION,
        "cancelled_direction_index": CANCELLED_DIRECTION,
        "base_raw_parallel_profile": raw_parallel_profile(base),
        "designs": tuple(catalog),
    }
    if args.catalog_output:
        write_json_atomic(args.catalog_output, catalog_payload)
    if args.catalog_only:
        return catalog_payload

    selected = searchable
    if args.design_parity != "all":
        residue = 0 if args.design_parity == "even" else 1
        selected = [row for row in selected if row[0] % 2 == residue]
    selected = selected[args.design_start :]
    if args.max_designs is not None:
        selected = selected[: args.max_designs]
    searches = []
    global_best: dict[str, object] | None = None
    total_evaluated = 0
    total_gpu_seconds = 0.0
    for ordinal, (design_index, design, cache, collision) in enumerate(selected):
        tables = build_design_tables(design, collision, cache)
        search = gpu_center_search(
            tables,
            collision,
            chains=args.chains,
            restarts=args.restarts,
            sweeps=args.sweeps,
            seed=args.seed + 104729 * design_index,
            primary=args.primary,
        )
        winner = search["best"]
        assert isinstance(winner, dict)
        replay = exact_graph_replay(
            design,
            collision,
            winner["centers"],
            int(winner["fixed_edge_index"]),
            tables,
        )
        if (
            replay["row_l1"] != tuple(winner["row_l1"])
            or replay["total_positive_l1_excess"]
            != winner["total_positive_l1_excess"]
            or replay["maximum_l1_excess"] != winner["maximum_l1_excess"]
        ):
            raise ArithmeticError("V100 winner failed exact CPU replay")
        row = {
            "design_index": design_index,
            "design_sha256": design_sha256(design),
            "scalar_collision_seed_count": len(clean_collision_seeds(design, cache)),
            "canonical_collision_seed": _collision_record(collision),
            "gpu_search": search,
            "best_exact_replay": replay,
        }
        searches.append(row)
        total_evaluated += int(search["objective_states_evaluated"])
        total_gpu_seconds += float(search["elapsed_seconds"])
        key = (
            (replay["total_positive_l1_excess"], replay["maximum_l1_excess"])
            if args.primary == "total"
            else (replay["maximum_l1_excess"], replay["total_positive_l1_excess"])
        )
        old_key = None if global_best is None else (
            (
                global_best["best_exact_replay"]["total_positive_l1_excess"],
                global_best["best_exact_replay"]["maximum_l1_excess"],
            )
            if args.primary == "total"
            else (
                global_best["best_exact_replay"]["maximum_l1_excess"],
                global_best["best_exact_replay"]["total_positive_l1_excess"],
            )
        )
        if old_key is None or key < old_key:
            global_best = row
            print(
                json.dumps(
                    {
                        "new_global_best": {
                            "design_index": design_index,
                            "design_sha256": row["design_sha256"],
                            "total_positive_l1_excess": replay[
                                "total_positive_l1_excess"
                            ],
                            "maximum_l1_excess": replay[
                                "maximum_l1_excess"
                            ],
                            "graph_sha256": replay["graph_sha256"],
                        },
                        "search_ordinal": ordinal,
                    }
                ),
                flush=True,
            )
        if replay["passes_all_l1_rows"]:
            break

    result = {
        "schema": "residual_branch_c_auxiliary_transverse_gpu_v1",
        "classification": "bounded exact-objective construction search",
        "host": socket.gethostname(),
        "p": P,
        "base_design_sha256": design_sha256(base),
        "base_raw_parallel_profile": raw_parallel_profile(base),
        "exact_pair_moved_design_count": len(moved),
        "moved_designs_with_clean_collision": len(searchable),
        "physical_clean_collision_seed_count": sum(
            row["clean_collision_seed_count"] for row in catalog
        ),
        "scalar_orbit_reduction": (
            "30 physical seeds per searchable design reduce to one objective "
            "fibre by simultaneous nonzero dilation"
        ),
        "design_catalog": tuple(catalog),
        "searched_design_start": args.design_start,
        "searched_design_parity": args.design_parity,
        "searched_design_count": len(searches),
        "searches": tuple(searches),
        "aggregate_objective_states_evaluated": total_evaluated,
        "aggregate_gpu_seconds": total_gpu_seconds,
        "aggregate_objective_states_per_second": (
            total_evaluated / total_gpu_seconds if total_gpu_seconds else 0.0
        ),
        "global_best": global_best,
        "l1_pass_is_sufficient_for_atom_decomposition": False,
        "pair_move_catalog_exhausted": True,
        "center_fibres_exhausted": False,
        "residual_ii_closed": False,
    }
    if args.output:
        write_json_atomic(args.output, result)
    return result


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chains", type=int, default=64)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--sweeps", type=int, default=6)
    parser.add_argument("--seed", type=int, default=2026090401)
    parser.add_argument("--primary", choices=("total", "max"), default="total")
    parser.add_argument("--design-start", type=int, default=0)
    parser.add_argument("--design-parity", choices=("all", "even", "odd"), default="all")
    parser.add_argument("--max-designs", type=int)
    parser.add_argument("--catalog-output", type=Path)
    parser.add_argument("--catalog-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    if (
        args.chains < 1
        or args.restarts < 1
        or args.sweeps < 1
        or args.design_start < 0
        or (args.max_designs is not None and args.max_designs < 1)
    ):
        parser.error("search counts must be positive and design-start nonnegative")
    return args


def main(argv: Iterable[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    result = run(args)
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


if __name__ == "__main__":
    main()
