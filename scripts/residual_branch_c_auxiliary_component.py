#!/usr/bin/env python3
"""Bounded/resumable exact two-half component enumerator for the p31 design.

Vertices are canonical sixteen-tuples ``(auxiliary direction, relative
scale)`` in the fixed hard-target order.  An edge changes exactly two halves
and preserves their full 32-entry parallel-profile sum.  The parity-forced
auxiliary direction set is invariant, so a neighbor can only retain or swap
the two released auxiliary directions; this is used as an exact acceleration,
not as a heuristic filter.

The component is large.  ``--max-designs`` is therefore a hard discovery
bound, and the JSON stores every canonical key plus the remaining BFS queue
so a later run can resume without revisiting work.  Search candidates retain
the published clean collision pair (halves 3 and 4 in canonical target order),
whose thirty physical center seeds are factored into one deduplicated catalog.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sys
import time
from typing import Iterable, Iterator, Sequence


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from io_atomic import write_json_atomic  # noqa: E402
import residual_branch_c_auxiliary_transverse_gpu as auxiliary  # noqa: E402


Key = tuple[tuple[int, int], ...]
COLLISION_PAIR = (3, 4)


def design_key(design: Sequence[auxiliary.HalfChoice]) -> Key:
    return tuple(
        (choice.auxiliary_index, choice.relative_scale) for choice in design
    )


def design_from_key(key: Key) -> tuple[auxiliary.HalfChoice, ...]:
    if len(key) != len(auxiliary.HARD):
        raise ValueError("a component key needs sixteen entries")
    design = tuple(
        auxiliary.HalfChoice(target, int(row[0]), int(row[1]))
        for target, row in zip(auxiliary.HARD, key, strict=True)
    )
    if len({choice.auxiliary_index for choice in design}) != len(design):
        raise ValueError("a component key repeated an auxiliary direction")
    return design


class PairMoveIndex:
    """Exact profile lookup for fast two-half neighbor generation."""

    def __init__(self) -> None:
        self.profile: list[dict[tuple[int, int], tuple[int, ...]]] = [
            {} for _ in auxiliary.HARD
        ]
        self.scales_by_profile: list[
            dict[int, dict[tuple[int, ...], tuple[int, ...]]]
        ] = [{} for _ in auxiliary.HARD]
        for position, rows in enumerate(auxiliary._all_options()):
            for choice, profile in rows:
                self.profile[position][
                    (choice.auxiliary_index, choice.relative_scale)
                ] = profile
            for direction in range(auxiliary.P + 1):
                if direction == auxiliary.HARD[position]:
                    continue
                grouped: defaultdict[tuple[int, ...], list[int]] = defaultdict(list)
                for scale in range(1, auxiliary.P):
                    grouped[self.profile[position][(direction, scale)]].append(scale)
                self.scales_by_profile[position][direction] = {
                    profile: tuple(scales) for profile, scales in grouped.items()
                }

    def neighbors(self, key: Key) -> Iterator[tuple[Key, tuple[int, int]]]:
        for first, second in combinations(range(len(auxiliary.HARD)), 2):
            first_aux, first_scale = key[first]
            second_aux, second_scale = key[second]
            old_sum = tuple(
                a + b
                for a, b in zip(
                    self.profile[first][(first_aux, first_scale)],
                    self.profile[second][(second_aux, second_scale)],
                    strict=True,
                )
            )
            for new_first_aux, new_second_aux in (
                (first_aux, second_aux),
                (second_aux, first_aux),
            ):
                if (
                    new_first_aux == auxiliary.HARD[first]
                    or new_second_aux == auxiliary.HARD[second]
                ):
                    continue
                second_lookup = self.scales_by_profile[second][new_second_aux]
                for new_first_scale in range(1, auxiliary.P):
                    first_profile = self.profile[first][
                        (new_first_aux, new_first_scale)
                    ]
                    needed = tuple(
                        total - value
                        for total, value in zip(old_sum, first_profile, strict=True)
                    )
                    for new_second_scale in second_lookup.get(needed, ()):
                        if (
                            new_first_aux,
                            new_first_scale,
                            new_second_aux,
                            new_second_scale,
                        ) == (
                            first_aux,
                            first_scale,
                            second_aux,
                            second_scale,
                        ):
                            continue
                        moved = list(key)
                        moved[first] = new_first_aux, new_first_scale
                        moved[second] = new_second_aux, new_second_scale
                        yield tuple(moved), (first, second)


def _key_sha256(key: Key) -> str:
    return hashlib.sha256(
        json.dumps(key, separators=(",", ":")).encode()
    ).hexdigest()


def _load_state(path: Path | None) -> tuple[
    list[Key], list[int | None], list[int], list[tuple[int, int] | None], deque[int]
]:
    if path is None:
        base = design_key(auxiliary.canonical_frozen_design())
        return [base], [None], [0], [None], deque([0])
    payload = json.loads(path.read_text())
    if payload.get("schema") != "residual_branch_c_auxiliary_component_v1":
        raise ValueError("resume file has the wrong schema")
    keys = [
        tuple((int(row[0]), int(row[1])) for row in key)
        for key in payload["all_design_keys"]
    ]
    parents = [
        None if value is None else int(value)
        for value in payload["parent_component_ids"]
    ]
    depths = [int(value) for value in payload["depths"]]
    moves = [
        None if value is None else (int(value[0]), int(value[1]))
        for value in payload["parent_moves"]
    ]
    queue = deque(int(value) for value in payload["frontier_component_ids"])
    if not (len(keys) == len(parents) == len(depths) == len(moves)):
        raise ValueError("resume arrays have inconsistent lengths")
    return keys, parents, depths, moves, queue


def enumerate_component(
    *, max_designs: int, resume: Path | None = None
) -> dict[str, object]:
    keys, parents, depths, moves, queue = _load_state(resume)
    if len(keys) > max_designs:
        raise ValueError("max-designs is below the resume checkpoint size")
    seen = {key: index for index, key in enumerate(keys)}
    index = PairMoveIndex()
    base = keys[0]
    base_auxiliary_set = {row[0] for row in base}
    raw = auxiliary.raw_parallel_profile(design_from_key(base))
    started = time.perf_counter()
    traversed_edges = 0
    while queue and len(keys) < max_designs:
        parent_id = queue.popleft()
        parent = keys[parent_id]
        for child, move in index.neighbors(parent):
            traversed_edges += 1
            if child in seen:
                continue
            if {row[0] for row in child} != base_auxiliary_set:
                raise ArithmeticError("a pair move changed the forced auxiliary set")
            child_id = len(keys)
            seen[child] = child_id
            keys.append(child)
            parents.append(parent_id)
            depths.append(depths[parent_id] + 1)
            moves.append(move)
            queue.append(child_id)
            if len(keys) >= max_designs:
                break

    collision_design = auxiliary.canonical_frozen_design()
    collision_cache = auxiliary._orbit_cache(collision_design)
    collision_seeds = tuple(
        seed
        for seed in auxiliary.clean_collision_seeds(
            collision_design, collision_cache
        )
        if (seed.first_half, seed.second_half) == COLLISION_PAIR
    )
    if len(collision_seeds) != auxiliary.P - 1:
        raise ArithmeticError("the factored base collision catalog changed")
    canonical_collision = min(
        collision_seeds,
        key=lambda row: (
            row.first_center, row.second_center, row.orbit
        ),
    )
    searchable = []
    for component_id, key in enumerate(keys):
        if key[COLLISION_PAIR[0]] != base[COLLISION_PAIR[0]] or key[
            COLLISION_PAIR[1]
        ] != base[COLLISION_PAIR[1]]:
            continue
        design = design_from_key(key)
        if auxiliary.raw_parallel_profile(design) != raw:
            raise ArithmeticError("a BFS vertex changed the exact raw profile")
        searchable.append(
            {
                "component_id": component_id,
                "stable_shard_parity": component_id % 2,
                "parent_component_id": parents[component_id],
                "depth": depths[component_id],
                "parent_move_half_indices": moves[component_id],
                "design_key_sha256": _key_sha256(key),
                "design_sha256": auxiliary.design_sha256(design),
                "choices": tuple(
                    {
                        "target_direction_index": choice.target_index,
                        "auxiliary_direction_index": choice.auxiliary_index,
                        "relative_scale": choice.relative_scale,
                    }
                    for choice in design
                ),
                "halves": tuple(
                    (choice.target, choice.auxiliary) for choice in design
                ),
                "collision_catalog_ref": "base_normalized_halves_3_4",
                "canonical_collision_seed": auxiliary._collision_record(
                    canonical_collision
                ),
                "fixed_direction_index": auxiliary.FIXED_DIRECTION,
                "raw_parallel_profile": raw,
            }
        )

    elapsed = time.perf_counter() - started
    return {
        "schema": "residual_branch_c_auxiliary_component_v1",
        "classification": "bounded resumable exact BFS; not an exhaustive certificate",
        "p": auxiliary.P,
        "hard_target_order": auxiliary.HARD,
        "base_design_key_sha256": _key_sha256(base),
        "base_raw_parallel_profile": raw,
        "auxiliary_direction_set": tuple(sorted(base_auxiliary_set)),
        "all_design_keys": tuple(keys),
        "parent_component_ids": tuple(parents),
        "depths": tuple(depths),
        "parent_moves": tuple(moves),
        "frontier_component_ids": tuple(queue),
        "discovered_design_count": len(keys),
        "processed_design_count": len(keys) - len(queue),
        "frontier_design_count": len(queue),
        "component_exhausted": not queue,
        "max_designs": max_designs,
        "traversed_neighbor_occurrences_this_run": traversed_edges,
        "elapsed_seconds_this_run": elapsed,
        "depth_histogram": dict(sorted(Counter(depths).items())),
        "collision_catalogs": {
            "base_normalized_halves_3_4": {
                "half_indices": COLLISION_PAIR,
                "half_choices": (
                    {
                        "target_direction_index": collision_design[
                            COLLISION_PAIR[0]
                        ].target_index,
                        "auxiliary_direction_index": collision_design[
                            COLLISION_PAIR[0]
                        ].auxiliary_index,
                        "relative_scale": collision_design[
                            COLLISION_PAIR[0]
                        ].relative_scale,
                    },
                    {
                        "target_direction_index": collision_design[
                            COLLISION_PAIR[1]
                        ].target_index,
                        "auxiliary_direction_index": collision_design[
                            COLLISION_PAIR[1]
                        ].auxiliary_index,
                        "relative_scale": collision_design[
                            COLLISION_PAIR[1]
                        ].relative_scale,
                    },
                ),
                "physical_seed_count": len(collision_seeds),
                "seeds": tuple(
                    auxiliary._collision_record(seed)
                    for seed in collision_seeds
                ),
                "canonical_seed": auxiliary._collision_record(
                    canonical_collision
                ),
            }
        },
        "searchable_design_count": len(searchable),
        "searchable_even_component_ids": sum(
            row["stable_shard_parity"] == 0 for row in searchable
        ),
        "searchable_odd_component_ids": sum(
            row["stable_shard_parity"] == 1 for row in searchable
        ),
        "searchable_designs": tuple(searchable),
        "residual_ii_closed": False,
    }


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-designs", type=int, default=10_000)
    parser.add_argument("--resume", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.max_designs < 1:
        parser.error("max-designs must be positive")
    return args


def main(argv: Iterable[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    result = enumerate_component(
        max_designs=args.max_designs, resume=args.resume
    )
    write_json_atomic(args.output, result)
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "schema",
                    "discovered_design_count",
                    "processed_design_count",
                    "frontier_design_count",
                    "component_exhausted",
                    "depth_histogram",
                    "searchable_design_count",
                    "searchable_even_component_ids",
                    "searchable_odd_component_ids",
                    "elapsed_seconds_this_run",
                )
            },
            indent=2,
            sort_keys=True,
        )
    )
    return result


if __name__ == "__main__":
    main()
