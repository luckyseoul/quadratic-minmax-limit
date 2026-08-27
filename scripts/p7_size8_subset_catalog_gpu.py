#!/usr/bin/env python3
"""Exact subset-catalog V100 sieve for the post-15.664 size-eight remainder.

The preceding omission scans leave candidate ``(boundary, stratum, leaf)``
triples that satisfy independently selected mod-3 and mod-7 equations after
all raised catalogs are omitted.  For each such triple this program restores
small, deliberately chosen subsets of the raised catalogs:

* the floor catalog paired with each raised positive catalog in stratum 11;
* every catalog pair in the three-raised-direction stratum-44 leaves; and
* every triple of raised positive catalogs in the five-direction stratum-44
  leaves.

For a variable support V and tested subset T, exact mod-7 left dependencies
are conditioned to vanish on V\\T.  Hence a zero join count for T is a
rigorous obstruction to the full catalog tuple.  Twenty-two independently
audited dependency rows are packed losslessly as one base-seven integer.  The GPU
performs exact pair/triple joins; an independent CPU prefix uses the same
integer tables.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import os
from pathlib import Path
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from p7_exceptional_omit_high_catalogs import (  # noqa: E402
    IncrementalRowBasis,
    modular_rank,
    modular_right_nullspace,
)
from p7_fixed_boundary_catalog_join import mapped_catalog  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_size8_floor_profile_gpu import direction_tables  # noqa: E402
from p7_size8_multi_elevation_tables import floor_for  # noqa: E402
from p7_size8_one_elevation_gpu import masks_for_boundary, unrank_lex  # noqa: E402
from p7_size8_remaining_gpu import allocation_patterns  # noqa: E402


MODULUS = 7
PROJECTION_ROWS = 22
MAX_TESTS = 10
EXPECTED_INTERSECTION = 181_104
EXPECTED_INTERSECTION_BY_STRATUM = (77_616, 0, 0, 103_488)


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def load_intersection(mod7_path: Path, mod3_path: Path) -> tuple[np.ndarray, dict, dict]:
    mod7 = json.loads(mod7_path.read_text())
    mod3 = json.loads(mod3_path.read_text())
    for payload, modulus in ((mod7, 7), (mod3, 3)):
        if (
            payload.get("experiment") != "p7_size8_remaining_gpu"
            or payload.get("status")
            != "complete_exact_post_15664_projected_exhaustion"
            or int(payload.get("checked_boundaries", -1)) != 450_978_066
            or int(payload.get("p", -1)) != 7
            or int(payload.get("c_H", 0)) != -1
        ):
            raise ValueError(f"invalid modulus-{modulus} omission source")
        recorded_modulus = payload.get("modulus")
        if recorded_modulus is not None and int(recorded_modulus) != modulus:
            raise ValueError(f"omission source modulus mismatch for {modulus}")
    arrays = []
    for payload in (mod7, mod3):
        value = np.asarray(
            payload["projected_survivor_rank_stratum_leaf"], dtype=np.int64
        )
        if value.ndim != 2 or value.shape[1] != 3:
            raise ValueError("malformed omission survivor triples")
        keys = (value[:, 0] << 8) | (value[:, 1] << 6) | value[:, 2]
        if len(np.unique(keys)) != len(keys):
            raise ValueError("duplicate omission survivor triple")
        arrays.append((value, keys))
    _keys, left, _right = np.intersect1d(
        arrays[0][1], arrays[1][1], assume_unique=True, return_indices=True
    )
    intersection = arrays[0][0][left]
    order = np.lexsort((intersection[:, 2], intersection[:, 1], intersection[:, 0]))
    intersection = intersection[order]
    histogram = tuple(
        int(np.count_nonzero(intersection[:, 1] == stratum)) for stratum in range(4)
    )
    if len(intersection) != EXPECTED_INTERSECTION or histogram != EXPECTED_INTERSECTION_BY_STRATUM:
        raise AssertionError(
            f"modular survivor intersection changed: {len(intersection)}, {histogram}"
        )
    return intersection, mod7, mod3


def restrict_to_candidate_source(
    base: np.ndarray,
    source_path: Path,
    mod7_path: Path,
    mod3_path: Path,
) -> tuple[np.ndarray, dict]:
    source = json.loads(source_path.read_text())
    if (
        source.get("experiment") != "p7_size8_subset_catalog_gpu"
        or source.get("status") != "complete_exact_subset_catalog_exhaustion"
        or source.get("source_mod7_sha256") != sha256(mod7_path)
        or source.get("source_mod3_sha256") != sha256(mod3_path)
        or int(source.get("p", -1)) != 7
        or int(source.get("c_H", 0)) != -1
    ):
        raise ValueError("invalid prior subset-catalog candidate source")
    candidates = np.asarray(
        source["subset_survivor_rank_stratum_leaf"], dtype=np.int64
    )
    if candidates.ndim != 2 or candidates.shape[1] != 3:
        raise ValueError("malformed prior subset-catalog survivors")
    base_keys = (base[:, 0] << 8) | (base[:, 1] << 6) | base[:, 2]
    candidate_keys = (
        (candidates[:, 0] << 8) | (candidates[:, 1] << 6) | candidates[:, 2]
    )
    if len(np.unique(candidate_keys)) != len(candidate_keys):
        raise ValueError("duplicate prior subset-catalog survivor")
    if np.any(~np.isin(candidate_keys, base_keys, assume_unique=True)):
        raise ValueError("prior subset-catalog survivor escaped the base intersection")
    order = np.lexsort((candidates[:, 2], candidates[:, 1], candidates[:, 0]))
    return candidates[order], source


def tested_subsets(
    stratum: int,
    leaf: int,
    variable_support: int,
    raised_support: int,
    epsilons: tuple[int, ...],
    test_level: str,
) -> tuple[int, ...]:
    directions = [
        direction for direction in range(8) if variable_support & (1 << direction)
    ]
    if test_level == "quad" and len(directions) == 5:
        positive_mask = sum(
            1 << direction for direction in directions if epsilons[direction] == 1
        )
        if positive_mask.bit_count() != 4:
            raise AssertionError("five-direction leaf lacks a four-positive subsystem")
        return (positive_mask,)
    if test_level == "triple" and (stratum == 0 or len(directions) == 5):
        return tuple(
            sum(1 << direction for direction in subset)
            for subset in itertools.combinations(directions, 3)
        )
    if stratum == 0:
        floor_support = variable_support ^ raised_support
        if floor_support.bit_count() != 1:
            raise AssertionError("stratum-11 candidate lacks one floor variable")
        return tuple(
            floor_support | (1 << direction)
            for direction, eps in enumerate(epsilons)
            if eps == 1
        )
    if stratum != 3:
        raise AssertionError(f"unexpected intersection stratum {stratum}")
    if len(directions) == 3:
        return tuple(
            (1 << directions[first]) | (1 << directions[second])
            for first, second in ((0, 1), (0, 2), (1, 2))
        )
    if len(directions) == 5:
        positive = [direction for direction in directions if epsilons[direction] == 1]
        if len(positive) != 4:
            raise AssertionError("five-direction leaf lacks four positive directions")
        return tuple(
            sum(1 << positive[index] for index in range(4) if index != omitted)
            for omitted in range(4)
        )
    raise AssertionError(f"unexpected stratum-44 support size {len(directions)}")


def prepare_candidates(
    triples: np.ndarray,
    labels: np.ndarray,
    epsilons: tuple[int, ...],
    allocation_supports: np.ndarray,
    test_level: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[tuple[int, int]]]:
    count = len(triples)
    masks = np.zeros((count, 8), dtype=np.uint8)
    variable_supports = np.zeros(count, dtype=np.uint8)
    context_pairs: set[tuple[int, int]] = set()
    candidate_subsets: list[tuple[int, ...]] = []
    previous_rank = -1
    previous_masks: tuple[int, ...] | None = None
    for index, (rank_value, stratum_value, leaf_value) in enumerate(triples):
        rank = int(rank_value)
        stratum = int(stratum_value)
        leaf = int(leaf_value)
        if rank != previous_rank:
            previous_masks = masks_for_boundary(unrank_lex(rank), labels)
            previous_rank = rank
        if previous_masks is None:
            raise AssertionError("boundary decoding failed")
        masks[index] = previous_masks
        raised = int(allocation_supports[stratum, leaf])
        floor_variables = [
            direction
            for direction, (eps, mask) in enumerate(zip(epsilons, previous_masks))
            if eps == -1 and mask.bit_count() == 4 and not raised & (1 << direction)
        ]
        if len(floor_variables) > 1:
            raise AssertionError("candidate has more than one floor variable catalog")
        variable = raised | sum(1 << direction for direction in floor_variables)
        variable_supports[index] = variable
        subsets = tested_subsets(
            stratum, leaf, variable, raised, epsilons, test_level
        )
        if len(subsets) > MAX_TESTS:
            raise AssertionError("candidate test count exceeds CUDA layout")
        candidate_subsets.append(subsets)
        context_pairs.update((variable, subset) for subset in subsets)

    contexts = sorted(context_pairs)
    context_index = {value: index for index, value in enumerate(contexts)}
    test_contexts = np.full((count, MAX_TESTS), -1, dtype=np.int16)
    test_masks = np.zeros((count, MAX_TESTS), dtype=np.uint8)
    for index, subsets in enumerate(candidate_subsets):
        variable = int(variable_supports[index])
        for slot, subset in enumerate(subsets):
            test_contexts[index, slot] = context_index[(variable, subset)]
            test_masks[index, slot] = subset
    return masks, variable_supports, test_contexts, test_masks, contexts


def choose_projection_rows(
    conditioned: np.ndarray, tested_mask: int
) -> tuple[np.ndarray, list[int]]:
    directions = [direction for direction in range(8) if tested_mask & (1 << direction)]
    blocks = {
        direction: conditioned[:, 2 + 35 * direction : 2 + 35 * (direction + 1)]
        for direction in directions
    }
    bases = {direction: IncrementalRowBasis(MODULUS) for direction in directions}
    selected: list[int] = []
    remaining = set(range(len(conditioned)))
    while len(selected) < PROJECTION_ROWS:
        best = None
        for index in remaining:
            gains = {
                direction: int(bases[direction].gains_rank(blocks[direction][index]))
                for direction in directions
            }
            score = (
                sum(gains.values()),
                min(
                    len(bases[direction].rows) + gains[direction]
                    for direction in directions
                ),
                sum(np.count_nonzero(blocks[direction][index]) for direction in directions),
                int(np.count_nonzero(conditioned[index])),
                -index,
            )
            if best is None or score > best[0]:
                best = (score, index)
        if best is None:
            raise AssertionError("projection row selection failed")
        index = best[1]
        selected.append(index)
        remaining.remove(index)
        for direction in directions:
            bases[direction].add(blocks[direction][index])
    return np.asarray(selected, dtype=np.int64), [
        len(bases[direction].rows) for direction in directions
    ]


def build_projections(
    contexts: list[tuple[int, int]], dependency: np.ndarray, matrix: np.ndarray
) -> tuple[np.ndarray, list[dict]]:
    projections = np.zeros((len(contexts), PROJECTION_ROWS, 282), dtype=np.uint8)
    records = []
    for context, (variable_mask, tested_mask) in enumerate(contexts):
        omitted_mask = variable_mask & ~tested_mask
        columns = np.asarray(
            [
                column
                for direction in range(8)
                if omitted_mask & (1 << direction)
                for column in range(2 + 35 * direction, 2 + 35 * (direction + 1))
            ],
            dtype=np.int64,
        )
        coefficients, block_rank = modular_right_nullspace(
            dependency[:, columns].astype(np.int64).T, MODULUS
        )
        conditioned = coefficients @ dependency.astype(np.int64) % MODULUS
        selected, block_ranks = choose_projection_rows(conditioned, tested_mask)
        projection = conditioned[selected].astype(np.uint8)
        if (
            modular_rank(projection, MODULUS) != PROJECTION_ROWS
            or np.any(projection[:, columns])
            or np.any(
                projection.astype(np.int64)
                @ (matrix.astype(np.int64) % MODULUS)
                % MODULUS
            )
        ):
            raise AssertionError(f"invalid subset projection {(variable_mask, tested_mask)}")
        projections[context] = projection
        records.append(
            {
                "context_index": context,
                "variable_support_mask": variable_mask,
                "tested_subset_mask": tested_mask,
                "tested_subset_size": tested_mask.bit_count(),
                "conditioned_dimension": int(len(conditioned)),
                "omitted_block_rank": int(block_rank),
                "selected_rank": PROJECTION_ROWS,
                "selected_tested_block_ranks": block_ranks,
            }
        )
    return projections, records


def pack_rows(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.uint64)
    if values.ndim != 2 or values.shape[1] != PROJECTION_ROWS:
        raise ValueError("signature rows have the wrong shape")
    powers = np.asarray(
        [MODULUS**index for index in range(PROJECTION_ROWS)], dtype=np.uint64
    )
    return (values @ powers).astype(np.uint64)


def increment_index(value: int) -> int:
    if value == 0:
        return 0
    if value == 4:
        return 1
    if value == 8:
        return 2
    raise AssertionError(f"unsupported subset increment {value}")


def build_tables(
    projections: np.ndarray,
    contexts: list[tuple[int, int]],
    masks: np.ndarray,
    triples: np.ndarray,
    test_contexts: np.ndarray,
    test_masks: np.ndarray,
    epsilons: tuple[int, ...],
    allocation_increments: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict]:
    context_count = len(contexts)
    base = np.zeros((context_count, PROJECTION_ROWS), dtype=np.uint8)
    singleton = np.full(
        (context_count, 8, 128, PROJECTION_ROWS), 255, dtype=np.uint8
    )
    for context, (variable_mask, _tested_mask) in enumerate(contexts):
        projection = projections[context].astype(np.int64)
        base[context] = (
            projection[:, :2] @ np.asarray([29, 1], dtype=np.int64) % MODULUS
        ).astype(np.uint8)
        for direction, eps in enumerate(epsilons):
            if variable_mask & (1 << direction):
                continue
            block = projection[:, 2 + 35 * direction : 2 + 35 * (direction + 1)]
            for mask in range(128):
                if mask.bit_count() not in (0, 2, 4, 6):
                    continue
                odd = {value for value in range(7) if mask & (1 << value)}
                values = mapped_catalog(
                    mask.bit_count(), int(eps == -1), floor_for(eps, mask.bit_count()), odd, None
                ).astype(np.int64)
                if len(values) == 1:
                    singleton[context, direction, mask] = (
                        block @ (13 - values[0]) % MODULUS
                    ).astype(np.uint8)

    active_keys: set[tuple[int, int, int, int]] = set()
    for candidate, (_rank, stratum_value, leaf_value) in enumerate(triples):
        stratum = int(stratum_value)
        leaf = int(leaf_value)
        for slot in range(MAX_TESTS):
            context = int(test_contexts[candidate, slot])
            if context < 0:
                continue
            tested_mask = int(test_masks[candidate, slot])
            variable_mask = contexts[context][0]
            for direction in range(8):
                if not tested_mask & (1 << direction):
                    continue
                increment = (
                    int(allocation_increments[stratum, leaf, direction])
                    if variable_mask & (1 << direction)
                    else 0
                )
                active_keys.add(
                    (context, direction, increment_index(increment), int(masks[candidate, direction]))
                )

    offsets = np.zeros((context_count, 8, 3, 128), dtype=np.uint64)
    counts = np.zeros((context_count, 8, 3, 128), dtype=np.uint32)
    signatures: list[np.ndarray] = []
    cursor = 0
    raw_histogram: Counter[int] = Counter()
    distinct_histogram: Counter[int] = Counter()
    started = time.time()
    for number, (context, direction, inc_index, mask) in enumerate(sorted(active_keys), 1):
        increment = (0, 4, 8)[inc_index]
        eps = epsilons[direction]
        odd = {value for value in range(7) if mask & (1 << value)}
        values = mapped_catalog(
            mask.bit_count(),
            int(eps == -1),
            floor_for(eps, mask.bit_count()) + increment,
            odd,
            None,
        ).astype(np.int64)
        block = projections[
            context, :, 2 + 35 * direction : 2 + 35 * (direction + 1)
        ].astype(np.int64)
        contribution = (block @ (13 - values).T % MODULUS).T.astype(np.uint8)
        packed = np.unique(pack_rows(contribution))
        offsets[context, direction, inc_index, mask] = cursor
        counts[context, direction, inc_index, mask] = len(packed)
        signatures.append(packed)
        cursor += len(packed)
        raw_histogram[len(values)] += 1
        distinct_histogram[len(packed)] += 1
        if number % 200 == 0 or number == len(active_keys):
            print(
                json.dumps(
                    {
                        "table_progress": [number, len(active_keys)],
                        "packed_signatures": cursor,
                        "elapsed_seconds": time.time() - started,
                    }
                ),
                flush=True,
            )
    flattened = (
        np.concatenate(signatures).astype(np.uint64, copy=False)
        if signatures
        else np.zeros(0, dtype=np.uint64)
    )
    arrays = {
        "base": base,
        "singleton": singleton,
        "offsets": offsets,
        "counts": counts,
        "signatures": flattened,
    }
    summary = {
        "active_catalog_keys": len(active_keys),
        "raw_catalog_size_histogram": dict(sorted(raw_histogram.items())),
        "distinct_signature_size_histogram": dict(sorted(distinct_histogram.items())),
        "packed_signature_count": int(len(flattened)),
        "array_shapes": {key: list(value.shape) for key, value in arrays.items()},
        "array_sha256": {key: array_sha256(value) for key, value in arrays.items()},
        "elapsed_seconds": time.time() - started,
    }
    return arrays, summary


def modular_sum_signatures(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    left = np.repeat(first.astype(np.uint64, copy=False), len(second)).copy()
    right = np.tile(second.astype(np.uint64, copy=False), len(first)).copy()
    output = np.zeros(len(left), dtype=np.uint64)
    place = np.uint64(1)
    for _row in range(PROJECTION_ROWS):
        output += place * ((left % MODULUS + right % MODULUS) % MODULUS)
        left //= MODULUS
        right //= MODULUS
        place *= np.uint64(MODULUS)
    return np.unique(output)


def build_quad_pair_tables(
    triples: np.ndarray,
    masks: np.ndarray,
    test_contexts: np.ndarray,
    test_masks: np.ndarray,
    allocation_increments: np.ndarray,
    tables: dict[str, np.ndarray],
) -> tuple[np.ndarray, dict]:
    pair_keys: set[tuple[int, int, int, int, int, int, int]] = set()
    candidate_pair_keys: dict[tuple[int, int, int], tuple[int, ...]] = {}
    for candidate, (_rank, stratum_value, leaf_value) in enumerate(triples):
        stratum = int(stratum_value)
        leaf = int(leaf_value)
        for slot in range(MAX_TESTS):
            context = int(test_contexts[candidate, slot])
            tested_mask = int(test_masks[candidate, slot])
            if context < 0 or tested_mask.bit_count() != 4:
                continue
            directions = [
                direction for direction in range(8) if tested_mask & (1 << direction)
            ]
            for pair_slot, pair in enumerate((directions[:2], directions[2:])):
                fields: list[int] = [context]
                for direction in pair:
                    increment = int(
                        allocation_increments[stratum, leaf, direction]
                    )
                    fields.extend(
                        (
                            direction,
                            increment_index(increment),
                            int(masks[candidate, direction]),
                        )
                    )
                key = tuple(fields)
                if len(key) != 7:
                    raise AssertionError("quad pair key layout changed")
                pair_keys.add(key)
                candidate_pair_keys[(candidate, slot, pair_slot)] = key

    ordered_keys = sorted(pair_keys)
    key_index = {key: index for index, key in enumerate(ordered_keys)}
    candidate_indices = np.full((len(triples), MAX_TESTS, 2), -1, dtype=np.int32)
    for location, key in candidate_pair_keys.items():
        candidate_indices[location] = key_index[key]
    pair_offsets = np.zeros(len(ordered_keys), dtype=np.uint64)
    pair_counts = np.zeros(len(ordered_keys), dtype=np.uint32)
    pieces = []
    cursor = 0
    size_histogram: Counter[int] = Counter()
    started = time.time()
    for number, key in enumerate(ordered_keys, 1):
        context, d0, inc0, mask0, d1, inc1, mask1 = key
        ranges = []
        for direction, inc_index, mask in (
            (d0, inc0, mask0),
            (d1, inc1, mask1),
        ):
            offset = int(tables["offsets"][context, direction, inc_index, mask])
            count = int(tables["counts"][context, direction, inc_index, mask])
            if count == 0:
                raise AssertionError("quad pair references an empty catalog signature")
            ranges.append(tables["signatures"][offset : offset + count])
        sums = modular_sum_signatures(ranges[0], ranges[1])
        pair_offsets[number - 1] = cursor
        pair_counts[number - 1] = len(sums)
        pieces.append(sums)
        cursor += len(sums)
        size_histogram[len(sums)] += 1
        if number % 500 == 0 or number == len(ordered_keys):
            print(
                json.dumps(
                    {
                        "quad_pair_progress": [number, len(ordered_keys)],
                        "packed_pair_signatures": cursor,
                        "elapsed_seconds": time.time() - started,
                    }
                ),
                flush=True,
            )
    pair_signatures = (
        np.concatenate(pieces).astype(np.uint64, copy=False)
        if pieces
        else np.zeros(0, dtype=np.uint64)
    )
    tables.update(
        {
            "pair_offsets": pair_offsets,
            "pair_counts": pair_counts,
            "pair_signatures": pair_signatures,
        }
    )
    return candidate_indices, {
        "pair_key_count": len(ordered_keys),
        "pair_signature_count": int(len(pair_signatures)),
        "pair_signature_size_histogram": dict(sorted(size_histogram.items())),
        "candidate_pair_indices_sha256": array_sha256(candidate_indices),
        "pair_offsets_sha256": array_sha256(pair_offsets),
        "pair_counts_sha256": array_sha256(pair_counts),
        "pair_signatures_sha256": array_sha256(pair_signatures),
        "elapsed_seconds": time.time() - started,
    }


CUDA_SOURCE = r'''
extern "C" {

__device__ __forceinline__ unsigned long long needed2(
    unsigned long long target, unsigned long long a)
{
    unsigned long long output = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int row = 0; row < 22; ++row) {
        const unsigned int t = (unsigned int)(target % 7ULL);
        const unsigned int x = (unsigned int)(a % 7ULL);
        output += place * (unsigned long long)((t + 14U - x) % 7U);
        target /= 7ULL;
        a /= 7ULL;
        place *= 7ULL;
    }
    return output;
}

__device__ __forceinline__ unsigned long long needed3(
    unsigned long long target, unsigned long long a, unsigned long long b)
{
    unsigned long long output = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int row = 0; row < 22; ++row) {
        const unsigned int t = (unsigned int)(target % 7ULL);
        const unsigned int x = (unsigned int)(a % 7ULL);
        const unsigned int y = (unsigned int)(b % 7ULL);
        output += place * (unsigned long long)((t + 21U - x - y) % 7U);
        target /= 7ULL;
        a /= 7ULL;
        b /= 7ULL;
        place *= 7ULL;
    }
    return output;
}

__device__ __forceinline__ bool contains_sorted(
    const unsigned long long* signatures,
    unsigned long long offset,
    unsigned int count,
    unsigned long long target)
{
    unsigned int low = 0U;
    unsigned int high = count;
    while (low < high) {
        const unsigned int middle = low + ((high - low) >> 1);
        const unsigned long long value = signatures[offset + middle];
        if (value < target) low = middle + 1U;
        else high = middle;
    }
    return low < count && signatures[offset + low] == target;
}

__global__ void subset_catalog_scan(
    const unsigned char* masks,
    const unsigned char* variable_supports,
    const short* test_contexts,
    const unsigned char* test_masks,
    const unsigned char* strata,
    const unsigned char* leaves,
    const unsigned char* allocation_increments,
    const unsigned char* base,
    const unsigned char* singleton,
    const unsigned long long* offsets,
    const unsigned int* counts,
    const unsigned long long* signatures,
    const int* quad_pair_indices,
    const unsigned long long* pair_offsets,
    const unsigned int* pair_counts,
    const unsigned long long* pair_signatures,
    const unsigned long long candidate_count,
    unsigned long long* survivor_indices,
    unsigned long long* survivor_count,
    unsigned long long* failure_counts,
    const unsigned long long survivor_capacity)
{
    const unsigned long long candidate =
        (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (candidate >= candidate_count) return;
    const int variable_mask = (int)variable_supports[candidate];
    bool survives = true;
    for (int slot = 0; slot < 10 && survives; ++slot) {
        const int context = (int)test_contexts[candidate * 10 + slot];
        if (context < 0) continue;
        const int tested_mask = (int)test_masks[candidate * 10 + slot];
        int syndrome[22];
        #pragma unroll
        for (int row = 0; row < 22; ++row)
            syndrome[row] = (int)base[context * 22 + row];
        bool valid = true;
        #pragma unroll
        for (int direction = 0; direction < 8; ++direction) {
            if (variable_mask & (1 << direction)) continue;
            const int mask = (int)masks[candidate * 8 + direction];
            const unsigned long long key =
                ((unsigned long long)context * 8ULL + (unsigned long long)direction)
                * 128ULL + (unsigned long long)mask;
            const unsigned long long start = key * 22ULL;
            #pragma unroll
            for (int row = 0; row < 22; ++row) {
                const int value = (int)singleton[start + row];
                if (value == 255) valid = false;
                syndrome[row] += value;
            }
        }
        unsigned long long target = 0ULL;
        unsigned long long place = 1ULL;
        #pragma unroll
        for (int row = 0; row < 22; ++row) {
            const unsigned int digit =
                (unsigned int)((7 - (syndrome[row] % 7)) % 7);
            target += place * (unsigned long long)digit;
            place *= 7ULL;
        }
        int directions[3];
        unsigned long long starts[3];
        unsigned int sizes[3];
        int tested_count = 0;
        const int stratum = (int)strata[candidate];
        const int leaf = (int)leaves[candidate];
        #pragma unroll
        for (int direction = 0; direction < 8; ++direction) {
            if (!(tested_mask & (1 << direction))) continue;
            const int increment = (int)allocation_increments[
                ((stratum * 44 + leaf) * 8) + direction
            ];
            const int increment_index = increment == 0 ? 0 : (increment == 4 ? 1 : 2);
            const int mask = (int)masks[candidate * 8 + direction];
            const unsigned long long key =
                (((unsigned long long)context * 8ULL + (unsigned long long)direction)
                 * 3ULL + (unsigned long long)increment_index)
                * 128ULL + (unsigned long long)mask;
            directions[tested_count] = direction;
            starts[tested_count] = offsets[key];
            sizes[tested_count] = counts[key];
            if (sizes[tested_count] == 0U) valid = false;
            ++tested_count;
        }
        bool match = false;
        if (valid && tested_count == 2) {
            int enumerate = sizes[0] <= sizes[1] ? 0 : 1;
            int search = 1 - enumerate;
            for (unsigned int i = 0; i < sizes[enumerate] && !match; ++i) {
                const unsigned long long needed = needed2(
                    target, signatures[starts[enumerate] + i]
                );
                match = contains_sorted(
                    signatures, starts[search], sizes[search], needed
                );
            }
        } else if (valid && tested_count == 3) {
            int search = 0;
            if (sizes[1] > sizes[search]) search = 1;
            if (sizes[2] > sizes[search]) search = 2;
            const int first = search == 0 ? 1 : 0;
            const int second = search == 2 ? 1 : 2;
            for (unsigned int i = 0; i < sizes[first] && !match; ++i) {
                const unsigned long long a = signatures[starts[first] + i];
                for (unsigned int j = 0; j < sizes[second] && !match; ++j) {
                    const unsigned long long needed = needed3(
                        target, a, signatures[starts[second] + j]
                    );
                    match = contains_sorted(
                        signatures, starts[search], sizes[search], needed
                    );
                }
            }
        } else if (valid && tested_count == 4) {
            const unsigned long long pair_key =
                (candidate * 10ULL + (unsigned long long)slot) * 2ULL;
            const int pair0 = quad_pair_indices[pair_key];
            const int pair1 = quad_pair_indices[pair_key + 1ULL];
            if (pair0 >= 0 && pair1 >= 0) {
                const int enumerate =
                    pair_counts[pair0] <= pair_counts[pair1] ? pair0 : pair1;
                const int search = enumerate == pair0 ? pair1 : pair0;
                for (unsigned int i = 0; i < pair_counts[enumerate] && !match; ++i) {
                    const unsigned long long needed = needed2(
                        target, pair_signatures[pair_offsets[enumerate] + i]
                    );
                    match = contains_sorted(
                        pair_signatures,
                        pair_offsets[search],
                        pair_counts[search],
                        needed
                    );
                }
            }
        }
        if (!match) {
            survives = false;
            atomicAdd(failure_counts + context, 1ULL);
        }
    }
    if (survives) {
        const unsigned long long slot = atomicAdd(survivor_count, 1ULL);
        if (slot < survivor_capacity) survivor_indices[slot] = candidate;
    }
}

}
'''


def needed2(target: int, value: int) -> int:
    output = 0
    place = 1
    for row in range(PROJECTION_ROWS):
        output += place * (((target % MODULUS) - (value % MODULUS)) % MODULUS)
        target //= MODULUS
        value //= MODULUS
        place *= MODULUS
    return output


def needed3(target: int, first: int, second: int) -> int:
    output = 0
    place = 1
    for row in range(PROJECTION_ROWS):
        output += place * (
            (
                (target % MODULUS)
                - (first % MODULUS)
                - (second % MODULUS)
            )
            % MODULUS
        )
        target //= MODULUS
        first //= MODULUS
        second //= MODULUS
        place *= MODULUS
    return output


def cpu_candidate_passes(
    candidate: int,
    masks: np.ndarray,
    variable_supports: np.ndarray,
    test_contexts: np.ndarray,
    test_masks: np.ndarray,
    strata: np.ndarray,
    leaves: np.ndarray,
    allocation_increments: np.ndarray,
    tables: dict[str, np.ndarray],
    quad_pair_indices: np.ndarray,
) -> bool:
    variable_mask = int(variable_supports[candidate])
    for slot in range(MAX_TESTS):
        context = int(test_contexts[candidate, slot])
        if context < 0:
            continue
        syndrome = tables["base"][context].astype(np.int16).copy()
        for direction in range(8):
            if variable_mask & (1 << direction):
                continue
            value = tables["singleton"][context, direction, int(masks[candidate, direction])]
            if np.any(value == 255):
                return False
            syndrome += value
        target = int(pack_rows(((-syndrome % MODULUS).astype(np.uint8))[None, :])[0])
        catalog_ranges = []
        tested_mask = int(test_masks[candidate, slot])
        for direction in range(8):
            if not tested_mask & (1 << direction):
                continue
            increment = int(
                allocation_increments[
                    int(strata[candidate]), int(leaves[candidate]), direction
                ]
            )
            inc_index = increment_index(increment)
            mask = int(masks[candidate, direction])
            offset = int(tables["offsets"][context, direction, inc_index, mask])
            count = int(tables["counts"][context, direction, inc_index, mask])
            catalog_ranges.append(tables["signatures"][offset : offset + count])
        match = False
        if len(catalog_ranges) == 2:
            first, second = sorted(catalog_ranges, key=len)
            for value in first:
                needed = needed2(target, int(value))
                index = int(np.searchsorted(second, needed))
                if index < len(second) and int(second[index]) == needed:
                    match = True
                    break
        elif len(catalog_ranges) == 3:
            search = max(range(3), key=lambda index: len(catalog_ranges[index]))
            others = [index for index in range(3) if index != search]
            target_catalog = catalog_ranges[search]
            for first in catalog_ranges[others[0]]:
                for second in catalog_ranges[others[1]]:
                    needed = needed3(target, int(first), int(second))
                    index = int(np.searchsorted(target_catalog, needed))
                    if index < len(target_catalog) and int(target_catalog[index]) == needed:
                        match = True
                        break
                if match:
                    break
        elif len(catalog_ranges) == 4:
            pair_indices = quad_pair_indices[candidate, slot]
            if np.any(pair_indices < 0):
                raise AssertionError("quad test lacks pair-signature tables")
            pair_ranges = []
            for pair_index in pair_indices:
                offset = int(tables["pair_offsets"][pair_index])
                count = int(tables["pair_counts"][pair_index])
                pair_ranges.append(
                    tables["pair_signatures"][offset : offset + count]
                )
            first, second = sorted(pair_ranges, key=len)
            for value in first:
                needed = needed2(target, int(value))
                index = int(np.searchsorted(second, needed))
                if index < len(second) and int(second[index]) == needed:
                    match = True
                    break
        else:
            raise AssertionError("subset test is not a pair, triple, or quadruple")
        if not match:
            return False
    return True


def gpu_scan(
    kernel,
    triples: np.ndarray,
    masks: np.ndarray,
    variable_supports: np.ndarray,
    test_contexts: np.ndarray,
    test_masks: np.ndarray,
    allocation_increments: np.ndarray,
    tables: dict[str, np.ndarray],
    quad_pair_indices: np.ndarray,
    context_count: int,
    blocks: int,
    threads: int,
) -> dict:
    import cupy as cp

    count = len(triples)
    gpu = {
        "masks": cp.asarray(masks),
        "variable_supports": cp.asarray(variable_supports),
        "test_contexts": cp.asarray(test_contexts),
        "test_masks": cp.asarray(test_masks),
        "strata": cp.asarray(triples[:, 1].astype(np.uint8)),
        "leaves": cp.asarray(triples[:, 2].astype(np.uint8)),
        "allocation_increments": cp.asarray(allocation_increments),
        "quad_pair_indices": cp.asarray(quad_pair_indices),
        **{key: cp.asarray(value) for key, value in tables.items()},
    }
    survivor_indices = cp.zeros(count, dtype=cp.uint64)
    survivor_count = cp.zeros(1, dtype=cp.uint64)
    failure_counts = cp.zeros(context_count, dtype=cp.uint64)
    started = time.perf_counter()
    kernel(
        (blocks,),
        (threads,),
        (
            gpu["masks"],
            gpu["variable_supports"],
            gpu["test_contexts"],
            gpu["test_masks"],
            gpu["strata"],
            gpu["leaves"],
            gpu["allocation_increments"],
            gpu["base"],
            gpu["singleton"],
            gpu["offsets"],
            gpu["counts"],
            gpu["signatures"],
            gpu["quad_pair_indices"],
            gpu["pair_offsets"],
            gpu["pair_counts"],
            gpu["pair_signatures"],
            np.uint64(count),
            survivor_indices,
            survivor_count,
            failure_counts,
            np.uint64(count),
        ),
    )
    cp.cuda.runtime.deviceSynchronize()
    elapsed = time.perf_counter() - started
    total = int(cp.asnumpy(survivor_count)[0])
    indices = cp.asnumpy(survivor_indices[:total]).astype(np.int64)
    failures = cp.asnumpy(failure_counts).astype(np.int64)
    return {
        "survivor_indices": np.sort(indices),
        "failure_counts": failures,
        "elapsed_seconds": elapsed,
    }


def run(args: argparse.Namespace) -> dict:
    import cupy as cp

    started = time.time()
    base_triples, mod7_source, mod3_source = load_intersection(
        args.mod7_source, args.mod3_source
    )
    prior_source = None
    triples = base_triples
    if args.candidate_source is not None:
        triples, prior_source = restrict_to_candidate_source(
            base_triples,
            args.candidate_source,
            args.mod7_source,
            args.mod3_source,
        )
    labels, epsilon_array = direction_tables()
    epsilons = tuple(int(value) for value in epsilon_array)
    allocation_counts, allocation_supports, allocation_increments = allocation_patterns(
        epsilons
    )
    masks, variable_supports, test_contexts, test_masks, contexts = prepare_candidates(
        triples, labels, epsilons, allocation_supports, args.test_level
    )
    matrix, dependencies, linear_rows = linear_data((MODULUS,))
    dependency = dependencies[MODULUS].astype(np.uint8)
    projections, projection_records = build_projections(contexts, dependency, matrix)
    tables, table_summary = build_tables(
        projections,
        contexts,
        masks,
        triples,
        test_contexts,
        test_masks,
        epsilons,
        allocation_increments,
    )
    quad_pair_indices, quad_pair_summary = build_quad_pair_tables(
        triples,
        masks,
        test_contexts,
        test_masks,
        allocation_increments,
        tables,
    )
    table_summary["quad_pair_tables"] = quad_pair_summary
    table_summary["array_shapes"] = {
        key: list(value.shape) for key, value in tables.items()
    }
    table_summary["array_sha256"] = {
        key: array_sha256(value) for key, value in tables.items()
    }
    kernel = cp.RawKernel(CUDA_SOURCE, "subset_catalog_scan", options=("--std=c++11",))
    prefix_count = min(args.cpu_prefix, len(triples))
    cpu_prefix_survivors = [
        index
        for index in range(prefix_count)
        if cpu_candidate_passes(
            index,
            masks,
            variable_supports,
            test_contexts,
            test_masks,
            triples[:, 1].astype(np.uint8),
            triples[:, 2].astype(np.uint8),
            allocation_increments,
            tables,
            quad_pair_indices,
        )
    ]
    prefix_gpu = gpu_scan(
        kernel,
        triples[:prefix_count],
        masks[:prefix_count],
        variable_supports[:prefix_count],
        test_contexts[:prefix_count],
        test_masks[:prefix_count],
        allocation_increments,
        tables,
        quad_pair_indices[:prefix_count],
        len(contexts),
        args.blocks,
        args.threads,
    )
    if cpu_prefix_survivors != prefix_gpu["survivor_indices"].tolist():
        gpu_prefix_survivors = prefix_gpu["survivor_indices"].tolist()
        cpu_set = set(cpu_prefix_survivors)
        gpu_set = set(gpu_prefix_survivors)
        print(
            json.dumps(
                {
                    "cpu_gpu_prefix_mismatch": True,
                    "cpu_survivor_count": len(cpu_prefix_survivors),
                    "gpu_survivor_count": len(gpu_prefix_survivors),
                    "cpu_only": sorted(cpu_set - gpu_set)[:20],
                    "gpu_only": sorted(gpu_set - cpu_set)[:20],
                }
            ),
            flush=True,
        )
        raise AssertionError("CPU/GPU subset-catalog prefix mismatch")
    full = gpu_scan(
        kernel,
        triples,
        masks,
        variable_supports,
        test_contexts,
        test_masks,
        allocation_increments,
        tables,
        quad_pair_indices,
        len(contexts),
        args.blocks,
        args.threads,
    )
    survivors = triples[full["survivor_indices"]]
    survivor_histogram = [
        int(np.count_nonzero(survivors[:, 1] == stratum)) for stratum in range(4)
    ]
    failure_rows = [
        {
            **projection_records[index],
            "first_failed_candidate_count": int(full["failure_counts"][index]),
        }
        for index in range(len(contexts))
    ]
    out = {
        "experiment": "p7_size8_subset_catalog_gpu",
        "status": "complete_exact_subset_catalog_exhaustion",
        "p": 7,
        "c_H": -1,
        "modulus": MODULUS,
        "projection_rows": PROJECTION_ROWS,
        "source_mod7": str(args.mod7_source),
        "source_mod7_sha256": sha256(args.mod7_source),
        "source_mod3": str(args.mod3_source),
        "source_mod3_sha256": sha256(args.mod3_source),
        "candidate_source": (
            str(args.candidate_source) if args.candidate_source is not None else None
        ),
        "candidate_source_sha256": (
            sha256(args.candidate_source) if args.candidate_source is not None else None
        ),
        "candidate_source_survivors": (
            int(prior_source["subset_survivor_count"])
            if prior_source is not None
            else None
        ),
        "source_mod7_survivors": len(
            mod7_source["projected_survivor_rank_stratum_leaf"]
        ),
        "source_mod3_survivors": len(
            mod3_source["projected_survivor_rank_stratum_leaf"]
        ),
        "base_intersected_candidate_count": len(base_triples),
        "intersected_candidate_count": len(triples),
        "intersection_by_stratum": list(EXPECTED_INTERSECTION_BY_STRATUM),
        "test_level": args.test_level,
        "mathematical_reduction": {
            "conditioned_equations": "left dependencies vanish on V\\T",
            "stratum_11_tests": (
                "the four-raised-positive-catalog subsystem"
                if args.test_level == "quad"
                else "all ten catalog triples"
                if args.test_level == "triple"
                else "floor variable paired with each of four raised positive catalogs"
            ),
            "stratum_44_three_support_tests": "all three catalog pairs",
            "stratum_44_five_support_tests": (
                "the four-raised-positive-catalog subsystem"
                if args.test_level == "quad"
                else "all ten catalog triples"
                if args.test_level == "triple"
                else "all four triples of raised positive catalogs"
            ),
            "packed_signature": "22 exact mod-7 digits in one base-seven uint64",
            "logical_implication": "zero subset join implies no full catalog tuple",
        },
        "linear_system": linear_rows,
        "dependency_sha256": array_sha256(dependency),
        "projection_sha256": array_sha256(projections),
        "context_count": len(contexts),
        "contexts": failure_rows,
        "table_build": table_summary,
        "verification": {
            "cpu_prefix_candidates": prefix_count,
            "cpu_prefix_survivors": len(cpu_prefix_survivors),
            "gpu_prefix_survivors": len(prefix_gpu["survivor_indices"]),
            "cpu_gpu_prefix_exact_match": True,
            "projection_left_null_audit": True,
            "projection_omitted_blocks_zero_audit": True,
            "projection_rank_audit": True,
        },
        "device": {
            "name": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
            "compute_capability": [
                cp.cuda.Device(0).compute_capability[0],
                cp.cuda.Device(0).compute_capability[1],
            ],
        },
        "launch": {"blocks": args.blocks, "threads": args.threads},
        "gpu_scan_seconds": full["elapsed_seconds"],
        "subset_survivor_count": len(survivors),
        "subset_survivors_by_stratum": survivor_histogram,
        "subset_survivor_rank_stratum_leaf": [list(map(int, row)) for row in survivors],
        "all_intersected_candidates_excluded": len(survivors) == 0,
        "closes_cminus1_post_15664_scope": len(survivors) == 0,
        "elapsed_seconds": time.time() - started,
    }
    atomic_json(args.output, out)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mod7-source", type=Path, required=True)
    parser.add_argument("--mod3-source", type=Path, required=True)
    parser.add_argument("--candidate-source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--test-level", choices=("local", "triple", "quad"), default="local"
    )
    parser.add_argument("--cpu-prefix", type=int, default=256)
    parser.add_argument("--blocks", type=int, default=640)
    parser.add_argument("--threads", type=int, default=256)
    args = parser.parse_args()
    output = run(args)
    print(
        json.dumps(
            {
                "status": output["status"],
                "intersected_candidate_count": output["intersected_candidate_count"],
                "subset_survivor_count": output["subset_survivor_count"],
                "subset_survivors_by_stratum": output["subset_survivors_by_stratum"],
                "all_intersected_candidates_excluded": output[
                    "all_intersected_candidates_excluded"
                ],
                "gpu_scan_seconds": output["gpu_scan_seconds"],
                "elapsed_seconds": output["elapsed_seconds"],
                "output": str(args.output),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
