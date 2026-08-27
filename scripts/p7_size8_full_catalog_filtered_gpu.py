#!/usr/bin/env python3
"""Exact filtered full-catalog V100 join for residual p=7 size-eight leaves.

Each variable catalog is first filtered by exact mod-seven dependencies that
vanish on every other variable direction.  The surviving rows retain their
signatures under one common 22-row full dependency projection.  A CUDA block
then performs a complete meet-in-the-middle join of all three or five variable
catalogs in a shared-memory exact hash table.  Zero matches are a rigorous
obstruction to the complete catalog tuple; nonzero matches remain unresolved.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
import os
from pathlib import Path
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from p7_exceptional_omit_high_catalogs import (  # noqa: E402
    modular_rank,
    modular_right_nullspace,
)
from p7_fixed_boundary_catalog_join import mapped_catalog  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_size8_floor_profile_gpu import direction_tables  # noqa: E402
from p7_size8_multi_elevation_tables import floor_for  # noqa: E402
from p7_size8_one_elevation_gpu import masks_for_boundary, unrank_lex  # noqa: E402
from p7_size8_remaining_gpu import allocation_patterns  # noqa: E402
from p7_size8_subset_catalog_gpu import (  # noqa: E402
    MODULUS,
    PROJECTION_ROWS,
    array_sha256,
    atomic_json,
    choose_projection_rows,
    increment_index,
    load_intersection,
    pack_rows,
    restrict_to_candidate_source,
    sha256,
)


HASH_CAPACITY = 4096
MAX_HASH_BUILD = 1800
MAX_VARIABLES = 5


def prepare_candidates(
    triples: np.ndarray,
    labels: np.ndarray,
    epsilons: tuple[int, ...],
    allocation_supports: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list[tuple[int, int]]]:
    masks = np.zeros((len(triples), 8), dtype=np.uint8)
    variable_supports = np.zeros(len(triples), dtype=np.uint8)
    contexts: set[tuple[int, int]] = set()
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
            raise AssertionError("candidate has multiple non-raised floor catalogs")
        variable = raised | sum(1 << direction for direction in floor_variables)
        if variable.bit_count() not in (3, 5):
            raise AssertionError(f"unexpected full-join variable count {variable.bit_count()}")
        variable_supports[index] = variable
        contexts.update(
            (variable, direction)
            for direction in range(8)
            if variable & (1 << direction)
        )
    return masks, variable_supports, sorted(contexts)


def build_projections(
    contexts: list[tuple[int, int]], dependency: np.ndarray, matrix: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict]]:
    selected, full_block_ranks = choose_projection_rows(dependency, 255)
    full = dependency[selected].astype(np.uint8)
    if (
        modular_rank(full, MODULUS) != PROJECTION_ROWS
        or np.any(full.astype(np.int64) @ (matrix.astype(np.int64) % MODULUS) % MODULUS)
    ):
        raise AssertionError("full dependency projection audit failed")
    isolates = np.zeros((len(contexts), PROJECTION_ROWS, 282), dtype=np.uint8)
    context_index = np.full((256, 8), -1, dtype=np.int16)
    records = []
    for context, (variable_mask, direction) in enumerate(contexts):
        context_index[variable_mask, direction] = context
        omitted_mask = variable_mask & ~(1 << direction)
        columns = np.asarray(
            [
                column
                for omitted in range(8)
                if omitted_mask & (1 << omitted)
                for column in range(2 + 35 * omitted, 2 + 35 * (omitted + 1))
            ],
            dtype=np.int64,
        )
        coefficients, block_rank = modular_right_nullspace(
            dependency[:, columns].astype(np.int64).T, MODULUS
        )
        conditioned = coefficients @ dependency.astype(np.int64) % MODULUS
        rows, tested_ranks = choose_projection_rows(conditioned, 1 << direction)
        projection = conditioned[rows].astype(np.uint8)
        if (
            modular_rank(projection, MODULUS) != PROJECTION_ROWS
            or np.any(projection[:, columns])
            or np.any(
                projection.astype(np.int64)
                @ (matrix.astype(np.int64) % MODULUS)
                % MODULUS
            )
        ):
            raise AssertionError(f"isolate projection audit failed for {contexts[context]}")
        isolates[context] = projection
        records.append(
            {
                "context_index": context,
                "variable_support_mask": variable_mask,
                "isolated_direction": direction,
                "conditioned_dimension": int(len(conditioned)),
                "omitted_block_rank": int(block_rank),
                "selected_rank": PROJECTION_ROWS,
                "isolated_block_rank": int(tested_ranks[0]),
            }
        )
    return full, isolates, context_index, [
        {
            "full_projection_rank": PROJECTION_ROWS,
            "full_direction_block_ranks": full_block_ranks,
        },
        *records,
    ]


def singleton_tables(
    full: np.ndarray,
    isolates: np.ndarray,
    contexts: list[tuple[int, int]],
    epsilons: tuple[int, ...],
) -> dict[str, np.ndarray]:
    full_base = (
        full[:, :2].astype(np.int64) @ np.asarray([29, 1], dtype=np.int64) % MODULUS
    ).astype(np.uint8)
    full_singleton = np.full((8, 128, PROJECTION_ROWS), 255, dtype=np.uint8)
    for direction, eps in enumerate(epsilons):
        block = full[:, 2 + 35 * direction : 2 + 35 * (direction + 1)].astype(
            np.int64
        )
        for mask in range(128):
            if mask.bit_count() not in (0, 2, 4, 6):
                continue
            odd = {value for value in range(7) if mask & (1 << value)}
            values = mapped_catalog(
                mask.bit_count(), int(eps == -1), floor_for(eps, mask.bit_count()), odd, None
            ).astype(np.int64)
            if len(values) == 1:
                full_singleton[direction, mask] = (
                    block @ (13 - values[0]) % MODULUS
                ).astype(np.uint8)

    isolate_base = np.zeros((len(contexts), PROJECTION_ROWS), dtype=np.uint8)
    isolate_singleton = np.full(
        (len(contexts), 8, 128, PROJECTION_ROWS), 255, dtype=np.uint8
    )
    for context, (variable_mask, _isolated_direction) in enumerate(contexts):
        projection = isolates[context].astype(np.int64)
        isolate_base[context] = (
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
                    mask.bit_count(),
                    int(eps == -1),
                    floor_for(eps, mask.bit_count()),
                    odd,
                    None,
                ).astype(np.int64)
                if len(values) == 1:
                    isolate_singleton[context, direction, mask] = (
                        block @ (13 - values[0]) % MODULUS
                    ).astype(np.uint8)
    return {
        "full_base": full_base,
        "full_singleton": full_singleton,
        "isolate_base": isolate_base,
        "isolate_singleton": isolate_singleton,
    }


def build_catalog_tables(
    triples: np.ndarray,
    masks: np.ndarray,
    variable_supports: np.ndarray,
    contexts: list[tuple[int, int]],
    context_index: np.ndarray,
    full: np.ndarray,
    isolates: np.ndarray,
    epsilons: tuple[int, ...],
    allocation_supports: np.ndarray,
    allocation_increments: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict]:
    active: set[tuple[int, int, int, int]] = set()
    for candidate, (_rank, stratum_value, leaf_value) in enumerate(triples):
        stratum = int(stratum_value)
        leaf = int(leaf_value)
        raised = int(allocation_supports[stratum, leaf])
        variable = int(variable_supports[candidate])
        for direction in range(8):
            if not variable & (1 << direction):
                continue
            context = int(context_index[variable, direction])
            increment = (
                int(allocation_increments[stratum, leaf, direction])
                if raised & (1 << direction)
                else 0
            )
            active.add(
                (
                    context,
                    direction,
                    increment_index(increment),
                    int(masks[candidate, direction]),
                )
            )

    offsets = np.zeros((len(contexts), 8, 3, 128), dtype=np.uint64)
    counts = np.zeros((len(contexts), 8, 3, 128), dtype=np.uint32)
    isolate_pieces = []
    full_pieces = []
    cursor = 0
    raw_histogram: Counter[int] = Counter()
    filtered_pair_histogram: Counter[int] = Counter()
    started = time.time()
    for number, (context, direction, inc_index, mask) in enumerate(sorted(active), 1):
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
        bad = (13 - values).T
        isolate_block = isolates[
            context, :, 2 + 35 * direction : 2 + 35 * (direction + 1)
        ].astype(np.int64)
        full_block = full[
            :, 2 + 35 * direction : 2 + 35 * (direction + 1)
        ].astype(np.int64)
        isolate_signatures = pack_rows(
            (isolate_block @ bad % MODULUS).T.astype(np.uint8)
        )
        full_signatures = pack_rows(
            (full_block @ bad % MODULUS).T.astype(np.uint8)
        )
        pairs = np.column_stack((isolate_signatures, full_signatures))
        order = np.lexsort((pairs[:, 1], pairs[:, 0]))
        pairs = pairs[order]
        if len(pairs) > 1:
            keep = np.ones(len(pairs), dtype=bool)
            keep[1:] = np.any(pairs[1:] != pairs[:-1], axis=1)
            pairs = pairs[keep]
        offsets[context, direction, inc_index, mask] = cursor
        counts[context, direction, inc_index, mask] = len(pairs)
        isolate_pieces.append(pairs[:, 0].astype(np.uint64, copy=False))
        full_pieces.append(pairs[:, 1].astype(np.uint64, copy=False))
        cursor += len(pairs)
        raw_histogram[len(values)] += 1
        filtered_pair_histogram[len(pairs)] += 1
        if number % 200 == 0 or number == len(active):
            print(
                json.dumps(
                    {
                        "catalog_progress": [number, len(active)],
                        "signature_pairs": cursor,
                        "elapsed_seconds": time.time() - started,
                    }
                ),
                flush=True,
            )
    isolate_keys = np.concatenate(isolate_pieces).astype(np.uint64, copy=False)
    full_signatures = np.concatenate(full_pieces).astype(np.uint64, copy=False)
    arrays = {
        "offsets": offsets,
        "counts": counts,
        "isolate_keys": isolate_keys,
        "full_signatures": full_signatures,
    }
    return arrays, {
        "active_catalog_keys": len(active),
        "signature_pair_count": int(cursor),
        "raw_catalog_size_histogram": dict(sorted(raw_histogram.items())),
        "distinct_signature_pair_size_histogram": dict(
            sorted(filtered_pair_histogram.items())
        ),
        "array_shapes": {key: list(value.shape) for key, value in arrays.items()},
        "array_sha256": {key: array_sha256(value) for key, value in arrays.items()},
        "elapsed_seconds": time.time() - started,
    }


CUDA_SOURCE = r'''
extern "C" {

#define HASH_CAPACITY 4096
#define HASH_MASK 4095
#define MAX_HASH_BUILD 1800ULL
#define EMPTY_KEY 0xffffffffffffffffULL

__device__ __forceinline__ unsigned long long mix64(unsigned long long x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

__device__ __forceinline__ unsigned long long pack_negative_syndrome(
    const int* syndrome)
{
    unsigned long long output = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int row = 0; row < 22; ++row) {
        const unsigned int digit = (unsigned int)((7 - syndrome[row] % 7) % 7);
        output += place * (unsigned long long)digit;
        place *= 7ULL;
    }
    return output;
}

__device__ __forceinline__ unsigned long long sum_combination(
    unsigned long long combination,
    int group_mask,
    int variable_count,
    const unsigned long long* range_starts,
    const unsigned int* range_counts,
    const unsigned long long* full_signatures)
{
    unsigned long long values[5] = {0ULL, 0ULL, 0ULL, 0ULL, 0ULL};
    int selected = 0;
    #pragma unroll
    for (int variable = 0; variable < 5; ++variable) {
        if (variable >= variable_count || !(group_mask & (1 << variable))) continue;
        const unsigned int count = range_counts[variable];
        const unsigned int choice = (unsigned int)(combination % count);
        combination /= count;
        values[selected++] = full_signatures[range_starts[variable] + choice];
    }
    unsigned long long output = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int row = 0; row < 22; ++row) {
        unsigned int digit = 0U;
        #pragma unroll
        for (int index = 0; index < 5; ++index) {
            if (index < selected) {
                digit += (unsigned int)(values[index] % 7ULL);
                values[index] /= 7ULL;
            }
        }
        output += place * (unsigned long long)(digit % 7U);
        place *= 7ULL;
    }
    return output;
}

__device__ __forceinline__ unsigned long long needed_key(
    unsigned long long target, unsigned long long value)
{
    unsigned long long output = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int row = 0; row < 22; ++row) {
        const unsigned int digit = (unsigned int)(
            (target % 7ULL + 14ULL - value % 7ULL) % 7ULL
        );
        output += place * (unsigned long long)digit;
        target /= 7ULL;
        value /= 7ULL;
        place *= 7ULL;
    }
    return output;
}

__device__ __forceinline__ unsigned long long lower_bound_key(
    const unsigned long long* keys,
    unsigned long long offset,
    unsigned int count,
    unsigned long long target)
{
    unsigned int low = 0U;
    unsigned int high = count;
    while (low < high) {
        const unsigned int middle = low + ((high - low) >> 1);
        if (keys[offset + middle] < target) low = middle + 1U;
        else high = middle;
    }
    return offset + low;
}

__device__ __forceinline__ unsigned long long upper_bound_key(
    const unsigned long long* keys,
    unsigned long long offset,
    unsigned int count,
    unsigned long long target)
{
    unsigned int low = 0U;
    unsigned int high = count;
    while (low < high) {
        const unsigned int middle = low + ((high - low) >> 1);
        if (keys[offset + middle] <= target) low = middle + 1U;
        else high = middle;
    }
    return offset + low;
}

__device__ __forceinline__ void hash_insert(
    unsigned long long* table, unsigned long long key)
{
    unsigned int slot = (unsigned int)(mix64(key) & HASH_MASK);
    #pragma unroll 1
    for (int step = 0; step < HASH_CAPACITY; ++step) {
        const unsigned long long old = atomicCAS(table + slot, EMPTY_KEY, key);
        if (old == EMPTY_KEY || old == key) return;
        slot = (slot + 1U) & HASH_MASK;
    }
}

__device__ __forceinline__ bool hash_contains(
    const unsigned long long* table, unsigned long long key)
{
    unsigned int slot = (unsigned int)(mix64(key) & HASH_MASK);
    #pragma unroll 1
    for (int step = 0; step < HASH_CAPACITY; ++step) {
        const unsigned long long value = table[slot];
        if (value == key) return true;
        if (value == EMPTY_KEY) return false;
        slot = (slot + 1U) & HASH_MASK;
    }
    return false;
}

__global__ void full_catalog_join(
    const unsigned char* masks,
    const unsigned char* variable_supports,
    const unsigned char* strata,
    const unsigned char* leaves,
    const short* context_index,
    const unsigned char* allocation_supports,
    const unsigned char* allocation_increments,
    const unsigned char* full_base,
    const unsigned char* full_singleton,
    const unsigned char* isolate_base,
    const unsigned char* isolate_singleton,
    const unsigned long long* offsets,
    const unsigned int* counts,
    const unsigned long long* isolate_keys,
    const unsigned long long* full_signatures,
    const unsigned long long candidate_count,
    unsigned long long* survivor_indices,
    unsigned long long* survivor_count,
    unsigned long long* rejection_counts,
    unsigned long long* workload_maxima)
{
    const unsigned long long candidate = (unsigned long long)blockIdx.x;
    if (candidate >= candidate_count) return;
    __shared__ unsigned long long hash_table[HASH_CAPACITY];
    __shared__ unsigned long long range_starts[5];
    __shared__ unsigned int range_counts[5];
    __shared__ int variable_directions[5];
    __shared__ int variable_count;
    __shared__ int valid;
    __shared__ int capacity_failure;
    __shared__ int build_mask;
    __shared__ unsigned long long build_product;
    __shared__ unsigned long long probe_product;
    __shared__ unsigned long long global_target;
    __shared__ int match;

    if (threadIdx.x == 0) {
        valid = 1;
        capacity_failure = 0;
        match = 0;
        variable_count = 0;
        const int variable_mask = (int)variable_supports[candidate];
        const int stratum = (int)strata[candidate];
        const int leaf = (int)leaves[candidate];
        const int raised_mask = (int)allocation_supports[stratum * 44 + leaf];
        int global_syndrome[22];
        #pragma unroll
        for (int row = 0; row < 22; ++row)
            global_syndrome[row] = (int)full_base[row];
        #pragma unroll
        for (int direction = 0; direction < 8; ++direction) {
            const int mask = (int)masks[candidate * 8 + direction];
            if (!(variable_mask & (1 << direction))) {
                const unsigned long long start =
                    ((unsigned long long)direction * 128ULL + mask) * 22ULL;
                #pragma unroll
                for (int row = 0; row < 22; ++row) {
                    const int value = (int)full_singleton[start + row];
                    if (value == 255) valid = 0;
                    global_syndrome[row] += value;
                }
                continue;
            }
            const int local = variable_count++;
            variable_directions[local] = direction;
            const int context = (int)context_index[variable_mask * 8 + direction];
            if (context < 0) {
                valid = 0;
                continue;
            }
            int syndrome[22];
            #pragma unroll
            for (int row = 0; row < 22; ++row)
                syndrome[row] = (int)isolate_base[context * 22 + row];
            #pragma unroll
            for (int fixed_direction = 0; fixed_direction < 8; ++fixed_direction) {
                if (variable_mask & (1 << fixed_direction)) continue;
                const int fixed_mask =
                    (int)masks[candidate * 8 + fixed_direction];
                const unsigned long long start =
                    ((unsigned long long)context * 8ULL + fixed_direction)
                    * 128ULL * 22ULL
                    + (unsigned long long)fixed_mask * 22ULL;
                #pragma unroll
                for (int row = 0; row < 22; ++row) {
                    const int value = (int)isolate_singleton[start + row];
                    if (value == 255) valid = 0;
                    syndrome[row] += value;
                }
            }
            const unsigned long long target = pack_negative_syndrome(syndrome);
            const int increment = raised_mask & (1 << direction)
                ? (int)allocation_increments[((stratum * 44 + leaf) * 8) + direction]
                : 0;
            const int inc_index = increment == 0 ? 0 : (increment == 4 ? 1 : 2);
            const unsigned long long catalog_key =
                (((unsigned long long)context * 8ULL + direction) * 3ULL + inc_index)
                * 128ULL + mask;
            const unsigned long long offset = offsets[catalog_key];
            const unsigned int count = counts[catalog_key];
            const unsigned long long first = lower_bound_key(
                isolate_keys, offset, count, target
            );
            const unsigned long long stop = upper_bound_key(
                isolate_keys, offset, count, target
            );
            range_starts[local] = first;
            range_counts[local] = (unsigned int)(stop - first);
            if (first == stop) valid = 0;
        }
        global_target = pack_negative_syndrome(global_syndrome);
        build_mask = 0;
        build_product = 0ULL;
        probe_product = 0ULL;
        unsigned long long best_work = 0xffffffffffffffffULL;
        const int all_mask = (1 << variable_count) - 1;
        for (int subset = 1; subset < all_mask; ++subset) {
            if (!(subset & 1)) continue;
            unsigned long long left = 1ULL;
            unsigned long long right = 1ULL;
            #pragma unroll
            for (int variable = 0; variable < 5; ++variable) {
                if (variable >= variable_count) continue;
                if (subset & (1 << variable)) left *= range_counts[variable];
                else right *= range_counts[variable];
            }
            if (left <= MAX_HASH_BUILD && left + right < best_work) {
                best_work = left + right;
                build_mask = subset;
                build_product = left;
                probe_product = right;
            }
            if (right <= MAX_HASH_BUILD && left + right < best_work) {
                best_work = left + right;
                build_mask = all_mask ^ subset;
                build_product = right;
                probe_product = left;
            }
        }
        if (valid && (build_product == 0ULL || build_product > MAX_HASH_BUILD)) {
            valid = 0;
            capacity_failure = 1;
        }
        atomicMax(workload_maxima, build_product);
        atomicMax(workload_maxima + 1, probe_product);
    }
    __syncthreads();
    if (!valid) {
        if (threadIdx.x == 0)
            atomicAdd(rejection_counts + (capacity_failure ? 1 : 0), 1ULL);
        return;
    }
    for (int slot = threadIdx.x; slot < HASH_CAPACITY; slot += blockDim.x)
        hash_table[slot] = EMPTY_KEY;
    __syncthreads();
    for (unsigned long long combination = threadIdx.x;
         combination < build_product;
         combination += blockDim.x)
    {
        const unsigned long long key = sum_combination(
            combination,
            build_mask,
            variable_count,
            range_starts,
            range_counts,
            full_signatures
        );
        hash_insert(hash_table, key);
    }
    __syncthreads();
    const int all_mask = (1 << variable_count) - 1;
    const int probe_mask = all_mask ^ build_mask;
    for (unsigned long long combination = threadIdx.x;
         combination < probe_product && !match;
         combination += blockDim.x)
    {
        const unsigned long long value = sum_combination(
            combination,
            probe_mask,
            variable_count,
            range_starts,
            range_counts,
            full_signatures
        );
        if (hash_contains(hash_table, needed_key(global_target, value)))
            atomicExch(&match, 1);
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        if (match) {
            const unsigned long long slot = atomicAdd(survivor_count, 1ULL);
            survivor_indices[slot] = candidate;
        } else {
            atomicAdd(rejection_counts + 2, 1ULL);
        }
    }
}

}
'''


def sum_signatures(values: tuple[int, ...]) -> int:
    work = list(values)
    output = 0
    place = 1
    for _row in range(PROJECTION_ROWS):
        output += place * (sum(value % MODULUS for value in work) % MODULUS)
        work = [value // MODULUS for value in work]
        place *= MODULUS
    return output


def needed_signature(target: int, value: int) -> int:
    output = 0
    place = 1
    for _row in range(PROJECTION_ROWS):
        output += place * ((target % MODULUS - value % MODULUS) % MODULUS)
        target //= MODULUS
        value //= MODULUS
        place *= MODULUS
    return output


def choose_partition(counts: list[int]) -> tuple[int, int, int, int]:
    all_mask = (1 << len(counts)) - 1
    best = None
    for subset in range(1, all_mask):
        if not subset & 1:
            continue
        products = [
            math.prod(counts[index] for index in range(len(counts)) if mask & (1 << index))
            for mask in (subset, all_mask ^ subset)
        ]
        for build_slot in range(2):
            build = products[build_slot]
            probe = products[1 - build_slot]
            if build > MAX_HASH_BUILD:
                continue
            score = (build + probe, build, subset, build_slot)
            if best is None or score < best[0]:
                build_mask = subset if build_slot == 0 else all_mask ^ subset
                best = (score, build_mask, build, probe)
    if best is None:
        return 0, 0, 0, 0
    return best[1], all_mask ^ best[1], best[2], best[3]


def combination_sums(
    ranges: list[np.ndarray], group_mask: int
):
    indices = [index for index in range(len(ranges)) if group_mask & (1 << index)]
    for choices in __import__("itertools").product(*(ranges[index] for index in indices)):
        yield sum_signatures(tuple(int(value) for value in choices))


def cpu_candidate_passes(
    candidate: int,
    triples: np.ndarray,
    masks: np.ndarray,
    variable_supports: np.ndarray,
    context_index: np.ndarray,
    epsilons: tuple[int, ...],
    allocation_supports: np.ndarray,
    allocation_increments: np.ndarray,
    singleton: dict[str, np.ndarray],
    catalogs: dict[str, np.ndarray],
) -> tuple[bool, tuple[int, ...]]:
    stratum = int(triples[candidate, 1])
    leaf = int(triples[candidate, 2])
    variable = int(variable_supports[candidate])
    raised = int(allocation_supports[stratum, leaf])
    global_syndrome = singleton["full_base"].astype(np.int16).copy()
    for direction in range(8):
        if variable & (1 << direction):
            continue
        value = singleton["full_singleton"][direction, int(masks[candidate, direction])]
        if np.any(value == 255):
            return False, ()
        global_syndrome += value
    global_target = int(pack_rows(((-global_syndrome % MODULUS).astype(np.uint8))[None])[0])
    ranges = []
    for direction in range(8):
        if not variable & (1 << direction):
            continue
        context = int(context_index[variable, direction])
        syndrome = singleton["isolate_base"][context].astype(np.int16).copy()
        for fixed in range(8):
            if variable & (1 << fixed):
                continue
            value = singleton["isolate_singleton"][
                context, fixed, int(masks[candidate, fixed])
            ]
            if np.any(value == 255):
                return False, ()
            syndrome += value
        target = int(pack_rows(((-syndrome % MODULUS).astype(np.uint8))[None])[0])
        increment = (
            int(allocation_increments[stratum, leaf, direction])
            if raised & (1 << direction)
            else 0
        )
        inc_index = increment_index(increment)
        mask = int(masks[candidate, direction])
        offset = int(catalogs["offsets"][context, direction, inc_index, mask])
        count = int(catalogs["counts"][context, direction, inc_index, mask])
        keys = catalogs["isolate_keys"][offset : offset + count]
        first = int(np.searchsorted(keys, target, side="left"))
        stop = int(np.searchsorted(keys, target, side="right"))
        if first == stop:
            return False, ()
        ranges.append(catalogs["full_signatures"][offset + first : offset + stop])
    counts = [len(value) for value in ranges]
    build_mask, probe_mask, _build, _probe = choose_partition(counts)
    if not build_mask:
        return False, tuple(counts)
    lookup = set(combination_sums(ranges, build_mask))
    for value in combination_sums(ranges, probe_mask):
        if needed_signature(global_target, value) in lookup:
            return True, tuple(counts)
    return False, tuple(counts)


def gpu_scan(
    kernel,
    triples: np.ndarray,
    masks: np.ndarray,
    variable_supports: np.ndarray,
    context_index: np.ndarray,
    allocation_supports: np.ndarray,
    allocation_increments: np.ndarray,
    singleton: dict[str, np.ndarray],
    catalogs: dict[str, np.ndarray],
    threads: int,
) -> dict:
    import cupy as cp

    gpu = {
        "masks": cp.asarray(masks),
        "variable_supports": cp.asarray(variable_supports),
        "strata": cp.asarray(triples[:, 1].astype(np.uint8)),
        "leaves": cp.asarray(triples[:, 2].astype(np.uint8)),
        "context_index": cp.asarray(context_index),
        "allocation_supports": cp.asarray(allocation_supports),
        "allocation_increments": cp.asarray(allocation_increments),
        **{key: cp.asarray(value) for key, value in singleton.items()},
        **{key: cp.asarray(value) for key, value in catalogs.items()},
    }
    survivors = cp.zeros(len(triples), dtype=cp.uint64)
    survivor_count = cp.zeros(1, dtype=cp.uint64)
    rejection_counts = cp.zeros(3, dtype=cp.uint64)
    workload_maxima = cp.zeros(2, dtype=cp.uint64)
    started = time.perf_counter()
    kernel(
        (len(triples),),
        (threads,),
        (
            gpu["masks"],
            gpu["variable_supports"],
            gpu["strata"],
            gpu["leaves"],
            gpu["context_index"],
            gpu["allocation_supports"],
            gpu["allocation_increments"],
            gpu["full_base"],
            gpu["full_singleton"],
            gpu["isolate_base"],
            gpu["isolate_singleton"],
            gpu["offsets"],
            gpu["counts"],
            gpu["isolate_keys"],
            gpu["full_signatures"],
            np.uint64(len(triples)),
            survivors,
            survivor_count,
            rejection_counts,
            workload_maxima,
        ),
    )
    cp.cuda.runtime.deviceSynchronize()
    elapsed = time.perf_counter() - started
    total = int(cp.asnumpy(survivor_count)[0])
    return {
        "survivor_indices": np.sort(
            cp.asnumpy(survivors[:total]).astype(np.int64)
        ),
        "rejection_counts": [int(value) for value in cp.asnumpy(rejection_counts)],
        "workload_maxima": [int(value) for value in cp.asnumpy(workload_maxima)],
        "elapsed_seconds": elapsed,
    }


def run(args: argparse.Namespace) -> dict:
    import cupy as cp

    started = time.time()
    base, _mod7, _mod3 = load_intersection(args.mod7_source, args.mod3_source)
    triples, source = restrict_to_candidate_source(
        base, args.candidate_source, args.mod7_source, args.mod3_source
    )
    labels, epsilon_array = direction_tables()
    epsilons = tuple(int(value) for value in epsilon_array)
    _allocation_counts, allocation_supports, allocation_increments = allocation_patterns(
        epsilons
    )
    masks, variable_supports, contexts = prepare_candidates(
        triples, labels, epsilons, allocation_supports
    )
    matrix, dependencies, linear_rows = linear_data((MODULUS,))
    dependency = dependencies[MODULUS].astype(np.uint8)
    full, isolates, context_index, projection_records = build_projections(
        contexts, dependency, matrix
    )
    singleton = singleton_tables(full, isolates, contexts, epsilons)
    catalogs, catalog_summary = build_catalog_tables(
        triples,
        masks,
        variable_supports,
        contexts,
        context_index,
        full,
        isolates,
        epsilons,
        allocation_supports,
        allocation_increments,
    )
    kernel = cp.RawKernel(CUDA_SOURCE, "full_catalog_join", options=("--std=c++11",))
    prefix_count = min(args.cpu_prefix, len(triples))
    cpu_survivors = []
    prefix_allowed_histogram: Counter[tuple[int, ...]] = Counter()
    for candidate in range(prefix_count):
        passes, counts = cpu_candidate_passes(
            candidate,
            triples,
            masks,
            variable_supports,
            context_index,
            epsilons,
            allocation_supports,
            allocation_increments,
            singleton,
            catalogs,
        )
        prefix_allowed_histogram[counts] += 1
        if passes:
            cpu_survivors.append(candidate)
    prefix_gpu = gpu_scan(
        kernel,
        triples[:prefix_count],
        masks[:prefix_count],
        variable_supports[:prefix_count],
        context_index,
        allocation_supports,
        allocation_increments,
        singleton,
        catalogs,
        args.threads,
    )
    if cpu_survivors != prefix_gpu["survivor_indices"].tolist():
        raise AssertionError(
            f"CPU/GPU full-join prefix mismatch: {cpu_survivors[:20]} vs "
            f"{prefix_gpu['survivor_indices'][:20].tolist()}"
        )
    full_gpu = gpu_scan(
        kernel,
        triples,
        masks,
        variable_supports,
        context_index,
        allocation_supports,
        allocation_increments,
        singleton,
        catalogs,
        args.threads,
    )
    if full_gpu["rejection_counts"][1] != 0:
        raise AssertionError(
            "one or more candidates lacked an admissible exact hash partition"
        )
    if (
        sum(full_gpu["rejection_counts"])
        + len(full_gpu["survivor_indices"])
        != len(triples)
    ):
        raise AssertionError("full-join outcome accounting is incomplete")
    survivors = triples[full_gpu["survivor_indices"]]
    out = {
        "experiment": "p7_size8_full_catalog_filtered_gpu",
        "status": "complete_exact_filtered_full_catalog_join",
        "p": 7,
        "c_H": -1,
        "modulus": MODULUS,
        "projection_rows": PROJECTION_ROWS,
        "candidate_source": str(args.candidate_source),
        "candidate_source_sha256": sha256(args.candidate_source),
        "candidate_source_count": int(source["subset_survivor_count"]),
        "mod7_source_sha256": sha256(args.mod7_source),
        "mod3_source_sha256": sha256(args.mod3_source),
        "input_candidate_count": len(triples),
        "mathematical_reduction": {
            "single_catalog_filter": "dependencies vanish on every other variable direction",
            "full_join": "all variable catalogs joined under one common exact 22-row mod-7 projection",
            "hash_table": "lossless base-seven uint64 keys with exact open addressing",
            "logical_implication": "zero projected full join implies no exact catalog tuple",
        },
        "linear_system": linear_rows,
        "dependency_sha256": array_sha256(dependency),
        "full_projection_sha256": array_sha256(full),
        "isolate_projections_sha256": array_sha256(isolates),
        "projection_records": projection_records,
        "singleton_array_sha256": {
            key: array_sha256(value) for key, value in singleton.items()
        },
        "catalog_tables": catalog_summary,
        "verification": {
            "cpu_prefix_candidates": prefix_count,
            "cpu_prefix_survivors": len(cpu_survivors),
            "gpu_prefix_survivors": len(prefix_gpu["survivor_indices"]),
            "cpu_gpu_prefix_exact_match": True,
            "prefix_allowed_count_histogram": {
                ",".join(map(str, key)): value
                for key, value in sorted(prefix_allowed_histogram.items())
            },
            "left_null_audits": True,
            "isolate_omitted_block_audits": True,
            "rank_audits": True,
        },
        "device": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        "launch": {
            "blocks": len(triples),
            "threads": args.threads,
            "hash_capacity_per_block": HASH_CAPACITY,
            "maximum_hash_build": MAX_HASH_BUILD,
        },
        "maximum_observed_build_probe_products": full_gpu["workload_maxima"],
        "rejected_by_empty_single_filter": full_gpu["rejection_counts"][0],
        "rejected_by_hash_partition_capacity": full_gpu["rejection_counts"][1],
        "rejected_by_empty_single_filter_or_capacity": (
            full_gpu["rejection_counts"][0] + full_gpu["rejection_counts"][1]
        ),
        "rejected_by_full_join": full_gpu["rejection_counts"][2],
        "gpu_seconds": full_gpu["elapsed_seconds"],
        "full_join_survivor_count": len(survivors),
        "full_join_survivors_by_stratum": [
            int(np.count_nonzero(survivors[:, 1] == stratum)) for stratum in range(4)
        ],
        "full_join_survivor_rank_stratum_leaf": [
            list(map(int, row)) for row in survivors
        ],
        "all_input_candidates_excluded": len(survivors) == 0,
        "closes_cminus1_post_15664_scope": len(survivors) == 0,
        "elapsed_seconds": time.time() - started,
    }
    atomic_json(args.output, out)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mod7-source", type=Path, required=True)
    parser.add_argument("--mod3-source", type=Path, required=True)
    parser.add_argument("--candidate-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu-prefix", type=int, default=64)
    parser.add_argument("--threads", type=int, default=128)
    args = parser.parse_args()
    output = run(args)
    print(
        json.dumps(
            {
                "status": output["status"],
                "input_candidate_count": output["input_candidate_count"],
                "full_join_survivor_count": output["full_join_survivor_count"],
                "full_join_survivors_by_stratum": output[
                    "full_join_survivors_by_stratum"
                ],
                "all_input_candidates_excluded": output[
                    "all_input_candidates_excluded"
                ],
                "gpu_seconds": output["gpu_seconds"],
                "elapsed_seconds": output["elapsed_seconds"],
                "output": str(args.output),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
