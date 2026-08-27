#!/usr/bin/env python3
"""Complete projected V100 sieve for the post-15.664 p=7 size-eight scope.

Every remaining mean allocation raises a known support of one to five
direction catalogs.  Exact mod-seven dependencies conditioned to vanish on
that whole support reduce each leaf to singleton floor catalogs and at most
one 36-row floor catalog.  This script checks all C(49,8) boundaries and all
23,892,792 remaining allocation leaves, with an independent CPU prefix.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from p7_exceptional_omit_high_catalogs import modular_rank  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_size8_floor_profile_gpu import choose_table, direction_tables  # noqa: E402
from p7_size8_one_elevation_gpu import masks_for_boundary  # noqa: E402


PROJECTION_ROWS = 40
MAX_LEAVES = 44
STRATUM_FLOORS = ((32, 16), (24, 24), (32, 8), (24, 16))
EXPECTED_BOUNDARIES = (154_056, 1_194_816, 1_176, 69_384)
EXPECTED_ALLOCATION_COUNTS = (11, 16, 24, 44)
EXPECTED_LEAVES = tuple(
    boundaries * allocations
    for boundaries, allocations in zip(EXPECTED_BOUNDARIES, EXPECTED_ALLOCATION_COUNTS)
)


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


def allocation_patterns(
    epsilons: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    negative = tuple(index for index, eps in enumerate(epsilons) if eps == -1)
    positive = tuple(index for index, eps in enumerate(epsilons) if eps == 1)
    if len(negative) != 4 or len(positive) != 4:
        raise AssertionError("quadratic direction types changed")

    positive_16 = []
    for index in positive:
        row = [0] * 8
        row[index] = 16
        positive_16.append(row)
    for first, second in itertools.combinations(positive, 2):
        row = [0] * 8
        row[first] = row[second] = 8
        positive_16.append(row)
    row = [0] * 8
    for index in positive:
        row[index] = 4
    positive_16.append(row)
    if len(positive_16) != 11:
        raise AssertionError("positive budget-16 pattern count changed")

    patterns: list[list[list[int]]] = [positive_16]
    mixed_8_8 = []
    for negative_index in negative:
        for positive_index in positive:
            row = [0] * 8
            row[negative_index] = row[positive_index] = 8
            mixed_8_8.append(row)
    patterns.append(mixed_8_8)

    positive_24 = []
    for index in positive:
        row = [0] * 8
        row[index] = 24
        positive_24.append(row)
    for first in positive:
        for second in positive:
            if first == second:
                continue
            row = [0] * 8
            row[first] = 16
            row[second] = 8
            positive_24.append(row)
    for subset in itertools.combinations(positive, 3):
        row = [0] * 8
        for index in subset:
            row[index] = 8
        positive_24.append(row)
    for major in positive:
        row = [0] * 8
        for index in positive:
            row[index] = 12 if index == major else 4
        positive_24.append(row)
    if len(positive_24) != 24:
        raise AssertionError("positive budget-24 pattern count changed")
    patterns.append(positive_24)

    mixed_8_16 = []
    for negative_index in negative:
        for positive_row in positive_16:
            row = list(positive_row)
            row[negative_index] = 8
            mixed_8_16.append(row)
    patterns.append(mixed_8_16)

    counts = np.asarray([len(rows) for rows in patterns], dtype=np.uint8)
    if tuple(int(value) for value in counts) != EXPECTED_ALLOCATION_COUNTS:
        raise AssertionError("allocation counts changed")
    increments = np.zeros((4, MAX_LEAVES, 8), dtype=np.uint8)
    supports = np.zeros((4, MAX_LEAVES), dtype=np.uint8)
    for stratum, rows in enumerate(patterns):
        for leaf, row_values in enumerate(rows):
            increments[stratum, leaf] = row_values
            supports[stratum, leaf] = sum(
                1 << index for index, value in enumerate(row_values) if value
            )
    return counts, supports, increments


def fast_type_costs(
    profile: tuple[int, ...], epsilons: tuple[int, ...]
) -> tuple[int, int]:
    costs = []
    for selected_eps in (-1, 1):
        total = 0
        for b, eps in zip(profile, epsilons):
            if eps != selected_eps:
                continue
            total += (
                (0 if b == 0 else 8)
                if eps > 0
                else (14 if b in (0, 4) else 6)
            )
        costs.append(total)
    return tuple(costs)


def load_tables(
    table_path: Path,
    summary_path: Path,
    expected_labels: np.ndarray,
    expected_epsilons: np.ndarray,
    modulus: int,
) -> tuple[dict[str, np.ndarray], dict]:
    summary = json.loads(summary_path.read_text())
    if (
        summary.get("experiment") != "p7_size8_multi_elevation_tables"
        or summary.get("status") != "complete_exact_post_15664_omission_tables"
        or int(summary.get("projection_rows", 0)) != PROJECTION_ROWS
        or int(summary.get("modulus", 0)) != modulus
        or summary.get("output_sha256") != sha256(table_path)
    ):
        raise ValueError("multi-elevation table identity failed")
    with np.load(table_path, allow_pickle=False) as handle:
        arrays = {key: handle[key] for key in handle.files}
    if not np.array_equal(arrays["labels"], expected_labels):
        raise AssertionError("direction labels changed")
    if not np.array_equal(arrays["epsilons"], expected_epsilons):
        raise AssertionError("direction signs changed")
    for key, value in arrays.items():
        if summary["array_sha256"].get(key) != array_sha256(value):
            raise AssertionError(f"table hash changed for {key}")
    dependency = arrays["dependency"].astype(np.int64)
    rebuilt = (
        arrays["selected_coefficients"].astype(np.int64) @ dependency % modulus
    ).astype(np.uint8)
    if not np.array_equal(rebuilt, arrays["projected_dependencies"]):
        raise AssertionError("conditioned projections do not rebuild")
    matrix, dependencies, _rows = linear_data((modulus,))
    if not np.array_equal(arrays["dependency"], dependencies[modulus]):
        raise AssertionError("full dependency basis changed")
    for table_index, support_mask in enumerate(arrays["support_masks"]):
        projected = arrays["projected_dependencies"][table_index]
        if modular_rank(projected, modulus) != PROJECTION_ROWS:
            raise AssertionError("conditioned projection lost rank")
        if np.any(
            projected.astype(np.int64)
            @ (matrix.astype(np.int64) % modulus)
            % modulus
        ):
            raise AssertionError("conditioned projection is not left-null")
        for direction in range(8):
            if support_mask & (1 << direction):
                block = projected[:, 2 + 35 * direction : 2 + 35 * (direction + 1)]
                if np.any(block):
                    raise AssertionError("raised direction was not omitted")
    return arrays, summary


CUDA_SOURCE = r'''
extern "C" __global__
void remaining_scan(
    const signed char* labels,
    const signed char* epsilons,
    const unsigned long long* choose,
    const int choose_stride,
    const short* support_index,
    const unsigned char* allocation_counts,
    const unsigned char* allocation_supports,
    const unsigned char* base,
    const unsigned char* singleton,
    const unsigned char* variable,
    const unsigned char* variable_count,
    const int modulus,
    const unsigned long long start_rank,
    const unsigned long long stop_rank,
    unsigned long long* checked,
    unsigned long long* stratum_counts,
    unsigned long long* tested_leaves,
    unsigned long long* projected_leaves,
    unsigned long long* projected_boundaries,
    unsigned long long* odd_histogram,
    unsigned long long* survivor_ranks,
    unsigned char* survivor_strata,
    unsigned char* survivor_leaf_codes,
    const unsigned long long survivor_capacity)
{
    const unsigned long long logical_thread =
        (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned long long logical_threads =
        (unsigned long long)gridDim.x * blockDim.x;
    unsigned long long local_checked = 0;
    unsigned long long local_strata[4] = {0, 0, 0, 0};

    for (unsigned long long rank0 = start_rank + logical_thread;
         rank0 < stop_rank;
         rank0 += logical_threads)
    {
        int index[8];
        unsigned long long rank = rank0;
        int next = 0;
        bool valid = true;
        #pragma unroll
        for (int position = 0; position < 8; ++position)
        {
            const int remaining = 7 - position;
            bool selected = false;
            const int last = 49 - (remaining + 1);
            for (int candidate = next; candidate <= last; ++candidate)
            {
                const unsigned long long ways =
                    choose[(49 - candidate - 1) * choose_stride + remaining];
                if (rank < ways)
                {
                    index[position] = candidate;
                    next = candidate + 1;
                    selected = true;
                    break;
                }
                rank -= ways;
            }
            if (!selected) valid = false;
        }
        if (!valid) continue;
        ++local_checked;

        int masks[8];
        int type_minus = 0;
        int type_plus = 0;
        int odd_secants = 0;
        #pragma unroll
        for (int direction = 0; direction < 8; ++direction)
        {
            int mask = 0;
            #pragma unroll
            for (int column = 0; column < 8; ++column)
                mask ^= 1 << (int)labels[49 * direction + index[column]];
            masks[direction] = mask;
            const int b = __popc((unsigned int)mask);
            odd_secants += b;
            const int eps = (int)epsilons[direction];
            const int cost = eps > 0
                ? (b == 0 ? 0 : 8)
                : ((b == 0 || b == 4) ? 14 : 6);
            if (eps < 0) type_minus += cost;
            else type_plus += cost;
        }
        if (odd_secants == 8) continue;
        int stratum = -1;
        if (type_minus == 32 && type_plus == 16) stratum = 0;
        else if (type_minus == 24 && type_plus == 24) stratum = 1;
        else if (type_minus == 32 && type_plus == 8) stratum = 2;
        else if (type_minus == 24 && type_plus == 16) stratum = 3;
        else continue;

        ++local_strata[stratum];
        atomicAdd(odd_histogram + 57 * stratum + odd_secants, 1ULL);
        const int leaf_count = (int)allocation_counts[stratum];
        atomicAdd(tested_leaves + stratum, (unsigned long long)leaf_count);
        bool boundary_pass = false;
        for (int leaf = 0; leaf < leaf_count; ++leaf)
        {
            const int support_mask =
                (int)allocation_supports[stratum * 44 + leaf];
            const int table_index = (int)support_index[support_mask];
            if (table_index < 0) continue;
            int syndrome[40];
            #pragma unroll
            for (int row = 0; row < 40; ++row)
                syndrome[row] = (int)base[40 * table_index + row];
            int variable_direction = -1;
            int variable_mask = -1;
            bool table_valid = true;
            #pragma unroll
            for (int direction = 0; direction < 8; ++direction)
            {
                if (support_mask & (1 << direction)) continue;
                const int mask = masks[direction];
                const int key = (table_index * 8 + direction) * 128 + mask;
                const int count = (int)variable_count[key];
                if (count)
                {
                    if (variable_direction >= 0) table_valid = false;
                    variable_direction = direction;
                    variable_mask = mask;
                }
                else
                {
                    const int offset = key * 40;
                    #pragma unroll
                    for (int row = 0; row < 40; ++row)
                    {
                        const int value = (int)singleton[offset + row];
                        if (value == 255) table_valid = false;
                        syndrome[row] += value;
                    }
                }
            }
            if (!table_valid) continue;
            bool passing = false;
            if (variable_direction < 0)
            {
                passing = true;
                #pragma unroll
                for (int row = 0; row < 40; ++row)
                    if (syndrome[row] % modulus != 0) passing = false;
            }
            else
            {
                const int key =
                    (table_index * 8 + variable_direction) * 128 + variable_mask;
                const int offset = key * 36 * 40;
                for (int catalog_row = 0; catalog_row < 36 && !passing; ++catalog_row)
                {
                    bool match = true;
                    #pragma unroll
                    for (int row = 0; row < 40; ++row)
                    {
                        const int value =
                            (int)variable[offset + catalog_row * 40 + row];
                        if ((syndrome[row] + value) % modulus != 0) match = false;
                    }
                    passing = match;
                }
            }
            if (passing)
            {
                boundary_pass = true;
                atomicAdd(projected_leaves + stratum, 1ULL);
                const unsigned long long slot = atomicAdd(projected_leaves + 4, 1ULL);
                if (slot < survivor_capacity)
                {
                    survivor_ranks[slot] = rank0;
                    survivor_strata[slot] = (unsigned char)stratum;
                    survivor_leaf_codes[slot] = (unsigned char)leaf;
                }
            }
        }
        if (boundary_pass) atomicAdd(projected_boundaries + stratum, 1ULL);
    }
    if (local_checked) atomicAdd(checked, local_checked);
    #pragma unroll
    for (int stratum = 0; stratum < 4; ++stratum)
        if (local_strata[stratum])
            atomicAdd(stratum_counts + stratum, local_strata[stratum]);
}
'''


def projected_passes(
    masks: tuple[int, ...],
    support_mask: int,
    arrays: dict[str, np.ndarray],
    modulus: int,
) -> bool:
    table_index = int(arrays["support_index"][support_mask])
    if table_index < 0:
        raise AssertionError(f"missing support table {support_mask}")
    syndrome = arrays["base"][table_index].astype(np.int16).copy()
    variable_key = None
    for direction, mask in enumerate(masks):
        if support_mask & (1 << direction):
            continue
        count = int(arrays["variable_count"][table_index, direction, mask])
        if count:
            if variable_key is not None:
                raise AssertionError("leaf has more than one variable floor catalog")
            variable_key = (direction, mask, count)
        else:
            value = arrays["singleton"][table_index, direction, mask]
            if np.any(value == 255):
                raise AssertionError("missing singleton floor contribution")
            syndrome += value
    if variable_key is None:
        return bool(np.all(syndrome % modulus == 0))
    direction, mask, count = variable_key
    values = arrays["variable"][table_index, direction, mask, :count]
    return bool(
        np.any(np.all((values.astype(np.int16) + syndrome) % modulus == 0, axis=1))
    )


def cpu_prefix(
    stop_rank: int,
    labels: np.ndarray,
    epsilons: tuple[int, ...],
    arrays: dict[str, np.ndarray],
    allocation_counts: np.ndarray,
    allocation_supports: np.ndarray,
    modulus: int,
) -> dict:
    floor_to_stratum = {value: index for index, value in enumerate(STRATUM_FLOORS)}
    stratum_counts = Counter()
    tested_leaves = Counter()
    passed_leaves = Counter()
    passed_boundaries = Counter()
    odd_histogram: dict[int, Counter] = {index: Counter() for index in range(4)}
    pairs = []
    for rank, boundary in enumerate(
        itertools.islice(itertools.combinations(range(49), 8), stop_rank)
    ):
        masks = masks_for_boundary(boundary, labels)
        profile = tuple(mask.bit_count() for mask in masks)
        if sum(profile) == 8:
            continue
        stratum = floor_to_stratum.get(fast_type_costs(profile, epsilons))
        if stratum is None:
            continue
        stratum_counts[stratum] += 1
        odd_histogram[stratum][sum(profile)] += 1
        boundary_pass = False
        for leaf in range(int(allocation_counts[stratum])):
            tested_leaves[stratum] += 1
            support = int(allocation_supports[stratum, leaf])
            if projected_passes(masks, support, arrays, modulus):
                passed_leaves[stratum] += 1
                pairs.append((rank, stratum, leaf))
                boundary_pass = True
        passed_boundaries[stratum] += int(boundary_pass)
    return {
        "checked": stop_rank,
        "stratum_counts": [stratum_counts[index] for index in range(4)],
        "tested_leaves": [tested_leaves[index] for index in range(4)],
        "projected_leaves": [passed_leaves[index] for index in range(4)],
        "projected_boundaries": [passed_boundaries[index] for index in range(4)],
        "odd_histogram": [dict(sorted(odd_histogram[index].items())) for index in range(4)],
        "survivor_triples": pairs,
    }


def gpu_scan(
    kernel,
    arrays_gpu: dict,
    choose_gpu,
    allocation_counts_gpu,
    allocation_supports_gpu,
    start_rank: int,
    stop_rank: int,
    blocks: int,
    threads: int,
    survivor_capacity: int,
    modulus: int,
) -> dict:
    import cupy as cp

    checked = cp.zeros(1, dtype=cp.uint64)
    stratum_counts = cp.zeros(4, dtype=cp.uint64)
    tested_leaves = cp.zeros(4, dtype=cp.uint64)
    projected_leaves = cp.zeros(5, dtype=cp.uint64)
    projected_boundaries = cp.zeros(4, dtype=cp.uint64)
    odd_histogram = cp.zeros((4, 57), dtype=cp.uint64)
    ranks = cp.zeros(survivor_capacity, dtype=cp.uint64)
    strata = cp.zeros(survivor_capacity, dtype=cp.uint8)
    leaves = cp.zeros(survivor_capacity, dtype=cp.uint8)
    started = time.perf_counter()
    kernel(
        (blocks,),
        (threads,),
        (
            arrays_gpu["labels"],
            arrays_gpu["epsilons"],
            choose_gpu,
            np.int32(choose_gpu.shape[1]),
            arrays_gpu["support_index"],
            allocation_counts_gpu,
            allocation_supports_gpu,
            arrays_gpu["base"],
            arrays_gpu["singleton"],
            arrays_gpu["variable"],
            arrays_gpu["variable_count"],
            np.int32(modulus),
            np.uint64(start_rank),
            np.uint64(stop_rank),
            checked,
            stratum_counts,
            tested_leaves,
            projected_leaves,
            projected_boundaries,
            odd_histogram,
            ranks,
            strata,
            leaves,
            np.uint64(survivor_capacity),
        ),
    )
    cp.cuda.runtime.deviceSynchronize()
    elapsed = time.perf_counter() - started
    total_passed = int(cp.asnumpy(projected_leaves)[4])
    if total_passed > survivor_capacity:
        raise RuntimeError(
            f"projection survivor capacity exceeded: {total_passed}>{survivor_capacity}"
        )
    rank_values = cp.asnumpy(ranks[:total_passed])
    stratum_values = cp.asnumpy(strata[:total_passed])
    leaf_values = cp.asnumpy(leaves[:total_passed])
    histogram_host = cp.asnumpy(odd_histogram)
    return {
        "checked": int(cp.asnumpy(checked)[0]),
        "stratum_counts": [int(value) for value in cp.asnumpy(stratum_counts)],
        "tested_leaves": [int(value) for value in cp.asnumpy(tested_leaves)],
        "projected_leaves": [int(value) for value in cp.asnumpy(projected_leaves)[:4]],
        "projected_boundaries": [int(value) for value in cp.asnumpy(projected_boundaries)],
        "odd_histogram": [
            {
                int(index): int(value)
                for index, value in enumerate(histogram_host[stratum])
                if value
            }
            for stratum in range(4)
        ],
        "survivor_triples": sorted(
            (int(rank), int(stratum), int(leaf))
            for rank, stratum, leaf in zip(rank_values, stratum_values, leaf_values)
        ),
        "elapsed_seconds": elapsed,
    }


def run(args: argparse.Namespace) -> dict:
    import cupy as cp

    started = time.time()
    structure = json.loads(args.structure.read_text())
    if (
        structure.get("experiment") != "p7_size8_remaining_allocation_structure"
        or structure.get("status") != "complete_exact_post_15664_structure"
        or int(structure.get("remaining_boundaries_per_sign", 0))
        != sum(EXPECTED_BOUNDARIES)
    ):
        raise ValueError("remaining-structure input failed identity checks")
    labels, epsilon_array = direction_tables()
    epsilons = tuple(int(value) for value in epsilon_array)
    allocation_counts, allocation_supports, allocation_increments = allocation_patterns(
        epsilons
    )
    arrays, table_summary = load_tables(
        args.tables, args.table_summary, labels, epsilon_array, args.modulus
    )
    for support in allocation_supports[allocation_supports != 0]:
        if int(arrays["support_index"][int(support)]) < 0:
            raise AssertionError(f"allocation support {support} has no table")

    if cp.cuda.runtime.getDeviceCount() < 1:
        raise RuntimeError("CUDA is unavailable")
    kernel = cp.RawKernel(CUDA_SOURCE, "remaining_scan", options=("--std=c++11",))
    arrays_gpu = {
        key: cp.asarray(arrays[key])
        for key in (
            "labels",
            "epsilons",
            "support_index",
            "base",
            "singleton",
            "variable",
            "variable_count",
        )
    }
    choose_gpu = cp.asarray(choose_table(49, 8))
    allocation_counts_gpu = cp.asarray(allocation_counts)
    allocation_supports_gpu = cp.asarray(allocation_supports)

    prefix_gpu = gpu_scan(
        kernel,
        arrays_gpu,
        choose_gpu,
        allocation_counts_gpu,
        allocation_supports_gpu,
        0,
        args.verify_prefix,
        args.blocks,
        args.threads,
        args.survivor_capacity,
        args.modulus,
    )
    prefix_cpu = cpu_prefix(
        args.verify_prefix,
        labels,
        epsilons,
        arrays,
        allocation_counts,
        allocation_supports,
        args.modulus,
    )
    prefix_gpu_comparable = {
        key: value for key, value in prefix_gpu.items() if key != "elapsed_seconds"
    }
    if prefix_gpu_comparable != prefix_cpu:
        raise AssertionError(
            f"independent CPU/GPU prefix mismatch: gpu={prefix_gpu_comparable}, "
            f"cpu={prefix_cpu}"
        )

    total = math.comb(49, 8)
    full = gpu_scan(
        kernel,
        arrays_gpu,
        choose_gpu,
        allocation_counts_gpu,
        allocation_supports_gpu,
        0,
        total,
        args.blocks,
        args.threads,
        args.survivor_capacity,
        args.modulus,
    )
    if full["checked"] != total:
        raise AssertionError("full CUDA rank interval was not exhausted")
    if tuple(full["stratum_counts"]) != EXPECTED_BOUNDARIES:
        raise AssertionError(f"stratum boundary counts changed: {full['stratum_counts']}")
    if tuple(full["tested_leaves"]) != EXPECTED_LEAVES:
        raise AssertionError(f"tested leaf counts changed: {full['tested_leaves']}")

    properties = cp.cuda.runtime.getDeviceProperties(0)
    device = properties["name"]
    if isinstance(device, bytes):
        device = device.decode()
    survivors = full["survivor_triples"]
    return {
        "experiment": "p7_size8_remaining_gpu",
        "status": "complete_exact_post_15664_projected_exhaustion",
        "p": 7,
        "c_H": -1,
        "modulus": args.modulus,
        "finite_boundary_size": 8,
        "source_structure": str(args.structure),
        "source_structure_sha256": sha256(args.structure),
        "tables": str(args.tables),
        "tables_sha256": sha256(args.tables),
        "table_summary": str(args.table_summary),
        "table_summary_sha256": sha256(args.table_summary),
        "table_array_sha256": table_summary["array_sha256"],
        "mathematical_reduction": {
            "stratum_floor_pairs": [list(value) for value in STRATUM_FLOORS],
            "allocations_per_boundary": list(EXPECTED_ALLOCATION_COUNTS),
            "raised_support_sizes": [1, 2, 3, 4, 5],
            "valid_raised_supports": int(len(arrays["support_masks"])),
            "minimum_conditioned_dependency_dimension": int(
                arrays["conditioned_dimensions"].min()
            ),
            "selected_conditioned_dependencies": PROJECTION_ROWS,
            "raised_catalog_blocks_identically_zero": True,
            "maximum_remaining_variable_floor_catalogs": 1,
        },
        "allocation_pattern_sha256": array_sha256(allocation_increments),
        "verification": {
            "method": "independent itertools CPU prefix versus direct-rank CUDA",
            "prefix_checked": args.verify_prefix,
            "prefix_exact_match": True,
        },
        "device": str(device),
        "launch": {"blocks": args.blocks, "threads": args.threads},
        "all_boundaries": total,
        "checked_boundaries": full["checked"],
        "stratum_boundaries": full["stratum_counts"],
        "tested_allocation_leaves": full["tested_leaves"],
        "odd_secant_histograms": full["odd_histogram"],
        "projected_survivor_leaves": full["projected_leaves"],
        "projected_survivor_boundaries": full["projected_boundaries"],
        "projected_survivor_rank_stratum_leaf": [list(value) for value in survivors],
        "all_remaining_leaves_mod7_infeasible": not survivors,
        "closes_cminus1_post_15664_scope": not survivors,
        "closes_all_p7_size_eight": False,
        "closes_residual_ii": False,
        "elapsed_seconds": time.time() - started,
        "cuda_scan_seconds": full["elapsed_seconds"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--table-summary", type=Path, required=True)
    parser.add_argument("--blocks", type=int, default=4096)
    parser.add_argument("--threads", type=int, default=256)
    parser.add_argument("--verify-prefix", type=int, default=100_000)
    parser.add_argument("--survivor-capacity", type=int, default=1_000_000)
    parser.add_argument("--modulus", type=int, choices=(3, 7), default=7)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 0 < args.verify_prefix <= math.comb(49, 8):
        raise ValueError("verify-prefix is outside the rank interval")
    output = run(args)
    atomic_json(args.output, output)
    print(
        json.dumps(
            {
                "status": output["status"],
                "checked_boundaries": output["checked_boundaries"],
                "tested_allocation_leaves": output["tested_allocation_leaves"],
                "projected_survivor_leaves": output["projected_survivor_leaves"],
                "projected_survivor_boundaries": output[
                    "projected_survivor_boundaries"
                ],
                "all_remaining_leaves_mod7_infeasible": output[
                    "all_remaining_leaves_mod7_infeasible"
                ],
                "output": str(args.output),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
