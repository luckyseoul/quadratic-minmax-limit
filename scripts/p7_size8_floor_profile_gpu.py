#!/usr/bin/env python3
"""Exact direct-rank CUDA floor census for eight finite points at p=7.

This is a reconnaissance/audit tool for the first open boundary size after
Proposition 15.661.  It scans all C(49,8) finite boundaries without building
the combinations on the CPU.  For each boundary it records the eight odd
fibre counts, their total (the number of affine odd secants), and whether the
two exact quadratic-type floor budgets from Proposition 15.632 are met.

At p=7 and even boundary size the only possible odd-fibre counts are
0,2,4,6.  The exact scaled floors are

    phase zero: (0,8,8,8),   phase one: (14,6,14,6),

and each four-direction quadratic type has budget 32.  No conclusion beyond
this necessary floor condition is claimed by the script.
"""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import islice
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def direction_tables() -> tuple[np.ndarray, np.ndarray]:
    labels = []
    epsilons = []
    for direction in projective_directions(7):
        epsilon, row = field_direction_data(7, direction)
        epsilons.append(int(epsilon))
        labels.append(tuple(int(value) for value in row))
    label_array = np.asarray(labels, dtype=np.int8)
    epsilon_array = np.asarray(epsilons, dtype=np.int8)
    if label_array.shape != (8, 49):
        raise AssertionError("unexpected direction-label table")
    if sorted(epsilon_array.tolist()) != [-1] * 4 + [1] * 4:
        raise AssertionError("unexpected quadratic-type split")
    return label_array, epsilon_array


def choose_table(n: int, support: int) -> np.ndarray:
    table = np.zeros((n + 1, support + 1), dtype=np.uint64)
    for a in range(n + 1):
        for b in range(min(a, support) + 1):
            value = math.comb(a, b)
            if value >= 1 << 64:
                raise OverflowError(f"C({a},{b}) does not fit uint64")
            table[a, b] = value
    return table


def decode_profile(code: int) -> list[int]:
    row = [0] * 8
    for direction in range(7, -1, -1):
        row[direction] = 2 * (code & 3)
        code >>= 2
    return row


def independent_cpu_prefix(c_h: int, stop_rank: int) -> dict:
    """Independently enumerate a short lexicographic prefix for GPU auditing."""
    from itertools import combinations

    labels, epsilons = direction_tables()
    all_odd: Counter[int] = Counter()
    survivor_odd: Counter[int] = Counter()
    survivor_profiles: Counter[int] = Counter()
    survivor_b_counts: Counter[int] = Counter()
    survivors = 0
    checked = 0
    for index in islice(combinations(range(49), 8), stop_rank):
        b_values = []
        type_cost = {-1: 0, 1: 0}
        profile = 0
        for direction in range(8):
            mask = 0
            for point in index:
                mask ^= 1 << int(labels[direction, point])
            b = mask.bit_count()
            b_values.append(b)
            profile = (profile << 2) | (b >> 1)
            phase = int(epsilons[direction]) == c_h
            if not phase:
                cost = 0 if b == 0 else 8
            else:
                cost = 14 if b in (0, 4) else 6
            type_cost[int(epsilons[direction])] += cost
        odd_secants = sum(b_values)
        all_odd[odd_secants] += 1
        checked += 1
        if type_cost[-1] <= 32 and type_cost[1] <= 32:
            survivors += 1
            survivor_odd[odd_secants] += 1
            survivor_profiles[profile] += 1
            survivor_b_counts.update(b_values)
    return {
        "checked": checked,
        "survivors": survivors,
        "all_odd": dict(all_odd),
        "survivor_odd": dict(survivor_odd),
        "survivor_profiles": dict(survivor_profiles),
        "survivor_b_counts": dict(survivor_b_counts),
    }


def verify_cpu_prefix(result: dict) -> None:
    start, stop = result["rank_interval"]
    if start != 0:
        raise ValueError("independent CPU verification only supports a rank-zero prefix")
    expected = independent_cpu_prefix(result["c_H"], stop)
    observed_profiles = {
        sum((b >> 1) << (2 * (7 - direction)) for direction, b in enumerate(row["b_by_direction"])):
        row["count"]
        for row in result["survivor_ordered_profiles"]
    }
    comparisons = {
        "checked": (result["checked_boundaries"], expected["checked"]),
        "survivors": (result["floor_surviving_boundaries"], expected["survivors"]),
        "all_odd": (
            {int(key): value for key, value in result["all_odd_secant_histogram"].items()},
            expected["all_odd"],
        ),
        "survivor_odd": (
            {int(key): value for key, value in result["survivor_odd_secant_histogram"].items()},
            expected["survivor_odd"],
        ),
        "survivor_profiles": (observed_profiles, expected["survivor_profiles"]),
        "survivor_b_counts": (
            {int(key): value for key, value in result["survivor_odd_fibre_count_histogram"].items()},
            expected["survivor_b_counts"],
        ),
    }
    mismatches = [name for name, pair in comparisons.items() if pair[0] != pair[1]]
    if mismatches:
        raise AssertionError(f"GPU/CPU prefix mismatch: {mismatches}")
    result["verification"] = {
        "method": "independent_itertools_prefix_census",
        "checked_boundaries": expected["checked"],
        "all_histograms_and_ordered_profiles_match": True,
    }


KERNEL = r'''
extern "C" __global__
void size8_floor_profiles(
    const signed char* labels,
    const signed char* epsilons,
    const unsigned long long* choose,
    const int choose_stride,
    const int c_h,
    const unsigned long long start_rank,
    const unsigned long long stop_rank,
    unsigned long long* checked,
    unsigned long long* survivors,
    unsigned long long* all_odd_secants,
    unsigned long long* survivor_odd_secants,
    unsigned long long* survivor_profiles,
    unsigned long long* survivor_b_counts,
    unsigned long long* minimum_rank_count,
    unsigned long long* minimum_ranks,
    unsigned long long* survivor_minimum_rank_count,
    unsigned long long* survivor_minimum_ranks,
    const unsigned long long minimum_rank_capacity)
{
    const unsigned long long logical_thread =
        (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned long long logical_threads =
        (unsigned long long)gridDim.x * blockDim.x;

    unsigned long long local_checked = 0;
    unsigned long long local_survivors = 0;
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

        int type_cost_minus = 0;
        int type_cost_plus = 0;
        int odd_secants = 0;
        unsigned int profile = 0;
        int b_values[8];
        #pragma unroll
        for (int direction = 0; direction < 8; ++direction)
        {
            unsigned int mask = 0;
            #pragma unroll
            for (int column = 0; column < 8; ++column)
                mask ^= 1U << (unsigned int)labels[49 * direction + index[column]];
            const int b = __popc(mask);
            b_values[direction] = b;
            odd_secants += b;
            profile = (profile << 2) | (unsigned int)(b >> 1);

            const int phase = epsilons[direction] == c_h;
            int cost;
            if (!phase)
                cost = b == 0 ? 0 : 8;
            else
                cost = (b == 0 || b == 4) ? 14 : 6;
            if (epsilons[direction] < 0) type_cost_minus += cost;
            else type_cost_plus += cost;
        }

        atomicAdd(all_odd_secants + odd_secants, 1ULL);
        if (odd_secants == 8)
        {
            const unsigned long long slot = atomicAdd(minimum_rank_count, 1ULL);
            if (slot < minimum_rank_capacity) minimum_ranks[slot] = rank0;
        }
        ++local_checked;
        if (type_cost_minus <= 32 && type_cost_plus <= 32)
        {
            ++local_survivors;
            atomicAdd(survivor_odd_secants + odd_secants, 1ULL);
            atomicAdd(survivor_profiles + profile, 1ULL);
            #pragma unroll
            for (int direction = 0; direction < 8; ++direction)
                atomicAdd(survivor_b_counts + (b_values[direction] >> 1), 1ULL);
            if (odd_secants == 8)
            {
                const unsigned long long slot =
                    atomicAdd(survivor_minimum_rank_count, 1ULL);
                if (slot < minimum_rank_capacity) survivor_minimum_ranks[slot] = rank0;
            }
        }
    }
    if (local_checked) atomicAdd(checked, local_checked);
    if (local_survivors) atomicAdd(survivors, local_survivors);
}
'''


def run(c_h: int, blocks: int, threads: int, start_rank: int, stop_rank: int | None) -> dict:
    import cupy as cp

    labels, epsilons = direction_tables()
    floors = {
        phase: {b: scaled_direction_floor(7, b, phase) for b in (0, 2, 4, 6)}
        for phase in (0, 1)
    }
    expected_floors = {
        0: {0: 0, 2: 8, 4: 8, 6: 8},
        1: {0: 14, 2: 6, 4: 14, 6: 6},
    }
    if floors != expected_floors:
        raise AssertionError(f"unexpected floor table: {floors}")

    total = math.comb(49, 8)
    stop = total if stop_rank is None else stop_rank
    if not 0 <= start_rank <= stop <= total:
        raise ValueError("rank interval is outside C(49,8)")

    kernel = cp.RawKernel(KERNEL, "size8_floor_profiles", options=("--std=c++11",))
    labels_gpu = cp.asarray(labels)
    epsilons_gpu = cp.asarray(epsilons)
    choose_gpu = cp.asarray(choose_table(49, 8))
    checked_gpu = cp.zeros(1, dtype=cp.uint64)
    survivors_gpu = cp.zeros(1, dtype=cp.uint64)
    all_odd_gpu = cp.zeros(57, dtype=cp.uint64)
    survivor_odd_gpu = cp.zeros(57, dtype=cp.uint64)
    profiles_gpu = cp.zeros(1 << 16, dtype=cp.uint64)
    b_counts_gpu = cp.zeros(4, dtype=cp.uint64)
    minimum_rank_capacity = 8192
    minimum_rank_count_gpu = cp.zeros(1, dtype=cp.uint64)
    minimum_ranks_gpu = cp.zeros(minimum_rank_capacity, dtype=cp.uint64)
    survivor_minimum_rank_count_gpu = cp.zeros(1, dtype=cp.uint64)
    survivor_minimum_ranks_gpu = cp.zeros(minimum_rank_capacity, dtype=cp.uint64)

    started = time.perf_counter()
    kernel(
        (blocks,),
        (threads,),
        (
            labels_gpu,
            epsilons_gpu,
            choose_gpu,
            np.int32(choose_gpu.shape[1]),
            np.int32(c_h),
            np.uint64(start_rank),
            np.uint64(stop),
            checked_gpu,
            survivors_gpu,
            all_odd_gpu,
            survivor_odd_gpu,
            profiles_gpu,
            b_counts_gpu,
            minimum_rank_count_gpu,
            minimum_ranks_gpu,
            survivor_minimum_rank_count_gpu,
            survivor_minimum_ranks_gpu,
            np.uint64(minimum_rank_capacity),
        ),
    )
    cp.cuda.runtime.deviceSynchronize()
    elapsed = time.perf_counter() - started

    checked = int(cp.asnumpy(checked_gpu)[0])
    survivor_count = int(cp.asnumpy(survivors_gpu)[0])
    all_odd = cp.asnumpy(all_odd_gpu)
    survivor_odd = cp.asnumpy(survivor_odd_gpu)
    profiles = cp.asnumpy(profiles_gpu)
    b_counts = cp.asnumpy(b_counts_gpu)
    minimum_rank_count = int(cp.asnumpy(minimum_rank_count_gpu)[0])
    survivor_minimum_rank_count = int(cp.asnumpy(survivor_minimum_rank_count_gpu)[0])
    if minimum_rank_count > minimum_rank_capacity:
        raise AssertionError("minimum-rank output capacity was exceeded")
    if survivor_minimum_rank_count > minimum_rank_capacity:
        raise AssertionError("surviving minimum-rank output capacity was exceeded")
    minimum_ranks = sorted(
        int(value)
        for value in cp.asnumpy(minimum_ranks_gpu[:minimum_rank_count])
    )
    survivor_minimum_ranks = sorted(
        int(value)
        for value in cp.asnumpy(
            survivor_minimum_ranks_gpu[:survivor_minimum_rank_count]
        )
    )
    properties = cp.cuda.runtime.getDeviceProperties(0)
    device = properties["name"]
    if isinstance(device, bytes):
        device = device.decode()

    profile_rows = [
        {"b_by_direction": decode_profile(code), "count": int(count)}
        for code, count in enumerate(profiles)
        if count
    ]
    profile_rows.sort(key=lambda row: (-row["count"], row["b_by_direction"]))
    return {
        "experiment": "p7_size8_floor_profile_gpu",
        "status": "complete_exact_floor_profile_census" if checked == stop - start_rank else "incomplete",
        "p": 7,
        "c_H": c_h,
        "infinity_value": 0,
        "finite_boundary_size": 8,
        "all_boundaries": total,
        "rank_interval": [start_rank, stop],
        "checked_boundaries": checked,
        "floor_surviving_boundaries": survivor_count,
        "floor_rejected_boundaries": checked - survivor_count,
        "device": str(device),
        "launch": {"blocks": blocks, "threads": threads},
        "elapsed_seconds": elapsed,
        "boundaries_per_second": checked / elapsed,
        "criterion": {
            "phase_one_quadratic_type": c_h,
            "floors": {
                str(phase): {str(b): value for b, value in row.items()}
                for phase, row in floors.items()
            },
            "budget_per_quadratic_type": 32,
        },
        "all_odd_secant_histogram": {
            str(index): int(count) for index, count in enumerate(all_odd) if count
        },
        "survivor_odd_secant_histogram": {
            str(index): int(count) for index, count in enumerate(survivor_odd) if count
        },
        "survivor_odd_fibre_count_histogram": {
            str(2 * index): int(count) for index, count in enumerate(b_counts) if count
        },
        "survivor_ordered_profile_count": len(profile_rows),
        "survivor_ordered_profiles": profile_rows,
        "minimum_odd_secant_ranks": minimum_ranks,
        "survivor_minimum_odd_secant_ranks": survivor_minimum_ranks,
        "proved_residual_ii": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), default=-1)
    parser.add_argument("--blocks", type=int, default=4096)
    parser.add_argument("--threads", type=int, default=256)
    parser.add_argument("--start-rank", type=int, default=0)
    parser.add_argument("--stop-rank", type=int)
    parser.add_argument("--verify-cpu-prefix", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.c_h, args.blocks, args.threads, args.start_rank, args.stop_rank)
    if args.verify_cpu_prefix:
        verify_cpu_prefix(result)
    atomic_write(args.output, result)
    compact = {
        key: value
        for key, value in result.items()
        if key
        not in (
            "survivor_ordered_profiles",
            "minimum_odd_secant_ranks",
            "survivor_minimum_odd_secant_ranks",
        )
    }
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
