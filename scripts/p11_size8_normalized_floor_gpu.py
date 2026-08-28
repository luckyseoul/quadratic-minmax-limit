#!/usr/bin/env python3
"""Exact normalized CUDA census for the p=11 eight-finite survivor.

Every finite eight-set has an affine image containing field points 0 and 1.
The normalizing scalar may swap quadratic direction types, so the kernel
tests both values of c_H.  It therefore suffices, for exclusion, to scan the
``C(119,6)`` sets containing 0 and 1.  This is an exact boundary-level test
of Proposition 15.632's split parity floors, not an edge-lift model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from itertools import islice
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)


P = 11
Q = P * P
BOUNDARY_SIZE = 8
FREE_POINTS = Q - 2
FREE_CHOICE = BOUNDARY_SIZE - 2
TYPE_BUDGET = (P + 1) ** 2 // 2
COST_STRIDE = (P + 1) * 2 * P // 2 + 1


def atomic_write(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def direction_tables() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    labels = []
    epsilons = []
    for direction in projective_directions(P):
        epsilon, row = field_direction_data(P, direction)
        epsilons.append(int(epsilon))
        labels.append([int(value) for value in row])
    label_array = np.asarray(labels, dtype=np.int8)
    epsilon_array = np.asarray(epsilons, dtype=np.int8)
    floors = np.asarray(
        [
            [scaled_direction_floor(P, b, phase) for b in range(P + 1)]
            for phase in (0, 1)
        ],
        dtype=np.int16,
    )
    if label_array.shape != (P + 1, Q):
        raise AssertionError("unexpected direction-label shape")
    if sorted(epsilon_array.tolist()) != [-1] * 6 + [1] * 6:
        raise AssertionError("unexpected quadratic direction split")
    return label_array, epsilon_array, floors


def choose_table(n: int, support: int) -> np.ndarray:
    table = np.zeros((n + 1, support + 1), dtype=np.uint64)
    for a in range(n + 1):
        for b in range(min(a, support) + 1):
            table[a, b] = math.comb(a, b)
    return table


def unrank_free(rank0: int) -> tuple[int, ...]:
    rank = rank0
    out = [0, 1]
    next_value = 0
    for position in range(FREE_CHOICE):
        remaining = FREE_CHOICE - 1 - position
        last = FREE_POINTS - (remaining + 1)
        for candidate in range(next_value, last + 1):
            ways = math.comb(FREE_POINTS - candidate - 1, remaining)
            if rank < ways:
                out.append(candidate + 2)
                next_value = candidate + 1
                break
            rank -= ways
        else:
            raise ArithmeticError("combination rank did not decode")
    return tuple(out)


def boundary_costs(
    boundary: tuple[int, ...],
    labels: np.ndarray,
    epsilons: np.ndarray,
    floors: np.ndarray,
    c_h: int,
) -> tuple[int, int, list[int]]:
    totals = {-1: 0, 1: 0}
    b_values = []
    for direction in range(P + 1):
        mask = 0
        for point in boundary:
            mask ^= 1 << int(labels[direction, point])
        b = mask.bit_count()
        phase = int(int(epsilons[direction]) == c_h)
        totals[int(epsilons[direction])] += int(floors[phase, b])
        b_values.append(b)
    return totals[-1], totals[1], b_values


KERNEL = r'''
extern "C" __global__
void p11_size8_normalized_floor(
    const signed char* labels,
    const signed char* epsilons,
    const short* floors,
    const unsigned long long* choose,
    const int choose_stride,
    const unsigned long long start_rank,
    const unsigned long long stop_rank,
    const int cost_stride,
    unsigned long long* checked,
    unsigned long long* survivors,
    unsigned long long* histograms,
    unsigned long long* minimum_keys)
{
    const unsigned long long logical_thread =
        (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned long long logical_threads =
        (unsigned long long)gridDim.x * blockDim.x;
    unsigned long long local_checked = 0;
    unsigned long long local_survivors[2] = {0ULL, 0ULL};

    for (unsigned long long rank0 = start_rank + logical_thread;
         rank0 < stop_rank;
         rank0 += logical_threads)
    {
        int index[8];
        index[0] = 0;
        index[1] = 1;
        unsigned long long rank = rank0;
        int next = 0;
        bool valid = true;
        #pragma unroll
        for (int position = 0; position < 6; ++position)
        {
            const int remaining = 5 - position;
            const int last = 119 - (remaining + 1);
            bool selected = false;
            for (int candidate = next; candidate <= last; ++candidate)
            {
                const unsigned long long ways =
                    choose[(119 - candidate - 1) * choose_stride + remaining];
                if (rank < ways)
                {
                    index[position + 2] = candidate + 2;
                    next = candidate + 1;
                    selected = true;
                    break;
                }
                rank -= ways;
            }
            if (!selected) valid = false;
        }
        if (!valid) continue;

        int b_values[12];
        #pragma unroll
        for (int direction = 0; direction < 12; ++direction)
        {
            unsigned int mask = 0U;
            #pragma unroll
            for (int column = 0; column < 8; ++column)
                mask ^= 1U << (unsigned int)labels[121 * direction + index[column]];
            b_values[direction] = __popc(mask);
        }

        #pragma unroll
        for (int sign_index = 0; sign_index < 2; ++sign_index)
        {
            const int c_h = sign_index ? 1 : -1;
            int cost_minus = 0;
            int cost_plus = 0;
            #pragma unroll
            for (int direction = 0; direction < 12; ++direction)
            {
                const int phase = epsilons[direction] == c_h;
                const int cost = floors[12 * phase + b_values[direction]];
                if (epsilons[direction] < 0) cost_minus += cost;
                else cost_plus += cost;
            }
            atomicAdd(
                histograms
                    + sign_index * cost_stride * cost_stride
                    + cost_minus * cost_stride + cost_plus,
                1ULL);
            const int maximum_cost = cost_minus > cost_plus ? cost_minus : cost_plus;
            const unsigned long long key =
                ((unsigned long long)maximum_cost << 32) | rank0;
            atomicMin(minimum_keys + sign_index, key);
            if (cost_minus <= 72 && cost_plus <= 72)
                ++local_survivors[sign_index];
        }
        ++local_checked;
    }
    if (local_checked) atomicAdd(checked, local_checked);
    if (local_survivors[0]) atomicAdd(survivors, local_survivors[0]);
    if (local_survivors[1]) atomicAdd(survivors + 1, local_survivors[1]);
}
'''


def launch(
    start_rank: int,
    stop_rank: int,
    blocks: int,
    threads: int,
) -> dict[str, object]:
    import cupy as cp

    labels, epsilons, floors = direction_tables()
    kernel = cp.RawKernel(
        KERNEL, "p11_size8_normalized_floor", options=("--std=c++11",)
    )
    checked_gpu = cp.zeros(1, dtype=cp.uint64)
    survivors_gpu = cp.zeros(2, dtype=cp.uint64)
    histograms_gpu = cp.zeros(
        (2, COST_STRIDE, COST_STRIDE), dtype=cp.uint64
    )
    minimum_keys_gpu = cp.full(2, np.uint64((1 << 64) - 1), dtype=cp.uint64)
    started = time.perf_counter()
    kernel(
        (blocks,),
        (threads,),
        (
            cp.asarray(labels),
            cp.asarray(epsilons),
            cp.asarray(floors),
            cp.asarray(choose_table(FREE_POINTS, FREE_CHOICE)),
            np.int32(FREE_CHOICE + 1),
            np.uint64(start_rank),
            np.uint64(stop_rank),
            np.int32(COST_STRIDE),
            checked_gpu,
            survivors_gpu,
            histograms_gpu,
            minimum_keys_gpu,
        ),
    )
    cp.cuda.runtime.deviceSynchronize()
    elapsed = time.perf_counter() - started
    histograms = cp.asnumpy(histograms_gpu)
    minimum_keys = cp.asnumpy(minimum_keys_gpu)
    properties = cp.cuda.runtime.getDeviceProperties(0)
    device = properties["name"]
    if isinstance(device, bytes):
        device = device.decode()
    return {
        "checked": int(cp.asnumpy(checked_gpu)[0]),
        "survivors": [int(value) for value in cp.asnumpy(survivors_gpu)],
        "histograms": histograms,
        "minimum_keys": [int(value) for value in minimum_keys],
        "elapsed_seconds": elapsed,
        "boundaries_per_second": (stop_rank - start_rank) / elapsed,
        "device": device,
    }


def independent_cpu_prefix(stop_rank: int) -> dict[str, object]:
    from itertools import combinations

    labels, epsilons, floors = direction_tables()
    histograms = np.zeros((2, COST_STRIDE, COST_STRIDE), dtype=np.uint64)
    survivors = [0, 0]
    minimum_keys = [(1 << 64) - 1, (1 << 64) - 1]
    combinations_iter = combinations(range(2, Q), FREE_CHOICE)
    for rank, free in enumerate(islice(combinations_iter, stop_rank)):
        boundary = (0, 1, *free)
        for sign_index, c_h in enumerate((-1, 1)):
            cost_minus, cost_plus, _b = boundary_costs(
                boundary, labels, epsilons, floors, c_h
            )
            histograms[sign_index, cost_minus, cost_plus] += 1
            key = (max(cost_minus, cost_plus) << 32) | rank
            minimum_keys[sign_index] = min(minimum_keys[sign_index], key)
            if cost_minus <= TYPE_BUDGET and cost_plus <= TYPE_BUDGET:
                survivors[sign_index] += 1
    return {
        "checked": stop_rank,
        "survivors": survivors,
        "histograms": histograms,
        "minimum_keys": minimum_keys,
    }


def histogram_rows(histogram: np.ndarray) -> list[dict[str, int]]:
    rows = []
    for cost_minus, cost_plus in zip(*np.nonzero(histogram)):
        rows.append(
            {
                "cost_minus": int(cost_minus),
                "cost_plus": int(cost_plus),
                "count": int(histogram[cost_minus, cost_plus]),
            }
        )
    rows.sort(key=lambda row: (max(row["cost_minus"], row["cost_plus"]), row["cost_minus"], row["cost_plus"]))
    return rows


def run(blocks: int, threads: int, prefix: int) -> dict[str, object]:
    total = math.comb(FREE_POINTS, FREE_CHOICE)
    prefix_stop = min(prefix, total)
    prefix_gpu = launch(0, prefix_stop, blocks, threads)
    prefix_cpu = independent_cpu_prefix(prefix_stop)
    prefix_match = bool(
        prefix_gpu["checked"] == prefix_cpu["checked"]
        and prefix_gpu["survivors"] == prefix_cpu["survivors"]
        and prefix_gpu["minimum_keys"] == prefix_cpu["minimum_keys"]
        and np.array_equal(prefix_gpu["histograms"], prefix_cpu["histograms"])
    )
    if not prefix_match:
        raise AssertionError("independent CPU/GPU prefix audit failed")

    full = launch(0, total, blocks, threads)
    labels, epsilons, floors = direction_tables()
    signs = []
    for sign_index, c_h in enumerate((-1, 1)):
        key = int(full["minimum_keys"][sign_index])
        minimum_cost = key >> 32
        minimum_rank = key & ((1 << 32) - 1)
        boundary = unrank_free(minimum_rank)
        cost_minus, cost_plus, b_values = boundary_costs(
            boundary, labels, epsilons, floors, c_h
        )
        if max(cost_minus, cost_plus) != minimum_cost:
            raise AssertionError("minimum-rank CPU audit failed")
        histogram = full["histograms"][sign_index]
        signs.append(
            {
                "c_H": c_h,
                "floor_survivors": int(full["survivors"][sign_index]),
                "minimum_maximum_type_cost": minimum_cost,
                "minimum_budget_excess": minimum_cost - TYPE_BUDGET,
                "first_minimum_rank": minimum_rank,
                "first_minimum_boundary": list(boundary),
                "first_minimum_coordinates": [
                    [point % P, point // P] for point in boundary
                ],
                "first_minimum_type_costs": {
                    "-1": cost_minus,
                    "1": cost_plus,
                },
                "first_minimum_b_by_direction": b_values,
                "cost_pair_histogram": histogram_rows(histogram),
                "histogram_sha256": hashlib.sha256(
                    histogram.tobytes()
                ).hexdigest(),
            }
        )
    proved = bool(
        full["checked"] == total
        and prefix_match
        and all(row["floor_survivors"] == 0 for row in signs)
    )
    return {
        "experiment": "p11_size8_normalized_floor_gpu",
        "status": "exact_normalized_boundary_census",
        "p": P,
        "finite_boundary_size": BOUNDARY_SIZE,
        "normalization": "contains field points 0 and 1",
        "coverage": (
            "every finite eight-set has an affine image containing 0 and 1; "
            "a nonsquare normalizing scalar swaps c_H, so both signs are tested"
        ),
        "normalized_boundaries": total,
        "checked_boundaries": int(full["checked"]),
        "type_budget": TYPE_BUDGET,
        "signs": signs,
        "prefix_verification": {
            "checked": prefix_stop,
            "independent_itertools_all_histograms_match": prefix_match,
        },
        "device": full["device"],
        "launch": {"blocks": blocks, "threads": threads},
        "elapsed_seconds": full["elapsed_seconds"],
        "boundaries_per_second": full["boundaries_per_second"],
        "all_p11_eight_finite_boundaries_excluded": proved,
        "closes_residual_ii": False,
        "L_status": "OPEN",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blocks", type=int, default=4096)
    parser.add_argument("--threads", type=int, default=256)
    parser.add_argument("--prefix", type=int, default=100_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = run(args.blocks, args.threads, args.prefix)
    if args.output is not None:
        atomic_write(args.output, out)
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
