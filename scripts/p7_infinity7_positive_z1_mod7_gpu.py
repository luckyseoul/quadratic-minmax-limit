#!/usr/bin/env python3
"""Exact projected-syndrome CUDA scan for positive p7 infinity+7 with z=1.

There are four exact mean allocations. Either the unique b=7 direction is
raised from mean 0 to 8, or it stays at zero and one of the other three
directions of its quadratic type is raised from 8 to 16. The sole variable
catalog has 1,764 rows (2,233 for b=3 at mean 16).

The GPU tests 23 direction-specific mod-seven dependencies by lossless
base-seven packed keys.
Every projected survivor is then checked on the host against all 135
dependencies and the complete catalog. Zero exact survivors excludes z=1.
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import field_direction_data, projective_directions  # noqa: E402
from p7_infinity7_positive_z0_mod7_gpu import choose_table  # noqa: E402
from p7_size6_positive_infinity_mod7_gpu import dependency_tables, finite_labels  # noqa: E402
from p7_unsaturated_mod7_batch import contribution_matrix  # noqa: E402
from p7_unsaturated_modular_catalog_filter import equation_matrix, left_dependencies  # noqa: E402


PROJECTED_DEPENDENCIES = 23
MAX_CATALOG = 2233
SURVIVOR_CAPACITY = 1_000_000
EXPECTED_Z1_BOUNDARIES = 6_324_528
EXPECTED_PROJECTED_SURVIVORS = 1_326
EXPECTED_PROJECTED_SHA256 = "23de1f85d34f641d06279e8cdbc17fe6615fcd98198885e18f695d1812982b4c"


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def independent_dependency_rows(dependencies: np.ndarray, direction: int) -> np.ndarray:
    """Choose dependency coordinates of full rank on one direction block."""
    block = dependencies[:, 2 + 35 * direction : 2 + 35 * (direction + 1)] % 7
    matrix = block.T.astype(np.int64).copy()
    pivot_columns = []
    rank = 0
    for column in range(matrix.shape[1]):
        candidates = np.flatnonzero(matrix[rank:, column])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        matrix[[rank, pivot]] = matrix[[pivot, rank]]
        matrix[rank] = matrix[rank] * pow(int(matrix[rank, column]), -1, 7) % 7
        for row in range(matrix.shape[0]):
            if row != rank and matrix[row, column]:
                matrix[row] = (matrix[row] - matrix[row, column] * matrix[rank]) % 7
        pivot_columns.append(column)
        rank += 1
        if rank == matrix.shape[0]:
            break
    if rank != 23:
        raise AssertionError(f"direction block rank changed: {rank}")
    return np.asarray(pivot_columns, dtype=np.int16)


def pack_rows(rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(rows, dtype=np.uint64).T
    first = np.zeros(len(values), dtype=np.uint64)
    for coordinate in range(22):
        first = 7 * first + values[:, coordinate]
    return first, values[:, 22].astype(np.uint8)


def catalog_tables(
    dependencies: np.ndarray, selected_rows: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    keys = np.full((8, 128, MAX_CATALOG), np.iinfo(np.uint64).max, dtype=np.uint64)
    tails = np.full((8, 128, MAX_CATALOG), 255, dtype=np.uint8)
    counts = np.zeros((8, 128), dtype=np.int16)
    count_histogram: dict[int, int] = {}
    for direction in range(8):
        for mask in range(1, 128):
            b = mask.bit_count()
            if b == 7:
                mean = 8
            elif b in (1, 3, 5):
                mean = 16
            else:
                continue
            contribution = contribution_matrix(
                dependencies[selected_rows[direction]],
                direction,
                b,
                0,
                mean,
                {index for index in range(7) if mask & (1 << index)},
            )
            packed, tail = pack_rows(contribution)
            order = np.lexsort((tail, packed))
            packed = packed[order]
            tail = tail[order]
            if len(packed) > MAX_CATALOG:
                raise AssertionError("projected catalog exceeds allocation")
            keys[direction, mask, : len(packed)] = packed
            tails[direction, mask, : len(packed)] = tail
            counts[direction, mask] = len(packed)
            count_histogram[len(packed)] = count_histogram.get(len(packed), 0) + 1
    expected = {1764: 8 * (7 + 21 + 1), 2233: 8 * 35}
    if count_histogram != expected:
        raise AssertionError(f"unexpected variable catalog histogram {count_histogram}")
    return keys, tails, counts, {str(key): value for key, value in sorted(count_histogram.items())}


def zero_contributions(dependencies: np.ndarray) -> np.ndarray:
    out = np.empty((8, 135), dtype=np.int16)
    for direction in range(8):
        block = dependencies[:, 2 + 35 * direction : 2 + 35 * (direction + 1)]
        out[direction] = (13 * np.sum(block, axis=1, dtype=np.int64) % 7).astype(np.int16)
    return out


KERNEL = r'''
__device__ bool contains_key(
    const unsigned long long* rows, const unsigned char* tails,
    int count, unsigned long long key, unsigned char tail) {
    int lo = 0, hi = count;
    while (lo < hi) {
        const int mid = (lo + hi) >> 1;
        const unsigned long long value = rows[mid];
        if (value < key || (value == key && tails[mid] < tail)) lo = mid + 1;
        else hi = mid;
    }
    return lo < count && rows[lo] == key && tails[lo] == tail;
}

extern "C" __global__ void scan_z1(
    const short* labels,
    const signed char* epsilons,
    const short* base,
    const short* floor_tables,
    const short* zero_tables,
    const short* selected_rows,
    const unsigned long long* catalog_keys,
    const unsigned char* catalog_tails,
    const short* catalog_counts,
    const unsigned long long* choose,
    const unsigned long long total,
    unsigned long long* checked,
    unsigned long long* z1_count,
    unsigned long long* projected_survivor_count,
    unsigned long long* projected_survivor_ranks,
    const unsigned long long survivor_capacity)
{
    const int coordinate = threadIdx.x;
    __shared__ int index[7];
    __shared__ int masks[8];
    __shared__ int undetermined;
    __shared__ int z;
    for (unsigned long long rank0 = blockIdx.x; rank0 < total; rank0 += gridDim.x) {
        if (coordinate == 0) {
            unsigned long long rank = rank0;
            int next = 0;
            #pragma unroll
            for (int position = 0; position < 7; ++position) {
                const int remaining = 6-position;
                const int last = 49-(remaining+1);
                for (int candidate = next; candidate <= last; ++candidate) {
                    const unsigned long long ways = choose[(49-candidate-1)*8+remaining];
                    if (rank < ways) { index[position]=candidate; next=candidate+1; break; }
                    rank -= ways;
                }
            }
            z = 0;
            undetermined = -1;
        }
        __syncthreads();
        if (coordinate < 8) {
            int mask = 0;
            #pragma unroll
            for (int column=0; column<7; ++column)
                mask ^= 1 << labels[49*coordinate+index[column]];
            masks[coordinate] = mask;
        }
        __syncthreads();
        if (coordinate == 0) {
            for (int direction=0; direction<8; ++direction) {
                if (__popc((unsigned int)masks[direction]) == 7) {
                    ++z;
                    undetermined = direction;
                }
            }
            atomicAdd(checked, 1ULL);
            if (z == 1) {
                atomicAdd(z1_count, 1ULL);
                bool passing = false;
                unsigned long long key = 0ULL;
                unsigned char tail = 0;
                // Allocation A: all means eight; b=7 is the variable catalog.
                for (int k=0; k<23; ++k) {
                    const int row = selected_rows[undetermined*23+k];
                    int target = -base[row];
                    for (int d=0; d<8; ++d)
                        if (d != undetermined)
                            target -= floor_tables[(d*128+masks[d])*135+row];
                    target %= 7;
                    if (target < 0) target += 7;
                    if (k < 22) key = 7ULL*key + (unsigned long long)target;
                    else tail = (unsigned char)target;
                }
                int count = catalog_counts[undetermined*128+masks[undetermined]];
                const unsigned long long* rows = catalog_keys +
                    ((undetermined*128+masks[undetermined])*2233);
                const unsigned char* row_tails = catalog_tails +
                    ((undetermined*128+masks[undetermined])*2233);
                if (contains_key(rows, row_tails, count, key, tail)) passing = true;

                // Allocations B: b=7 stays at zero and one same-type ordinary
                // direction rises from mean eight to sixteen.
                for (int elevated=0; elevated<8 && !passing; ++elevated) {
                    if (elevated == undetermined || epsilons[elevated] != epsilons[undetermined])
                        continue;
                    key = 0ULL;
                    tail = 0;
                    for (int k=0; k<23; ++k) {
                        const int row = selected_rows[elevated*23+k];
                        int target = -base[row] - zero_tables[undetermined*135+row];
                        for (int d=0; d<8; ++d)
                            if (d != undetermined && d != elevated)
                                target -= floor_tables[(d*128+masks[d])*135+row];
                        target %= 7;
                        if (target < 0) target += 7;
                        if (k < 22) key = 7ULL*key + (unsigned long long)target;
                        else tail = (unsigned char)target;
                    }
                    count = catalog_counts[elevated*128+masks[elevated]];
                    rows = catalog_keys + ((elevated*128+masks[elevated])*2233);
                    row_tails = catalog_tails + ((elevated*128+masks[elevated])*2233);
                    if (contains_key(rows, row_tails, count, key, tail)) passing = true;
                }
                if (passing) {
                    const unsigned long long slot = atomicAdd(projected_survivor_count, 1ULL);
                    if (slot < survivor_capacity) projected_survivor_ranks[slot] = rank0;
                }
            }
        }
        __syncthreads();
    }
}
'''


def unrank(rank: int) -> tuple[int, ...]:
    out = []
    next_value = 0
    for position in range(7):
        remaining = 6-position
        for candidate in range(next_value, 49-remaining):
            ways = math.comb(49-candidate-1, remaining)
            if rank < ways:
                out.append(candidate)
                next_value = candidate+1
                break
            rank -= ways
    return tuple(out)


def exact_survives(
    rank: int,
    labels: np.ndarray,
    epsilons: np.ndarray,
    base: np.ndarray,
    floor_tables: np.ndarray,
    zero_tables: np.ndarray,
    dependencies: np.ndarray,
) -> bool:
    boundary = unrank(rank)
    masks = []
    for direction in range(8):
        mask = 0
        for point in boundary:
            mask ^= 1 << int(labels[direction, point])
        masks.append(mask)
    undetermined = [d for d, mask in enumerate(masks) if mask.bit_count() == 7]
    if len(undetermined) != 1:
        raise AssertionError("projected survivor is not z=1")
    u = undetermined[0]

    fixed = base.astype(np.int16).copy()
    for d, mask in enumerate(masks):
        if d != u:
            fixed += floor_tables[d, mask]
    target = (-fixed) % 7
    contribution = contribution_matrix(dependencies, u, 7, 0, 8, set(range(7)))
    if any(np.array_equal(contribution[:, i], target) for i in range(contribution.shape[1])):
        return True

    for elevated in range(8):
        if elevated == u or epsilons[elevated] != epsilons[u]:
            continue
        fixed = base.astype(np.int16).copy() + zero_tables[u]
        for d, mask in enumerate(masks):
            if d != u and d != elevated:
                fixed += floor_tables[d, mask]
        target = (-fixed) % 7
        mask = masks[elevated]
        contribution = contribution_matrix(
            dependencies,
            elevated,
            mask.bit_count(),
            0,
            16,
            {index for index in range(7) if mask & (1 << index)},
        )
        if any(np.array_equal(contribution[:, i], target) for i in range(contribution.shape[1])):
            return True
    return False


def run(blocks: int) -> dict:
    import cupy as cp

    started = time.time()
    matrix = equation_matrix()
    rank, dependencies = left_dependencies(matrix, 7)
    base, floor_tables, linear = dependency_tables()
    labels = finite_labels()
    epsilons = np.asarray(
        [field_direction_data(7, direction)[0] for direction in projective_directions(7)],
        dtype=np.int8,
    )
    zero_tables = zero_contributions(dependencies)
    selected_rows = np.stack(
        [independent_dependency_rows(dependencies, direction) for direction in range(8)]
    )
    keys, tails, counts, catalog_histogram = catalog_tables(dependencies, selected_rows)
    total = math.comb(49, 7)
    kernel = cp.RawKernel(KERNEL, "scan_z1", options=("--std=c++11",))
    checked = cp.zeros(1, dtype=cp.uint64)
    z1_count = cp.zeros(1, dtype=cp.uint64)
    projected_count = cp.zeros(1, dtype=cp.uint64)
    projected_ranks = cp.zeros(SURVIVOR_CAPACITY, dtype=cp.uint64)
    kernel(
        (blocks,),
        (32,),
        (
            cp.asarray(labels, dtype=cp.int16),
            cp.asarray(epsilons),
            cp.asarray(base, dtype=cp.int16),
            cp.asarray(floor_tables, dtype=cp.int16),
            cp.asarray(zero_tables, dtype=cp.int16),
            cp.asarray(selected_rows, dtype=cp.int16),
            cp.asarray(keys),
            cp.asarray(tails),
            cp.asarray(counts),
            cp.asarray(choose_table()),
            total,
            checked,
            z1_count,
            projected_count,
            projected_ranks,
            SURVIVOR_CAPACITY,
        ),
    )
    cp.cuda.runtime.deviceSynchronize()
    checked_value = int(cp.asnumpy(checked)[0])
    z1_value = int(cp.asnumpy(z1_count)[0])
    projected_value = int(cp.asnumpy(projected_count)[0])
    if checked_value != total or projected_value > SURVIVOR_CAPACITY:
        raise AssertionError(
            f"CUDA z1 scan incomplete or overflowed: checked={checked_value} "
            f"total={total} projected={projected_value} capacity={SURVIVOR_CAPACITY}"
        )
    ranks = cp.asnumpy(projected_ranks)[:projected_value].astype(int).tolist()
    sorted_ranks = sorted(ranks)
    rank_bytes = json.dumps(sorted_ranks, separators=(",", ":")).encode()
    projected_sha256 = hashlib.sha256(rank_bytes).hexdigest()
    if (
        z1_value != EXPECTED_Z1_BOUNDARIES
        or projected_value != EXPECTED_PROJECTED_SURVIVORS
        or projected_sha256 != EXPECTED_PROJECTED_SHA256
    ):
        raise AssertionError(
            "z1 projected census changed: "
            f"z1={z1_value}, projected={projected_value}, sha256={projected_sha256}"
        )
    exact = [
        rank0
        for rank0 in sorted_ranks
        if exact_survives(rank0, labels, epsilons, base, floor_tables, zero_tables, dependencies)
    ]
    properties = cp.cuda.runtime.getDeviceProperties(0)
    name = properties["name"]
    if isinstance(name, bytes):
        name = name.decode()
    return {
        "experiment": "p7_infinity7_positive_z1_mod7_gpu",
        "status": "complete_projected_then_exact_mod_seven_z1_exhaustion",
        "p": 7,
        "c_H": 1,
        "all_boundaries": total,
        "checked_boundaries": checked_value,
        "z1_boundaries": z1_value,
        "mean_allocation_count_per_boundary": 4,
        "projected_dependency_count": PROJECTED_DEPENDENCIES,
        "projected_survivors": projected_value,
        "projected_survivor_rank_sha256": projected_sha256,
        "first_projected_survivor_ranks": sorted_ranks[:64],
        "all_dependency_survivors": len(exact),
        "all_dependency_survivor_ranks": exact,
        "z1_branch_excluded": len(exact) == 0,
        "linear_system": linear,
        "rank_mod_7": rank,
        "catalog_row_histogram_by_direction_mask": catalog_histogram,
        "backend": "CUDA/CuPy plus exact NumPy host validation",
        "device": str(name),
        "blocks": blocks,
        "elapsed_seconds": time.time()-started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blocks", type=int, default=65535)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = run(args.blocks)
    atomic_write(args.output, out)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
