#!/usr/bin/env python3
"""Complete CUDA mod-seven exclusion of positive p7 infinity+7 with z=0.

Here z counts undetermined directions (odd-fibre count seven). If z=0,
every direction has b in {1,3,5}, phase zero, and exact mean eight. The
Johnson slack catalog is unique for each mask, so every seven-point boundary
fixes the full right side of the common affine edge system. This scanner
unranks all C(49,7) boundaries on device and tests all 135 mod-seven left-null
dependencies. A zero survivor count is an exact exclusion of the z=0 branch.
"""
from __future__ import annotations

import argparse
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

from p7_size6_positive_infinity_mod7_gpu import dependency_tables, finite_labels  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def choose_table() -> np.ndarray:
    out = np.zeros((50, 8), dtype=np.uint64)
    for n in range(50):
        for k in range(8):
            if k <= n:
                out[n, k] = math.comb(n, k)
    return out


KERNEL = r'''
extern "C" __global__ void scan_z0(
    const short* labels,
    const short* base,
    const short* tables,
    const unsigned long long* choose,
    const unsigned long long total,
    unsigned long long* checked,
    unsigned long long* z_histogram,
    unsigned long long* z0_count,
    unsigned long long* survivor_count,
    unsigned long long* survivor_ranks,
    const unsigned long long survivor_capacity)
{
    const int coordinate = threadIdx.x;
    __shared__ int index[7];
    __shared__ int masks[8];
    __shared__ int valid;
    __shared__ int passing;
    __shared__ int z;
    for (unsigned long long rank0 = blockIdx.x; rank0 < total; rank0 += gridDim.x) {
        if (coordinate == 0) {
            unsigned long long rank = rank0;
            int next = 0;
            valid = 1;
            passing = 1;
            z = 0;
            #pragma unroll
            for (int position = 0; position < 7; ++position) {
                const int remaining = 6 - position;
                bool selected = false;
                const int last = 49 - (remaining + 1);
                for (int candidate = next; candidate <= last; ++candidate) {
                    const unsigned long long ways = choose[(49-candidate-1)*8+remaining];
                    if (rank < ways) {
                        index[position] = candidate;
                        next = candidate + 1;
                        selected = true;
                        break;
                    }
                    rank -= ways;
                }
                if (!selected) valid = 0;
            }
        }
        __syncthreads();
        if (coordinate < 8) {
            int mask = 0;
            #pragma unroll
            for (int column = 0; column < 7; ++column)
                mask ^= 1 << labels[49*coordinate + index[column]];
            masks[coordinate] = mask;
            const int b = __popc((unsigned int)mask);
            if (b == 7) {
                atomicAdd(&z, 1);
                atomicExch(&valid, 0);
            } else if (b != 1 && b != 3 && b != 5) atomicExch(&valid, 0);
        }
        __syncthreads();
        if (coordinate < 135 && valid) {
            int syndrome = base[coordinate];
            #pragma unroll
            for (int direction = 0; direction < 8; ++direction)
                syndrome += tables[(direction*128+masks[direction])*135+coordinate];
            if (syndrome % 7 != 0) atomicExch(&passing, 0);
        }
        __syncthreads();
        if (coordinate == 0) {
            atomicAdd(checked, 1ULL);
            atomicAdd(z_histogram + z, 1ULL);
            if (valid) {
                atomicAdd(z0_count, 1ULL);
                if (passing) {
                    const unsigned long long slot = atomicAdd(survivor_count, 1ULL);
                    if (slot < survivor_capacity) survivor_ranks[slot] = rank0;
                }
            }
        }
        __syncthreads();
    }
}
'''


def cpu_prefix(limit: int, base: np.ndarray, tables: np.ndarray, labels: np.ndarray) -> dict:
    checked = 0
    z0 = 0
    survivors = []
    for rank, boundary in enumerate(itertools.islice(itertools.combinations(range(49), 7), limit)):
        masks = []
        valid = True
        for direction in range(8):
            mask = 0
            for point in boundary:
                mask ^= 1 << int(labels[direction, point])
            if mask.bit_count() not in (1, 3, 5):
                valid = False
                break
            masks.append(mask)
        checked += 1
        if not valid:
            continue
        z0 += 1
        syndrome = base.astype(np.int64).copy()
        for direction, mask in enumerate(masks):
            syndrome += tables[direction, mask]
        if np.all(syndrome % 7 == 0):
            survivors.append(rank)
    return {"checked": checked, "z0": z0, "survivor_ranks": survivors}


def run(blocks: int, prefix: int) -> dict:
    import cupy as cp

    started = time.time()
    base, tables, linear = dependency_tables()
    labels = finite_labels()
    choose = choose_table()
    total = math.comb(49, 7)
    kernel = cp.RawKernel(KERNEL, "scan_z0", options=("--std=c++11",))
    checked = cp.zeros(1, dtype=cp.uint64)
    z_histogram = cp.zeros(9, dtype=cp.uint64)
    z0_count = cp.zeros(1, dtype=cp.uint64)
    survivor_count = cp.zeros(1, dtype=cp.uint64)
    survivor_ranks = cp.zeros(1024, dtype=cp.uint64)
    kernel(
        (blocks,),
        (160,),
        (
            cp.asarray(labels, dtype=cp.int16),
            cp.asarray(base, dtype=cp.int16),
            cp.asarray(tables, dtype=cp.int16),
            cp.asarray(choose),
            total,
            checked,
            z_histogram,
            z0_count,
            survivor_count,
            survivor_ranks,
            1024,
        ),
    )
    cp.cuda.runtime.deviceSynchronize()
    checked_value = int(cp.asnumpy(checked)[0])
    z_histogram_value = cp.asnumpy(z_histogram).astype(int).tolist()
    z0_value = int(cp.asnumpy(z0_count)[0])
    survivor_value = int(cp.asnumpy(survivor_count)[0])
    stored = cp.asnumpy(survivor_ranks)[: min(survivor_value, 1024)].astype(int).tolist()
    if checked_value != total:
        raise AssertionError("CUDA rank interval incomplete")
    verification = cpu_prefix(prefix, base, tables, labels)
    if survivor_value == 0 and verification["survivor_ranks"]:
        raise AssertionError("CPU prefix found a survivor missed by CUDA")
    properties = cp.cuda.runtime.getDeviceProperties(0)
    name = properties["name"]
    if isinstance(name, bytes):
        name = name.decode()
    return {
        "experiment": "p7_infinity7_positive_z0_mod7_gpu",
        "status": "complete_exact_mod_seven_z0_boundary_exhaustion",
        "p": 7,
        "c_H": 1,
        "finite_boundary_size": 7,
        "all_boundaries": total,
        "checked_boundaries": checked_value,
        "z0_boundaries": z0_value,
        "undetermined_direction_histogram": {
            str(index): count for index, count in enumerate(z_histogram_value) if count
        },
        "mod7_survivors": survivor_value,
        "stored_survivor_ranks": stored,
        "z0_branch_excluded": survivor_value == 0,
        "linear_system": linear,
        "base_sha256": hashlib.sha256(np.ascontiguousarray(base).tobytes()).hexdigest(),
        "tables_sha256": hashlib.sha256(np.ascontiguousarray(tables).tobytes()).hexdigest(),
        "cpu_prefix_verification": verification,
        "backend": "CUDA/CuPy",
        "device": str(name),
        "blocks": blocks,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blocks", type=int, default=65535)
    parser.add_argument("--cpu-prefix", type=int, default=100000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = run(args.blocks, args.cpu_prefix)
    atomic_write(args.output, out)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
