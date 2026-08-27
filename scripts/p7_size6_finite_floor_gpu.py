#!/usr/bin/env python3
"""Exact floor-budget sieve for six finite boundary points at p=7.

For a six-finite boundary the odd-fibre counts are in ``{0,2,4,6}``.
The phase-zero scaled floors are ``(0,8,8,8)`` and the phase-one floors
are ``(14,6,14,6)``.  Each quadratic type has four directions and budget
32.  The phase-zero half always fits, while the phase-one half fits exactly
when at most one of its directions has odd-fibre count zero or four.

For no-infinity boundaries the phase-one quadratic type is ``eps=c_H``.
The CUDA and NumPy backends apply this integer criterion to every
``C(49,6)`` finite boundary and retain the complete canonical survivor list.
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
        eps, row = field_direction_data(7, direction)
        epsilons.append(int(eps))
        labels.append(tuple(int(value) for value in row))
    label_array = np.asarray(labels, dtype=np.int16)
    epsilon_array = np.asarray(epsilons, dtype=np.int8)
    if label_array.shape != (8, 49):
        raise AssertionError("unexpected p=7 direction-label table")
    if sorted(epsilon_array.tolist()) != [-1] * 4 + [1] * 4:
        raise AssertionError("unexpected p=7 quadratic-type split")
    return label_array, epsilon_array


def batches(size: int):
    source = itertools.combinations(range(49), 6)
    while True:
        rows = list(itertools.islice(source, size))
        if not rows:
            return
        yield np.asarray(rows, dtype=np.int16)


def scan_numpy(
    labels: np.ndarray, epsilons: np.ndarray, c_h: int, batch_size: int
) -> tuple[int, list[list[int]], dict[int, int]]:
    checked = 0
    survivors: list[list[int]] = []
    histogram = {0: 0, 2: 0, 4: 0, 6: 0}
    popcount = np.asarray([value.bit_count() for value in range(128)], dtype=np.int8)
    phase_one = epsilons == c_h
    for batch in batches(batch_size):
        sizes = np.empty((len(batch), 8), dtype=np.int8)
        for direction in range(8):
            fibre_labels = labels[direction, batch]
            masks = np.bitwise_xor.reduce(np.left_shift(1, fibre_labels), axis=1)
            sizes[:, direction] = popcount[masks]
        if np.any((sizes < 0) | (sizes > 6) | ((sizes & 1) != 0)):
            raise AssertionError("six points produced an invalid odd-fibre count")
        for value in histogram:
            histogram[value] += int(np.count_nonzero(sizes == value))
        costly = (sizes[:, phase_one] == 0) | (sizes[:, phase_one] == 4)
        passing = np.count_nonzero(costly, axis=1) <= 1
        survivors.extend(batch[passing].astype(int).tolist())
        checked += len(batch)
    return checked, survivors, histogram


def scan_cuda(
    labels: np.ndarray, epsilons: np.ndarray, c_h: int, batch_size: int
) -> tuple[int, list[list[int]], dict[int, int], str]:
    import cupy as cp

    if cp.cuda.runtime.getDeviceCount() < 1:
        raise RuntimeError("CUDA is unavailable")
    labels_gpu = cp.asarray(labels, dtype=cp.int16)
    eps_gpu = cp.asarray(epsilons, dtype=cp.int8)
    kernel = cp.RawKernel(
        r'''
        extern "C" __global__ void floor_sieve(
            const short* combinations,
            const short* labels,
            const signed char* epsilons,
            const int c_h,
            unsigned char* passing,
            unsigned long long* histogram,
            const int count) {
          const int row = blockDim.x * blockIdx.x + threadIdx.x;
          if (row >= count) return;
          int costly_phase_one = 0;
          bool valid = true;
          #pragma unroll
          for (int direction = 0; direction < 8; ++direction) {
            int mask = 0;
            #pragma unroll
            for (int column = 0; column < 6; ++column) {
              const int point = combinations[6 * row + column];
              mask ^= 1 << labels[49 * direction + point];
            }
            const int b = __popc((unsigned int)mask);
            if (b == 0) atomicAdd(histogram + 0, 1ULL);
            else if (b == 2) atomicAdd(histogram + 1, 1ULL);
            else if (b == 4) atomicAdd(histogram + 2, 1ULL);
            else if (b == 6) atomicAdd(histogram + 3, 1ULL);
            else valid = false;
            if (epsilons[direction] == c_h && (b == 0 || b == 4)) {
              ++costly_phase_one;
            }
          }
          passing[row] = valid && costly_phase_one <= 1;
        }
        ''',
        "floor_sieve",
        options=("--std=c++11",),
    )
    checked = 0
    survivors: list[list[int]] = []
    histogram = {0: 0, 2: 0, 4: 0, 6: 0}
    for batch in batches(batch_size):
        batch_gpu = cp.asarray(batch, dtype=cp.int16)
        passing_gpu = cp.empty(len(batch), dtype=cp.uint8)
        batch_histogram_gpu = cp.zeros(4, dtype=cp.uint64)
        kernel(
            ((len(batch) + 255) // 256,),
            (256,),
            (
                batch_gpu,
                labels_gpu,
                eps_gpu,
                c_h,
                passing_gpu,
                batch_histogram_gpu,
                len(batch),
            ),
        )
        passing = cp.asnumpy(passing_gpu).astype(bool)
        batch_histogram = cp.asnumpy(batch_histogram_gpu)
        survivors.extend(batch[passing].astype(int).tolist())
        for index, value in enumerate((0, 2, 4, 6)):
            histogram[value] += int(batch_histogram[index])
        checked += len(batch)
    properties = cp.cuda.runtime.getDeviceProperties(0)
    name = properties["name"]
    if isinstance(name, bytes):
        name = name.decode()
    return checked, survivors, histogram, str(name)


def run(backend: str, c_h: int, batch_size: int) -> dict:
    started = time.time()
    labels, epsilons = direction_tables()
    floors = {
        phase: {b: scaled_direction_floor(7, b, phase) for b in (0, 2, 4, 6)}
        for phase in (0, 1)
    }
    if floors != {0: {0: 0, 2: 8, 4: 8, 6: 8}, 1: {0: 14, 2: 6, 4: 14, 6: 6}}:
        raise AssertionError("unexpected p=7 six-finite floor table")
    if backend == "cuda":
        checked, survivors, histogram, device = scan_cuda(
            labels, epsilons, c_h, batch_size
        )
    else:
        checked, survivors, histogram = scan_numpy(
            labels, epsilons, c_h, batch_size
        )
        device = "CPU/NumPy"
    expected = math.comb(49, 6)
    if checked != expected:
        raise AssertionError("combination enumeration was incomplete")
    canonical = json.dumps(survivors, separators=(",", ":")).encode()
    return {
        "experiment": "p7_size6_finite_floor_gpu",
        "status": "complete_exact_integer_floor_budget_sieve",
        "p": 7,
        "c_H": c_h,
        "boundary_size": 6,
        "infinity_value": 0,
        "finite_boundary_size": 6,
        "backend": backend,
        "device": device,
        "batch_size": batch_size,
        "criterion": {
            "phase_one_quadratic_type": c_h,
            "scaled_floor_by_phase_and_odd_fibre_count": {
                str(phase): {str(b): value for b, value in row.items()}
                for phase, row in floors.items()
            },
            "directions_per_type": 4,
            "budget_per_type": 32,
            "maximum_phase_one_b0_or_b4_directions": 1,
        },
        "all_boundaries": expected,
        "checked_boundaries": checked,
        "direction_odd_fibre_histogram": {
            str(key): value for key, value in histogram.items()
        },
        "floor_surviving_boundaries": len(survivors),
        "floor_rejected_boundaries": checked - len(survivors),
        "survivor_sha256": hashlib.sha256(canonical).hexdigest(),
        "survivors_finite_field": survivors,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", choices=("cuda", "numpy"), default="cuda")
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--batch-size", type=int, default=262144)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = run(args.backend, args.c_h, args.batch_size)
    atomic_write(args.output, out)
    print(
        json.dumps(
            {key: value for key, value in out.items() if key != "survivors_finite_field"},
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
