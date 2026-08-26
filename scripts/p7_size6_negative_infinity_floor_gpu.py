#!/usr/bin/env python3
"""Exact floor-budget sieve for p=7 negative infinity-plus-five boundaries.

For ``c_H=-1`` with infinity in the boundary, every one of the eight
projective directions has phase one.  A five-point finite boundary has odd
fibre count ``b`` in ``{1,3,5}``; its scaled floor is 14 for ``b=3`` and 6
otherwise.  Since each quadratic type has four directions and budget 32,
the boundary survives Proposition 15.632 exactly when each type has at most
one direction with ``b=3``.

The CUDA and NumPy backends evaluate that integer criterion for every
``C(49,5)`` boundary.  The retained list is intended as an implementation-
independent cross-check of the serial symmetry-orbit catalog.
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
    source = itertools.combinations(range(49), 5)
    while True:
        rows = list(itertools.islice(source, size))
        if not rows:
            return
        yield np.asarray(rows, dtype=np.int16)


def scan_numpy(
    labels: np.ndarray, epsilons: np.ndarray, batch_size: int
) -> tuple[int, list[list[int]], dict[int, int]]:
    checked = 0
    survivors: list[list[int]] = []
    histogram = {1: 0, 3: 0, 5: 0}
    popcount = np.asarray([value.bit_count() for value in range(128)], dtype=np.int8)
    for batch in batches(batch_size):
        sizes = np.empty((len(batch), 8), dtype=np.int8)
        for direction in range(8):
            fibre_labels = labels[direction, batch]
            masks = np.bitwise_xor.reduce(np.left_shift(1, fibre_labels), axis=1)
            sizes[:, direction] = popcount[masks]
        if np.any((sizes != 1) & (sizes != 3) & (sizes != 5)):
            raise AssertionError("five points produced an even odd-fibre count")
        for value in histogram:
            histogram[value] += int(np.count_nonzero(sizes == value))
        passing = np.ones(len(batch), dtype=bool)
        for eps in (-1, 1):
            passing &= np.count_nonzero(sizes[:, epsilons == eps] == 3, axis=1) <= 1
        survivors.extend(batch[passing].astype(int).tolist())
        checked += len(batch)
    return checked, survivors, histogram


def scan_cuda(
    labels: np.ndarray, epsilons: np.ndarray, batch_size: int
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
            unsigned char* passing,
            unsigned long long* histogram,
            const int count) {
          const int row = blockDim.x * blockIdx.x + threadIdx.x;
          if (row >= count) return;
          int b3_minus = 0;
          int b3_plus = 0;
          bool valid = true;
          #pragma unroll
          for (int direction = 0; direction < 8; ++direction) {
            int mask = 0;
            #pragma unroll
            for (int column = 0; column < 5; ++column) {
              const int point = combinations[5 * row + column];
              mask ^= 1 << labels[49 * direction + point];
            }
            const int b = __popc((unsigned int)mask);
            if (b == 1) atomicAdd(histogram + 0, 1ULL);
            else if (b == 3) {
              atomicAdd(histogram + 1, 1ULL);
              if (epsilons[direction] < 0) ++b3_minus;
              else ++b3_plus;
            } else if (b == 5) atomicAdd(histogram + 2, 1ULL);
            else valid = false;
          }
          passing[row] = valid && b3_minus <= 1 && b3_plus <= 1;
        }
        ''',
        "floor_sieve",
        options=("--std=c++11",),
    )
    checked = 0
    survivors: list[list[int]] = []
    histogram = {1: 0, 3: 0, 5: 0}
    for batch in batches(batch_size):
        batch_gpu = cp.asarray(batch, dtype=cp.int16)
        passing_gpu = cp.empty(len(batch), dtype=cp.uint8)
        batch_histogram_gpu = cp.zeros(3, dtype=cp.uint64)
        kernel(
            ((len(batch) + 255) // 256,),
            (256,),
            (
                batch_gpu,
                labels_gpu,
                eps_gpu,
                passing_gpu,
                batch_histogram_gpu,
                len(batch),
            ),
        )
        passing = cp.asnumpy(passing_gpu).astype(bool)
        batch_histogram = cp.asnumpy(batch_histogram_gpu)
        survivors.extend(batch[passing].astype(int).tolist())
        for index, value in enumerate((1, 3, 5)):
            histogram[value] += int(batch_histogram[index])
        checked += len(batch)
    properties = cp.cuda.runtime.getDeviceProperties(0)
    name = properties["name"]
    if isinstance(name, bytes):
        name = name.decode()
    return checked, survivors, histogram, str(name)


def run(backend: str, batch_size: int) -> dict:
    started = time.time()
    labels, epsilons = direction_tables()
    if backend == "cuda":
        checked, survivors, histogram, device = scan_cuda(
            labels, epsilons, batch_size
        )
    else:
        checked, survivors, histogram = scan_numpy(labels, epsilons, batch_size)
        device = "CPU/NumPy"
    expected = math.comb(49, 5)
    if checked != expected:
        raise AssertionError("combination enumeration was incomplete")
    canonical = json.dumps(survivors, separators=(",", ":")).encode()
    return {
        "experiment": "p7_size6_negative_infinity_floor_gpu",
        "status": "complete_exact_integer_floor_budget_sieve",
        "p": 7,
        "c_H": -1,
        "boundary_size": 6,
        "infinity_value": 1,
        "finite_boundary_size": 5,
        "backend": backend,
        "device": device,
        "batch_size": batch_size,
        "criterion": {
            "phase": 1,
            "scaled_floor_by_odd_fibre_count": {"1": 6, "3": 14, "5": 6},
            "directions_per_type": 4,
            "budget_per_type": 32,
            "maximum_b3_directions_per_type": 1,
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
    parser.add_argument("--batch-size", type=int, default=262144)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = run(args.backend, args.batch_size)
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
