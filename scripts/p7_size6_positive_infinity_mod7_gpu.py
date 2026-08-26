#!/usr/bin/env python3
"""Exact mod-seven filter for p=7 positive infinity-plus-five boundaries.

For ``c_H=+1`` and infinity in the six-point boundary, every direction has
phase zero and scaled mean eight.  The complete ``J(7,4)`` classification is
unique: ``A=t mod 2`` for odd-fibre size one or five, and
``A=(t-2)^2`` for odd-fibre size three.  Thus a finite five-set determines
all 280 affine score right sides.  This program checks all C(49,5) sets
against the 135 left-null dependencies of the common 282-by-1225 edge
system over F_7.  Rejection is an exact affine infeasibility certificate.
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
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)
from p7_unsaturated_modular_catalog_filter import (  # noqa: E402
    equation_matrix,
    left_dependencies,
)


POINTS = tuple(itertools.combinations(range(7), 4))
POINT_MASKS = tuple(sum(1 << value for value in point) for point in POINTS)


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def dependency_tables() -> tuple[np.ndarray, np.ndarray, dict]:
    matrix = equation_matrix()
    rank, dependencies = left_dependencies(matrix, 7)
    if rank != 147 or dependencies.shape != (135, 282):
        raise AssertionError("unexpected mod-seven dimensions")
    if np.any(dependencies @ (matrix.astype(np.int64) % 7) % 7):
        raise AssertionError("left-null dependency audit failed")

    base = (
        dependencies[:, :2] @ np.asarray([29, 1], dtype=np.int64) % 7
    ).astype(np.int16)
    tables = np.full((8, 128, 135), -1, dtype=np.int16)
    valid_masks = tuple(
        mask for mask in range(128) if mask.bit_count() in (1, 3, 5)
    )
    for direction_index in range(8):
        block = dependencies[
            :, 2 + 35 * direction_index : 2 + 35 * (direction_index + 1)
        ].astype(np.int64)
        for mask in valid_masks:
            b = mask.bit_count()
            slack = []
            for point_mask in POINT_MASKS:
                t = (mask & point_mask).bit_count()
                slack.append((t - 2) ** 2 if b == 3 else (t & 1))
            bad = 13 - np.asarray(slack, dtype=np.int64)
            tables[direction_index, mask] = (block @ bad % 7).astype(
                np.int16
            )
    return base, tables, {
        "equations": int(matrix.shape[0]),
        "edge_variables": int(matrix.shape[1]),
        "rank_mod_7": rank,
        "left_dependency_dimension": int(len(dependencies)),
        "left_null_audit": True,
        "valid_odd_fibre_masks": len(valid_masks),
    }


def finite_labels() -> np.ndarray:
    rows = []
    for direction in projective_directions(7):
        _eps, labels = field_direction_data(7, direction)
        rows.append(labels)
    out = np.asarray(rows, dtype=np.int16)
    if out.shape != (8, 49):
        raise AssertionError("unexpected affine label table")
    return out


def batches(size: int):
    source = itertools.combinations(range(49), 5)
    while True:
        rows = list(itertools.islice(source, size))
        if not rows:
            return
        yield np.asarray(rows, dtype=np.int16)


def boundary_masks_numpy(batch: np.ndarray, labels: np.ndarray) -> np.ndarray:
    out = np.empty((len(batch), 8), dtype=np.int16)
    for direction in range(8):
        fibre_labels = labels[direction, batch]
        out[:, direction] = np.bitwise_xor.reduce(
            np.left_shift(1, fibre_labels), axis=1
        )
    return out


def scan_numpy(
    base: np.ndarray, tables: np.ndarray, labels: np.ndarray, batch_size: int
) -> tuple[int, list[list[int]], dict[int, int]]:
    checked = 0
    survivors: list[list[int]] = []
    mask_size_histogram = {1: 0, 3: 0, 5: 0}
    for batch in batches(batch_size):
        masks = boundary_masks_numpy(batch, labels)
        syndrome = np.broadcast_to(base, (len(batch), 135)).copy()
        for direction in range(8):
            selected = tables[direction, masks[:, direction]]
            if np.any(selected < 0):
                raise AssertionError("boundary produced a non-odd fibre mask")
            syndrome += selected
            sizes = np.fromiter(
                (int(mask).bit_count() for mask in masks[:, direction]),
                dtype=np.int8,
                count=len(batch),
            )
            for size in mask_size_histogram:
                mask_size_histogram[size] += int(np.count_nonzero(sizes == size))
        passing = np.all(syndrome % 7 == 0, axis=1)
        survivors.extend(batch[passing].astype(int).tolist())
        checked += len(batch)
        print(
            f"checked={checked} survivors={len(survivors)}",
            flush=True,
        )
    return checked, survivors, mask_size_histogram


def scan_cuda(
    base: np.ndarray, tables: np.ndarray, labels: np.ndarray, batch_size: int
) -> tuple[int, list[list[int]], dict[int, int], str]:
    import cupy as cp

    if cp.cuda.runtime.getDeviceCount() < 1:
        raise RuntimeError("CUDA is unavailable")
    base_t = cp.asarray(base, dtype=cp.int16)
    tables_t = cp.asarray(tables, dtype=cp.int16)
    labels_t = cp.asarray(labels, dtype=cp.int16)
    kernel = cp.RawKernel(
        r'''
        extern "C" __global__ void filter_boundaries(
            const short* combinations,
            const short* labels,
            const short* base,
            const short* tables,
            unsigned char* passing,
            unsigned long long* histogram,
            const int count) {
          const int row = blockIdx.x;
          const int coordinate = threadIdx.x;
          if (row >= count) return;
          __shared__ int masks[8];
          if (coordinate < 8) {
            int mask = 0;
            #pragma unroll
            for (int column = 0; column < 5; ++column) {
              const int point = combinations[5 * row + column];
              mask ^= 1 << labels[49 * coordinate + point];
            }
            masks[coordinate] = mask;
            const int size = __popc((unsigned int)mask);
            if (size == 1) atomicAdd(histogram + 0, 1ULL);
            else if (size == 3) atomicAdd(histogram + 1, 1ULL);
            else if (size == 5) atomicAdd(histogram + 2, 1ULL);
            else passing[row] = 0;
          }
          if (coordinate == 0) passing[row] = 1;
          __syncthreads();
          if (coordinate < 135) {
            int syndrome = base[coordinate];
            #pragma unroll
            for (int direction = 0; direction < 8; ++direction) {
              const int mask = masks[direction];
              const short value = tables[
                  (direction * 128 + mask) * 135 + coordinate];
              if (value < 0) passing[row] = 0;
              syndrome += value;
            }
            if (syndrome % 7 != 0) passing[row] = 0;
          }
        }
        ''',
        "filter_boundaries",
        options=("--std=c++11",),
    )
    checked = 0
    survivors: list[list[int]] = []
    mask_size_histogram = {1: 0, 3: 0, 5: 0}
    for batch in batches(batch_size):
        batch_t = cp.asarray(batch, dtype=cp.int16)
        passing = cp.empty(len(batch), dtype=cp.uint8)
        histogram = cp.zeros(3, dtype=cp.uint64)
        kernel(
            (len(batch),),
            (160,),
            (batch_t, labels_t, base_t, tables_t, passing, histogram, len(batch)),
        )
        passing_cpu = cp.asnumpy(passing).astype(bool)
        batch_histogram = cp.asnumpy(histogram)
        for index, size in enumerate((1, 3, 5)):
            mask_size_histogram[size] += int(batch_histogram[index])
        survivors.extend(batch[passing_cpu].astype(int).tolist())
        checked += len(batch)
        print(
            f"checked={checked} survivors={len(survivors)}",
            flush=True,
        )
    return (
        checked,
        survivors,
        mask_size_histogram,
        cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
    )


def run(backend: str, batch_size: int) -> dict:
    started = time.time()
    base, tables, linear = dependency_tables()
    labels = finite_labels()
    if backend == "cuda":
        checked, survivors, histogram, device = scan_cuda(
            base, tables, labels, batch_size
        )
    else:
        checked, survivors, histogram = scan_numpy(
            base, tables, labels, batch_size
        )
        device = "CPU/NumPy"
    expected = math.comb(49, 5)
    if checked != expected:
        raise AssertionError("combination enumeration was incomplete")
    survivor_bytes = json.dumps(survivors, separators=(",", ":")).encode()
    return {
        "experiment": "p7_size6_positive_infinity_mod7_gpu",
        "status": "complete_exact_mod_seven_boundary_exhaustion",
        "p": 7,
        "c_H": 1,
        "boundary_size": 6,
        "infinity_value": 1,
        "finite_boundary_size": 5,
        "backend": backend,
        "device": device,
        "batch_size": batch_size,
        "linear_system": linear,
        "unique_slack_classification": {
            "b1_b5": "A(X)=|X cap B| mod 2",
            "b3": "A(X)=(|X cap B|-2)^2",
            "scaled_mean": 8,
        },
        "all_boundaries": expected,
        "checked_boundaries": checked,
        "direction_mask_size_histogram": {
            str(key): value for key, value in histogram.items()
        },
        "surviving_boundaries": len(survivors),
        "mod7_infeasible_boundaries": checked - len(survivors),
        "all_mod7_infeasible": not survivors,
        "survivor_sha256": hashlib.sha256(survivor_bytes).hexdigest(),
        "survivors_finite_field": survivors,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", choices=("numpy", "cuda"), default="cuda")
    parser.add_argument("--batch-size", type=int, default=65536)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = run(args.backend, args.batch_size)
    atomic_write(args.output, out)
    print(json.dumps({key: value for key, value in out.items() if key != "survivors_finite_field"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
