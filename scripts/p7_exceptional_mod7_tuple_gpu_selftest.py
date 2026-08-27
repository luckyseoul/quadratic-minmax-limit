#!/usr/bin/env python3
"""Independent small brute-force audit of the mod-7 tuple CUDA join."""
from __future__ import annotations

import itertools
import argparse
import json
import os
import time
from pathlib import Path

import cupy as cp
import numpy as np

from p7_exceptional_projected_join_gpu import CUDA_SOURCE


def pack(rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    words = []
    for half in (rows[:, :11], rows[:, 11:]):
        word = np.zeros(len(rows), dtype=np.uint64)
        for index in range(half.shape[1]):
            word |= half[:, index].astype(np.uint64) << (4 * index)
        words.append(word)
    return words[0], words[1]


def key(digits: np.ndarray) -> int:
    value = 0
    place = 1
    for digit in digits:
        value += place * int(digit)
        place *= 7
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = time.time()
    rng = np.random.default_rng(15661)
    sizes = (3, 4, 2, 5)
    projections = [
        [rng.integers(0, 7, size=(size, 22), dtype=np.uint8) for size in sizes]
        for _ in range(3)
    ]
    planted = (1, 2, 0, 3)
    targets = [
        sum((projection[index][planted[index]].astype(np.int16) for index in range(4)),
            start=np.zeros(22, dtype=np.int16)) % 7
        for projection in projections
    ]
    cpu_build = set()
    for indices in itertools.product(range(sizes[0]), range(sizes[1])):
        cpu_build.add(tuple(
            key((projection[0][indices[0]] + projection[1][indices[1]]) % 7)
            for projection in projections
        ))
    cpu_matches = 0
    for indices in itertools.product(range(sizes[2]), range(sizes[3])):
        needed = tuple(
            key((target - projection[2][indices[0]] - projection[3][indices[1]]) % 7)
            for target, projection in zip(targets, projections)
        )
        cpu_matches += int(needed in cpu_build)

    module = cp.RawModule(code=CUDA_SOURCE, options=("--std=c++14",))
    make_sum = module.get_function("make_sum_keys7_22")
    count = module.get_function("count_needed_sorted_triples7_22")
    gpu = [[tuple(cp.asarray(word) for word in pack(rows)) for rows in projection] for projection in projections]
    zero = (cp.asarray([0], dtype=cp.uint64), cp.asarray([0], dtype=cp.uint64))
    build_count = sizes[0] * sizes[1]
    build = []
    for projection in gpu:
        out = cp.empty(build_count, dtype=cp.uint64)
        make_sum((1,), (256,), (
            projection[0][0], projection[0][1], np.uint64(sizes[0]),
            projection[1][0], projection[1][1], np.uint64(sizes[1]),
            zero[0], zero[1], np.uint64(1),
            np.uint64(0), np.uint64(build_count), out,
        ))
        build.append(out)
    order = cp.lexsort(cp.stack((build[2], build[1], build[0])))
    sorted_keys = [row[order] for row in build]
    packed_targets = [pack(target.reshape(1, 22)) for target in targets]
    match_count = cp.zeros(1, dtype=cp.uint64)
    probe_count = sizes[2] * sizes[3]
    count((1,), (256,), (
        gpu[0][2][0], gpu[0][2][1], gpu[1][2][0], gpu[1][2][1], gpu[2][2][0], gpu[2][2][1], np.uint64(sizes[2]),
        gpu[0][3][0], gpu[0][3][1], gpu[1][3][0], gpu[1][3][1], gpu[2][3][0], gpu[2][3][1], np.uint64(sizes[3]),
        zero[0], zero[1], zero[0], zero[1], zero[0], zero[1], np.uint64(1),
        np.uint64(packed_targets[0][0][0]), np.uint64(packed_targets[0][1][0]),
        np.uint64(packed_targets[1][0][0]), np.uint64(packed_targets[1][1][0]),
        np.uint64(packed_targets[2][0][0]), np.uint64(packed_targets[2][1][0]),
        np.uint64(0), np.uint64(probe_count),
        sorted_keys[0], sorted_keys[1], sorted_keys[2], np.uint64(build_count), match_count,
    ))
    gpu_matches = int(match_count.get()[0])
    if gpu_matches != cpu_matches:
        raise AssertionError((cpu_matches, gpu_matches))
    if cpu_matches < 1:
        raise AssertionError("planted tuple was lost")
    result = {
        "experiment": "p7_exceptional_mod7_tuple_gpu_selftest",
        "status": "passed_independent_cpu_gpu_exact_match_count",
        "seed": 15661,
        "catalog_sizes": list(sizes),
        "planted_catalog_indices": list(planted),
        "cpu_matches": cpu_matches,
        "gpu_matches": gpu_matches,
        "gpu": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        "elapsed_seconds": time.time() - started,
    }
    if args.output is not None:
        temporary = args.output.with_name(f".{args.output.name}.{os.getpid()}.tmp")
        temporary.write_text(json.dumps(result, indent=2) + "\n")
        os.replace(temporary, args.output)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
