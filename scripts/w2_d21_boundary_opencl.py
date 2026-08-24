#!/usr/bin/env python3
"""OpenCL cross-check and search for the d=21 boundary criterion.

This is the Intel/portable counterpart of ``w2_d21_boundary_gpu.py``.  It
uses the same exact integer kernel and host-side field conventions, but runs
through OpenCL so selected records can be checked by a genuinely independent
GPU compiler and runtime.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from io_atomic import write_json_atomic  # noqa: E402
import w2_d21_boundary_gpu as common  # noqa: E402


def opencl_source() -> str:
    """Translate the deliberately CUDA/OpenCL-common integer kernel."""
    source = common.KERNEL_SOURCE
    source = source.replace("__device__ __forceinline__", "inline")
    source = source.replace('extern "C" __global__ void', "__kernel void")
    source = source.replace(
        "(unsigned long long)blockIdx.x * blockDim.x + threadIdx.x",
        "(unsigned long long)get_global_id(0)",
    )
    source = source.replace(
        "(unsigned long long)gridDim.x * blockDim.x",
        "(unsigned long long)get_global_size(0)",
    )
    source = source.replace(
        "const unsigned int* levels,", "__global const unsigned int* levels,"
    )
    source = source.replace(
        "const unsigned long long* roots,",
        "__global const unsigned long long* roots,",
    )
    source = source.replace("unsigned int* output)", "__global unsigned int* output)")
    source = source.replace("atomicXor", "atomic_xor")
    return source


def context_and_program():
    import pyopencl as cl

    devices = [
        device
        for platform in cl.get_platforms()
        for device in platform.get_devices()
        if device.type & cl.device_type.GPU
    ]
    if not devices:
        raise RuntimeError("no OpenCL GPU device found")
    device = devices[0]
    context = cl.Context([device])
    queue = cl.CommandQueue(context)
    # Intel IGC can otherwise over-inline the nested exact 64-bit modular
    # exponentiations and crash during optimization on Arc Alchemist.
    program = cl.Program(context, opencl_source()).build(options=["-cl-opt-disable"])
    return cl, device, context, queue, program


def record(p: int, runtime) -> dict:
    cl, _device, context, queue, program = runtime
    started = time.perf_counter()
    ib, sigma_character_inverse, roots = common.field_setup(p)
    pole_t = (-2 * pow(ib, p - 2, p)) % p
    line_levels = common.levels(p)
    output = np.zeros((12, 21), dtype=np.uint32)
    flags = cl.mem_flags
    level_buffer = cl.Buffer(
        context, flags.READ_ONLY | flags.COPY_HOST_PTR, hostbuf=line_levels
    )
    root_buffer = cl.Buffer(
        context, flags.READ_ONLY | flags.COPY_HOST_PTR, hostbuf=roots
    )
    output_buffer = cl.Buffer(
        context, flags.READ_WRITE | flags.COPY_HOST_PTR, hostbuf=output
    )
    local_size = 256
    useful = 12 * p
    global_size = min(65535 * local_size, ((useful + local_size - 1) // local_size) * local_size)
    program.d21_selected_lines(
        queue,
        (global_size,),
        (local_size,),
        np.uint32(p),
        np.uint32(ib),
        np.uint32(pole_t),
        np.uint64((p * p - 1) // 42),
        np.uint64(sigma_character_inverse),
        np.uint32(len(line_levels)),
        level_buffer,
        root_buffer,
        output_buffer,
    )
    cl.enqueue_copy(queue, output, output_buffer).wait()
    line_bins = (output & 1).astype(np.uint8)
    remainders = [[], []]
    polynomials = []
    for a in range(3):
        delta = np.bitwise_xor.reduce(line_bins[4 * a : 4 * a + 4], axis=0)
        polynomial = common.polynomial_bits(delta)
        polynomials.append(hex(polynomial))
        for factor_index, factor in enumerate(common.FACTORS):
            remainders[factor_index].append(common.remainder_bits(polynomial, factor))
    zero_triples = [all(value == 0 for value in side) for side in remainders]
    return {
        "p": p,
        "ib": ib,
        "pole_t": pole_t,
        "nonsquare_component_hex": polynomials,
        "remainders": remainders,
        "zero_triples": zero_triples,
        "simultaneous_zero_triples": all(zero_triples),
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--primes")
    parser.add_argument("--retain", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.primes:
        primes = sorted(set(map(int, args.primes.split(","))))
    else:
        if args.start is None or args.stop is None:
            raise ValueError("start and stop are required without explicit primes")
        primes = [
            int(p)
            for p in sp.primerange(max(29, args.start), args.stop + 1)
            if p % 84 == 29
        ]
    runtime = context_and_program()
    print(f"OpenCL device: {runtime[1].name}", flush=True)
    rows = []
    retained = []
    started = time.perf_counter()
    for index, p in enumerate(primes, 1):
        row = record(p, runtime)
        if args.retain or any(row["zero_triples"]):
            retained.append(row)
        rows.append(
            {
                "p": p,
                "zero_triples": row["zero_triples"],
                "simultaneous_zero_triples": row["simultaneous_zero_triples"],
                "elapsed_seconds": row["elapsed_seconds"],
            }
        )
        print(
            f"[{index}/{len(primes)}] p={p} zero={row['zero_triples']} "
            f"seconds={row['elapsed_seconds']:.4f}",
            flush=True,
        )
        if row["simultaneous_zero_triples"]:
            break
    result = {
        "backend": runtime[1].name,
        "range": None if args.primes else [args.start, args.stop],
        "explicit_primes": primes if args.primes else None,
        "congruence_class": "p == 29 (mod 84)",
        "n_primes_requested": len(primes),
        "n_primes_completed": len(rows),
        "one_sided_zero_triples": [row for row in retained if any(row["zero_triples"])],
        "counterexamples": [
            row for row in retained if row["simultaneous_zero_triples"]
        ],
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(json.dumps({key: result[key] for key in ("backend", "n_primes_completed", "elapsed_seconds")}), flush=True)


if __name__ == "__main__":
    main()
