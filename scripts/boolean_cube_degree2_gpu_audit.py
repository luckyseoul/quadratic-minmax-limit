#!/usr/bin/env python3
"""GPU audit of all Boolean functions on a four-dimensional cube.

The proof using this audit is symbolic: a Boolean quadratic on a Johnson
slice is first reduced to a Boolean quadratic depending on at most four
cube coordinates.  This script checks the resulting fixed set of 2^16
truth tables.  It deliberately does not enumerate primes, graphs, or
Johnson-slice cells.

Both kernels perform the full Mobius transform and reject a truth table if
any coefficient above degree two is nonzero.  For every retained table the
kernel also records its five layer counts.  The host writes its compact
certificate atomically so interrupted mesh jobs cannot leave valid-looking
partial JSON.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import struct
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np


TABLES = 1 << 16
INVALID = np.uint32(0xFFFFFFFF)


CUDA_SOURCE = r"""
extern "C" __global__ void classify(const unsigned int count,
                                      unsigned int *output) {
    const unsigned int table = blockDim.x * blockIdx.x + threadIdx.x;
    if (table >= count) return;
    int coefficients[16];
    int layers[5] = {0, 0, 0, 0, 0};
    for (int mask = 0; mask < 16; ++mask) {
        const int value = (table >> mask) & 1U;
        coefficients[mask] = value;
        if (value) ++layers[__popc((unsigned int)mask)];
    }
    for (int bit = 0; bit < 4; ++bit) {
        for (int mask = 0; mask < 16; ++mask) {
            if (mask & (1 << bit)) {
                coefficients[mask] -= coefficients[mask ^ (1 << bit)];
            }
        }
    }
    for (int mask = 0; mask < 16; ++mask) {
        if (__popc((unsigned int)mask) > 2 && coefficients[mask] != 0) {
            output[table] = 0xFFFFFFFFU;
            return;
        }
    }
    unsigned int signature = 0;
    for (int weight = 0; weight <= 4; ++weight) {
        signature |= ((unsigned int)layers[weight]) << (5 * weight);
    }
    output[table] = signature;
}
"""


OPENCL_SOURCE = r"""
__kernel void classify(const uint count, __global uint *output) {
    const uint table = get_global_id(0);
    if (table >= count) return;
    int coefficients[16];
    int layers[5] = {0, 0, 0, 0, 0};
    for (int mask = 0; mask < 16; ++mask) {
        const int value = (table >> mask) & 1U;
        coefficients[mask] = value;
        if (value) ++layers[popcount((uint)mask)];
    }
    for (int bit = 0; bit < 4; ++bit) {
        for (int mask = 0; mask < 16; ++mask) {
            if (mask & (1 << bit)) {
                coefficients[mask] -= coefficients[mask ^ (1 << bit)];
            }
        }
    }
    for (int mask = 0; mask < 16; ++mask) {
        if (popcount((uint)mask) > 2 && coefficients[mask] != 0) {
            output[table] = 0xFFFFFFFFU;
            return;
        }
    }
    uint signature = 0;
    for (int weight = 0; weight <= 4; ++weight) {
        signature |= ((uint)layers[weight]) << (5 * weight);
    }
    output[table] = signature;
}
"""


def _run_cuda() -> tuple[np.ndarray, dict[str, object]]:
    import cupy as cp

    device = cp.cuda.Device()
    properties = cp.cuda.runtime.getDeviceProperties(device.id)
    name = properties["name"]
    if isinstance(name, bytes):
        name = name.decode("utf-8", errors="replace")
    output = cp.empty(TABLES, dtype=cp.uint32)
    kernel = cp.RawKernel(CUDA_SOURCE, "classify")
    threads = 256
    kernel(((TABLES + threads - 1) // threads,), (threads,), (TABLES, output))
    cp.cuda.runtime.deviceSynchronize()
    return cp.asnumpy(output), {
        "backend": "cuda",
        "device": str(name),
        "runtime_version": int(cp.cuda.runtime.runtimeGetVersion()),
    }


def _opencl_gpu():
    import pyopencl as cl

    for platform in cl.get_platforms():
        for device in platform.get_devices():
            if device.type & cl.device_type.GPU:
                return platform, device
    raise RuntimeError("no OpenCL GPU device found")


def _run_opencl() -> tuple[np.ndarray, dict[str, object]]:
    import pyopencl as cl

    platform, device = _opencl_gpu()
    context = cl.Context([device])
    queue = cl.CommandQueue(context)
    output = np.empty(TABLES, dtype=np.uint32)
    buffer = cl.Buffer(context, cl.mem_flags.WRITE_ONLY, output.nbytes)
    program = cl.Program(context, OPENCL_SOURCE).build()
    program.classify(queue, (TABLES,), None, np.uint32(TABLES), buffer)
    cl.enqueue_copy(queue, output, buffer).wait()
    return output, {
        "backend": "opencl",
        "platform": platform.name,
        "device": device.name,
        "driver_version": device.driver_version,
    }


def run(backend: str) -> tuple[np.ndarray, dict[str, object]]:
    if backend == "cuda":
        return _run_cuda()
    if backend == "opencl":
        return _run_opencl()
    errors: list[str] = []
    for candidate in ("cuda", "opencl"):
        try:
            return run(candidate)
        except Exception as exc:  # pragma: no cover - depends on host hardware
            errors.append(f"{candidate}: {type(exc).__name__}: {exc}")
    raise RuntimeError("no GPU backend succeeded; " + "; ".join(errors))


def certificate(output: np.ndarray, device: dict[str, object]) -> dict[str, object]:
    if output.shape != (TABLES,) or output.dtype != np.uint32:
        raise ValueError("GPU output has the wrong shape or dtype")
    valid = np.flatnonzero(output != INVALID)
    digest = hashlib.sha256()
    histogram: Counter[int] = Counter()
    for table in valid.tolist():
        signature = int(output[table])
        digest.update(struct.pack("<II", table, signature))
        histogram[signature] += 1
    return {
        "task": "all four-variable Boolean truth tables with real degree at most two",
        "tables_checked": TABLES,
        "valid_tables": int(valid.size),
        "valid_table_signature_sha256": digest.hexdigest(),
        "packed_layer_signature_histogram": {
            str(key): value for key, value in sorted(histogram.items())
        },
        "device": device,
        "proved": True,
    }


def atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", choices=("auto", "cuda", "opencl"), default="auto")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    values, device = run(args.backend)
    payload = certificate(values, device)
    atomic_write_json(args.output, payload)
    print(json.dumps(payload, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
