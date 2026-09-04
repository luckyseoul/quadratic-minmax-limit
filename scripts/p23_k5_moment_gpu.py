#!/usr/bin/env python3
"""Cross-backend GPU replay for the p=23 K5 quartic-octic endpoint."""
from __future__ import annotations

import argparse
import json
import socket
import time


P = 23
TOTAL = P**5

KERNEL_BODY = r"""
inline uint mod23(long value) {
    long reduced = value % 23;
    return (uint)(reduced < 0 ? reduced + 23 : reduced);
}

inline void audit_index(ulong gid, volatile __global uint *counts) {
    uint x4 = (uint)(gid % 23); gid /= 23;
    uint x3 = (uint)(gid % 23); gid /= 23;
    uint x2 = (uint)(gid % 23); gid /= 23;
    uint x1 = (uint)(gid % 23); gid /= 23;
    uint x0 = (uint)(gid % 23);
    if (!(x0 < x1 && x1 < x2 && x2 < x3 && x3 < x4)) return;
    uint xs[5] = {x0, x1, x2, x3, x4};
    uint s2 = 0, s4 = 0, s6 = 0, s8 = 0;
    for (uint left = 0; left < 5; ++left) {
        for (uint right = left + 1; right < 5; ++right) {
            uint delta = mod23((long)xs[left] - (long)xs[right]);
            uint d2 = mod23((long)delta * delta);
            uint d4 = mod23((long)d2 * d2);
            uint d6 = mod23((long)d4 * d2);
            uint d8 = mod23((long)d4 * d4);
            s2 = mod23((long)s2 + d2);
            s4 = mod23((long)s4 + d4);
            s6 = mod23((long)s6 + d6);
            s8 = mod23((long)s8 + d8);
        }
    }
    uint g4 = mod23(-2L * s4 - (long)s2 * s2);
    uint s2sq = mod23((long)s2 * s2);
    uint s2four = mod23((long)s2sq * s2sq);
    uint g8 = mod23(-24L * s8 - 32L * s2 * s6 + 5L * s2four);
    atomic_inc(counts + 0);
    if (g4 == 0) atomic_inc(counts + 1);
    if (g8 == 0) atomic_inc(counts + 2);
    if (g4 == 0 && g8 == 0) atomic_inc(counts + 3);
}
"""


def run_opencl() -> tuple[str, list[int], float]:
    import numpy as np
    import pyopencl as cl

    devices = [device for platform in cl.get_platforms() for device in platform.get_devices(device_type=cl.device_type.GPU)]
    if not devices:
        raise RuntimeError("no OpenCL GPU")
    device = devices[0]
    context = cl.Context([device])
    queue = cl.CommandQueue(context)
    source = KERNEL_BODY + r"""
__kernel void audit(volatile __global uint *counts) {
    audit_index((ulong)get_global_id(0), counts);
}
"""
    program = cl.Program(context, source).build()
    host = np.zeros(4, dtype=np.uint32)
    buffer = cl.Buffer(context, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=host)
    start = time.monotonic()
    program.audit(queue, (TOTAL,), None, buffer)
    cl.enqueue_copy(queue, host, buffer).wait()
    return device.name.strip(), [int(value) for value in host], time.monotonic() - start


def run_cuda() -> tuple[str, list[int], float]:
    import cupy as cp

    source = KERNEL_BODY.replace("volatile __global uint *", "unsigned int *").replace("inline uint", "__device__ unsigned int").replace("ulong", "unsigned long long").replace("uint", "unsigned int").replace("atomic_inc", "atomicAdd")
    source = source.replace("atomicAdd(counts + 0);", "atomicAdd(counts + 0, 1);").replace("atomicAdd(counts + 1);", "atomicAdd(counts + 1, 1);").replace("atomicAdd(counts + 2);", "atomicAdd(counts + 2, 1);").replace("atomicAdd(counts + 3);", "atomicAdd(counts + 3, 1);")
    source += r"""
extern "C" __global__ void audit(unsigned int *counts, unsigned long long total) {
    unsigned long long gid = (unsigned long long)blockDim.x * blockIdx.x + threadIdx.x;
    if (gid < total) audit_index(gid, counts);
}
"""
    kernel = cp.RawKernel(source, "audit", options=("--std=c++11",))
    counts = cp.zeros(4, dtype=cp.uint32)
    threads = 256
    blocks = (TOTAL + threads - 1) // threads
    start = time.monotonic()
    kernel((blocks,), (threads,), (counts, TOTAL))
    cp.cuda.Stream.null.synchronize()
    device = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
    return device, [int(value) for value in cp.asnumpy(counts)], time.monotonic() - start


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("backend", choices=("opencl", "cuda"))
    args = parser.parse_args()
    if args.backend == "opencl":
        device, counts, elapsed = run_opencl()
    else:
        device, counts, elapsed = run_cuda()
    expected = [33649, 1518, 2024, 0]
    print(json.dumps({
        "host": socket.gethostname(),
        "backend": args.backend,
        "device": device,
        "ordered_base23_tuples_checked": TOTAL,
        "strict_five_sets": counts[0],
        "quartic_zero_sets": counts[1],
        "octic_zero_sets": counts[2],
        "simultaneous_zero_sets": counts[3],
        "expected": expected,
        "proved": counts == expected,
        "elapsed_seconds": elapsed,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
