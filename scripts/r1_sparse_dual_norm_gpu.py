#!/usr/bin/env python3
"""GPU reconnaissance for sparse representatives of the R1 dual lattice.

For ``L* = P Z^n`` with ``P=(I+C/p)/2`` and a sparse vector
``z in {0,+-1}^n``, the exact scaled norm is

    2p ||Pz||^2 = p ||z||^2 + z^T C z.

The CUDA kernel assigns one support to each logical thread, scans all sign
patterns modulo the global ``z -> -z`` symmetry, and accumulates only a
block-local norm histogram.  The host receives O(pk+k^2) counters rather
than billions of representatives.

This is a discovery tool, not a complete lattice enumerator: coefficients
outside ``{0,+-1}`` are not scanned.  It is nevertheless collision-free when
``2k < p+1``.  Indeed, if ``Pz=Pz'`` and ``h=z-z'`` is nonzero, then
``Ch=-ph``.  At a coordinate where ``|h_i|=M=max_j |h_j|``, the zero diagonal
and unit off-diagonal entries of ``C`` give

    p M <= sum_{j != i} |h_j| = ||h||_1-M,

so ``||h||_1 >= (p+1)M >= p+1``.  But two scanned support-``k`` vectors have
``||z-z'||_1 <= 2k``, a contradiction.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


KERNEL = r'''
extern "C" __global__
void sparse_norm_hist(
    const signed char* conference,
    const unsigned long long* choose,
    const int choose_stride,
    const int n,
    const int p,
    const int support,
    const unsigned long long combinations,
    unsigned long long* histogram,
    const int histogram_size)
{
    extern __shared__ unsigned int local_histogram[];
    for (int bin = threadIdx.x; bin < histogram_size; bin += blockDim.x)
        local_histogram[bin] = 0U;
    __syncthreads();

    const unsigned long long logical_thread =
        (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned long long logical_threads =
        (unsigned long long)gridDim.x * blockDim.x;

    for (unsigned long long rank0 = logical_thread;
         rank0 < combinations;
         rank0 += logical_threads)
    {
        int index[12];
        unsigned long long rank = rank0;
        int next = 0;
        bool valid = true;

        // Lexicographic unranking of a support from [0,n).
        for (int position = 0; position < support; ++position)
        {
            const int remaining = support - position - 1;
            bool selected = false;
            const int last = n - (remaining + 1);
            for (int candidate = next; candidate <= last; ++candidate)
            {
                const unsigned long long ways =
                    choose[(n - candidate - 1) * choose_stride + remaining];
                if (rank < ways)
                {
                    index[position] = candidate;
                    next = candidate + 1;
                    selected = true;
                    break;
                }
                rank -= ways;
            }
            if (!selected) { valid = false; break; }
        }
        if (!valid) continue;

        // Fix the first sign to +1, quotienting the global +/- symmetry.
        const unsigned int sign_patterns = 1U << (support - 1);
        for (unsigned int mask = 0; mask < sign_patterns; ++mask)
        {
            int cross = 0;
            for (int a = 0; a < support; ++a)
            {
                const int sign_a =
                    (a == 0 || ((mask >> (a - 1)) & 1U) == 0U) ? 1 : -1;
                for (int b = a + 1; b < support; ++b)
                {
                    const int sign_b =
                        (((mask >> (b - 1)) & 1U) == 0U) ? 1 : -1;
                    cross += sign_a * sign_b
                        * (int)conference[index[a] * n + index[b]];
                }
            }
            const int scaled_norm = p * support + 2 * cross;
            if (scaled_norm >= 0 && scaled_norm < histogram_size)
                atomicAdd(&local_histogram[scaled_norm], 1U);
        }
    }
    __syncthreads();

    for (int bin = threadIdx.x; bin < histogram_size; bin += blockDim.x)
    {
        const unsigned int value = local_histogram[bin];
        if (value)
            atomicAdd(&histogram[bin], (unsigned long long)value);
    }
}
'''


def conference_matrix(p: int) -> np.ndarray:
    matrix = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    n = p * p + 1
    if matrix.shape != (n, n):
        raise ArithmeticError(f"unexpected conference shape {matrix.shape}")
    if not np.array_equal(matrix, matrix.T):
        raise ArithmeticError("R1 scanner requires a symmetric conference matrix")
    if np.any(np.diag(matrix)) or not np.all(np.isin(matrix, (-1, 0, 1))):
        raise ArithmeticError("conference matrix is not zero-diagonal +/-1")
    check = matrix.astype(np.int64) @ matrix.astype(np.int64)
    if not np.array_equal(check, p * p * np.eye(n, dtype=np.int64)):
        raise ArithmeticError("conference identity C^2=p^2 I failed")
    return np.ascontiguousarray(matrix)


def choose_table(n: int, support: int) -> np.ndarray:
    table = np.zeros((n + 1, support + 1), dtype=np.uint64)
    for a in range(n + 1):
        for b in range(min(a, support) + 1):
            value = math.comb(a, b)
            if value >= 1 << 64:
                raise OverflowError(f"binomial C({a},{b}) does not fit uint64")
            table[a, b] = value
    return np.ascontiguousarray(table)


def cpu_histogram(matrix: np.ndarray, p: int, support: int) -> Counter[int]:
    out: Counter[int] = Counter()
    for indices in itertools.combinations(range(len(matrix)), support):
        pairs = [
            (a, b, int(matrix[indices[a], indices[b]]))
            for a in range(support)
            for b in range(a + 1, support)
        ]
        for mask in range(1 << (support - 1)):
            signs = [1] + [
                -1 if (mask >> (position - 1)) & 1 else 1
                for position in range(1, support)
            ]
            cross = sum(2 * signs[a] * signs[b] * value for a, b, value in pairs)
            out[p * support + cross] += 1
    return out


def gpu_histogram(
    kernel,
    matrix_gpu,
    table_gpu,
    p: int,
    support: int,
    blocks: int,
    threads: int,
) -> tuple[np.ndarray, float]:
    import cupy as cp

    n = p * p + 1
    combinations = math.comb(n, support)
    max_norm = p * support + support * (support - 1)
    histogram_size = max_norm + 1
    histogram_gpu = cp.zeros(histogram_size, dtype=cp.uint64)
    shared_bytes = histogram_size * np.dtype(np.uint32).itemsize

    start = cp.cuda.Event()
    stop = cp.cuda.Event()
    start.record()
    kernel(
        (blocks,),
        (threads,),
        (
            matrix_gpu,
            table_gpu,
            np.int32(support + 1),
            np.int32(n),
            np.int32(p),
            np.int32(support),
            np.uint64(combinations),
            histogram_gpu,
            np.int32(histogram_size),
        ),
        shared_mem=shared_bytes,
    )
    stop.record()
    stop.synchronize()
    elapsed = float(cp.cuda.get_elapsed_time(start, stop)) / 1000.0
    return cp.asnumpy(histogram_gpu), elapsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--supports", default="1,2,3,4")
    parser.add_argument("--blocks", type=int, default=4096)
    parser.add_argument("--threads", type=int, default=256)
    parser.add_argument(
        "--cpu-check-limit",
        type=int,
        default=2_000_000,
        help="CPU-check scans at most this many signed representatives",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    supports = [int(value) for value in args.supports.split(",")]
    if not supports or min(supports) < 1 or max(supports) > 12:
        raise ValueError("supports must lie in 1..12")
    p = args.p
    n = p * p + 1

    import cupy as cp

    matrix = conference_matrix(p)
    matrix_gpu = cp.asarray(matrix)
    kernel_start = time.perf_counter()
    kernel = cp.RawKernel(KERNEL, "sparse_norm_hist", options=("-std=c++11",))
    kernel.compile()
    compile_seconds = time.perf_counter() - kernel_start

    rows = []
    for support in supports:
        combinations = math.comb(n, support)
        representatives = combinations * (1 << (support - 1))
        table = choose_table(n, support)
        histogram, elapsed = gpu_histogram(
            kernel,
            matrix_gpu,
            cp.asarray(table),
            p,
            support,
            args.blocks,
            args.threads,
        )
        observed = {
            str(index): int(count)
            for index, count in enumerate(histogram)
            if count
        }
        if int(histogram.sum(dtype=np.uint64)) != representatives:
            raise ArithmeticError(
                f"GPU count mismatch at support {support}: "
                f"{int(histogram.sum())} != {representatives}"
            )

        cpu_checked = representatives <= args.cpu_check_limit
        if cpu_checked:
            control = cpu_histogram(matrix, p, support)
            if observed != {str(key): value for key, value in sorted(control.items())}:
                raise ArithmeticError(f"CPU/GPU histogram mismatch at k={support}")

        rows.append(
            {
                "support": support,
                "combinations": combinations,
                "half_representatives": representatives,
                "collision_free_by_2k_lt_p_plus_1": 2 * support < p + 1,
                "elapsed_seconds": elapsed,
                "representatives_per_second": representatives / elapsed,
                "cpu_checked": cpu_checked,
                "scaled_norm_histogram": observed,
            }
        )
        print(
            f"p={p} k={support} reps={representatives} time={elapsed:.6f}s "
            f"rate={representatives / elapsed:.3e}/s norms={observed}",
            flush=True,
        )

    result = {
        "experiment": "r1_sparse_dual_norm_gpu",
        "status": "finite_reconnaissance_not_proof",
        "p": p,
        "n": n,
        "device": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        "compute_capability": cp.cuda.Device().compute_capability,
        "cupy": cp.__version__,
        "kernel_compile_seconds": compile_seconds,
        "blocks": args.blocks,
        "threads": args.threads,
        "rows": rows,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
