#!/usr/bin/env python3
"""GPU signed degree-four harmonic sums on sparse R1 dual representatives.

This extends ``r1_sparse_dual_norm_gpu.py``.  For deterministic random
admissible matrices ``W`` (``PWP=W``, ``diag(W)=0``, ``||W||_F=1``), it
accumulates by exact scaled norm ``s=2p||Pz||^2`` the floating-point scout

    sum_z (-1)^s H_W(Pz/2).

One representative from each ``+/-z`` pair is scanned, so the kernel adds
``(-1)^s H_W(Pz)/8`` per representative (the factor is ``2/16``).  Counts
and harmonic sums use block-local atomics; only the small shell tables are
copied to the host.  This is finite sparse-shell reconnaissance, not a tail
bound or an R1 proof.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from r1_sparse_dual_norm_gpu import conference_matrix, choose_table  # noqa: E402


KERNEL = r'''
extern "C" __global__
void sparse_harmonic_hist(
    const signed char* conference,
    const double* matrices,
    const double* squares,
    const unsigned long long* choose,
    const int choose_stride,
    const int n,
    const int p,
    const int d,
    const int support,
    const int channels,
    const unsigned long long combinations,
    unsigned long long* histogram,
    double* harmonic,
    const int histogram_size,
    const int count_words)
{
    extern __shared__ unsigned char shared_bytes[];
    unsigned int* local_histogram = (unsigned int*)shared_bytes;
    double* local_harmonic = (double*)(local_histogram + count_words);

    for (int bin = threadIdx.x; bin < count_words; bin += blockDim.x)
        local_histogram[bin] = 0U;
    for (int item = threadIdx.x; item < channels * histogram_size;
         item += blockDim.x)
        local_harmonic[item] = 0.0;
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

        const unsigned int sign_patterns = 1U << (support - 1);
        for (unsigned int mask = 0; mask < sign_patterns; ++mask)
        {
            int sign[12];
            sign[0] = 1;
            for (int a = 1; a < support; ++a)
                sign[a] = ((mask >> (a - 1)) & 1U) ? -1 : 1;

            int cross = 0;
            for (int a = 0; a < support; ++a)
                for (int b = a + 1; b < support; ++b)
                    cross += sign[a] * sign[b]
                        * (int)conference[index[a] * n + index[b]];
            const int scaled_norm = p * support + 2 * cross;
            if (scaled_norm < 0 || scaled_norm >= histogram_size) continue;
            atomicAdd(&local_histogram[scaled_norm], 1U);

            const double norm = (double)scaled_norm / (2.0 * (double)p);
            const double phase_scale = (scaled_norm & 1) ? -0.125 : 0.125;
            for (int channel = 0; channel < channels; ++channel)
            {
                const long long base = (long long)channel * n * n;
                double quadratic = 0.0;
                double quadratic_square = 0.0;
                for (int a = 0; a < support; ++a)
                {
                    const int ia = index[a];
                    quadratic += matrices[base + (long long)ia * n + ia];
                    quadratic_square += squares[base + (long long)ia * n + ia];
                    for (int b = a + 1; b < support; ++b)
                    {
                        const int ib = index[b];
                        const double signed_twice = 2.0 * sign[a] * sign[b];
                        quadratic += signed_twice
                            * matrices[base + (long long)ia * n + ib];
                        quadratic_square += signed_twice
                            * squares[base + (long long)ia * n + ib];
                    }
                }
                const double H = quadratic * quadratic
                    - 4.0 * norm * quadratic_square / (double)(d + 4)
                    + 2.0 * norm * norm
                        / ((double)(d + 2) * (double)(d + 4));
                atomicAdd(
                    &local_harmonic[channel * histogram_size + scaled_norm],
                    phase_scale * H);
            }
        }
    }
    __syncthreads();

    for (int bin = threadIdx.x; bin < histogram_size; bin += blockDim.x)
    {
        const unsigned int value = local_histogram[bin];
        if (value)
            atomicAdd(&histogram[bin], (unsigned long long)value);
    }
    for (int item = threadIdx.x; item < channels * histogram_size;
         item += blockDim.x)
    {
        const double value = local_harmonic[item];
        if (value != 0.0) atomicAdd(&harmonic[item], value);
    }
}
'''


def symmetric_coordinate_map(eigenspace: np.ndarray) -> tuple[np.ndarray, list[tuple[int, int]]]:
    n, d = eigenspace.shape
    labels: list[tuple[int, int]] = []
    columns = []
    for a in range(d):
        labels.append((a, a))
        columns.append(eigenspace[:, a] ** 2)
        for b in range(a + 1, d):
            labels.append((a, b))
            columns.append(math.sqrt(2.0) * eigenspace[:, a] * eigenspace[:, b])
    return np.column_stack(columns).reshape(n, -1), labels


def random_admissible_matrices(
    conference: np.ndarray,
    p: int,
    channels: int,
    seed: int,
) -> tuple[np.ndarray, dict]:
    values, vectors = np.linalg.eigh(conference.astype(np.float64))
    eigenspace = vectors[:, np.isclose(values, p, atol=1e-8)]
    n, d = eigenspace.shape
    if d != n // 2:
        raise ArithmeticError(f"unexpected +p eigenspace dimension {d}")
    diagonal_map, labels = symmetric_coordinate_map(eigenspace)
    _left, singular, vh = np.linalg.svd(diagonal_map, full_matrices=False)
    tolerance = max(diagonal_map.shape) * singular[0] * np.finfo(float).eps
    rank = int(np.count_nonzero(singular > tolerance))
    row_basis = vh[:rank]
    rng = np.random.default_rng(seed)

    matrices = []
    diagnostics = []
    for channel in range(channels):
        coefficients = rng.standard_normal(len(labels))
        coefficients -= row_basis.T @ (row_basis @ coefficients)
        coefficients /= np.linalg.norm(coefficients)
        intrinsic = np.zeros((d, d), dtype=np.float64)
        for coefficient, (a, b) in zip(coefficients, labels, strict=True):
            if a == b:
                intrinsic[a, a] = coefficient
            else:
                intrinsic[a, b] = intrinsic[b, a] = coefficient / math.sqrt(2.0)
        ambient = eigenspace @ intrinsic @ eigenspace.T
        ambient = (ambient + ambient.T) / 2.0
        ambient /= np.linalg.norm(ambient)
        projection_error = np.linalg.norm(
            ambient
            - ((np.eye(n) + conference.astype(np.float64) / p) / 2.0)
            @ ambient
            @ ((np.eye(n) + conference.astype(np.float64) / p) / 2.0)
        )
        diagonal_error = float(np.max(np.abs(np.diag(ambient))))
        trace_error = abs(float(np.trace(ambient)))
        if max(projection_error, diagonal_error, trace_error) > 2e-12:
            raise ArithmeticError(
                f"admissibility failure channel={channel}: "
                f"{projection_error=}, {diagonal_error=}, {trace_error=}"
            )
        matrices.append(ambient)
        diagnostics.append(
            {
                "channel": channel,
                "projection_error": float(projection_error),
                "diagonal_error": diagonal_error,
                "trace_error": trace_error,
                "frobenius_norm": float(np.linalg.norm(ambient)),
            }
        )
    return np.ascontiguousarray(matrices), {
        "eigenspace_dimension": d,
        "diagonal_map_rank": rank,
        "admissible_dimension": len(labels) - rank,
        "channels": diagnostics,
    }


def scan_support(
    kernel,
    conference_gpu,
    matrices_gpu,
    squares_gpu,
    choose_gpu,
    p: int,
    support: int,
    channels: int,
    blocks: int,
    threads: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    import cupy as cp

    n = p * p + 1
    d = n // 2
    combinations = math.comb(n, support)
    histogram_size = p * support + support * (support - 1) + 1
    count_words = (histogram_size + 1) & ~1
    histogram_gpu = cp.zeros(histogram_size, dtype=cp.uint64)
    harmonic_gpu = cp.zeros((channels, histogram_size), dtype=cp.float64)
    shared_bytes = count_words * 4 + channels * histogram_size * 8
    start = cp.cuda.Event(); stop = cp.cuda.Event()
    start.record()
    kernel(
        (blocks,),
        (threads,),
        (
            conference_gpu,
            matrices_gpu,
            squares_gpu,
            choose_gpu,
            np.int32(support + 1),
            np.int32(n),
            np.int32(p),
            np.int32(d),
            np.int32(support),
            np.int32(channels),
            np.uint64(combinations),
            histogram_gpu,
            harmonic_gpu,
            np.int32(histogram_size),
            np.int32(count_words),
        ),
        shared_mem=shared_bytes,
    )
    stop.record(); stop.synchronize()
    elapsed = float(cp.cuda.get_elapsed_time(start, stop)) / 1000.0
    return cp.asnumpy(histogram_gpu), cp.asnumpy(harmonic_gpu), elapsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--supports", default="1,2,3,4,5")
    parser.add_argument("--channels", type=int, default=4)
    parser.add_argument("--seed", type=int, default=15631)
    parser.add_argument("--blocks", type=int, default=4096)
    parser.add_argument("--threads", type=int, default=256)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    supports = [int(value) for value in args.supports.split(",")]
    if min(supports) < 1 or max(supports) > 12:
        raise ValueError("supports must lie in 1..12")
    if not 1 <= args.channels <= 8:
        raise ValueError("channels must lie in 1..8")

    import cupy as cp

    p = args.p
    n = p * p + 1
    conference = conference_matrix(p)
    matrices, diagnostics = random_admissible_matrices(
        conference, p, args.channels, args.seed
    )
    squares = np.ascontiguousarray(np.matmul(matrices, matrices))
    conference_gpu = cp.asarray(conference)
    matrices_gpu = cp.asarray(matrices)
    squares_gpu = cp.asarray(squares)
    compile_start = time.perf_counter()
    kernel = cp.RawKernel(KERNEL, "sparse_harmonic_hist", options=("-std=c++11",))
    kernel.compile()
    compile_seconds = time.perf_counter() - compile_start

    aggregate_counts: dict[int, int] = {}
    aggregate_harmonic: dict[int, np.ndarray] = {}
    rows = []
    for support in supports:
        representatives = math.comb(n, support) * (1 << (support - 1))
        histogram, harmonic, elapsed = scan_support(
            kernel,
            conference_gpu,
            matrices_gpu,
            squares_gpu,
            cp.asarray(choose_table(n, support)),
            p,
            support,
            args.channels,
            args.blocks,
            args.threads,
        )
        if int(histogram.sum(dtype=np.uint64)) != representatives:
            raise ArithmeticError(f"count mismatch at support {support}")
        counts = {str(s): int(value) for s, value in enumerate(histogram) if value}
        sums = {
            str(s): [float(value) for value in harmonic[:, s]]
            for s, count in enumerate(histogram)
            if count
        }
        for s, count in enumerate(histogram):
            if not count:
                continue
            aggregate_counts[s] = aggregate_counts.get(s, 0) + int(count)
            aggregate_harmonic[s] = aggregate_harmonic.get(
                s, np.zeros(args.channels, dtype=np.float64)
            ) + harmonic[:, s]
        rows.append(
            {
                "support": support,
                "half_representatives": representatives,
                "collision_free_by_2k_lt_p_plus_1": 2 * support < p + 1,
                "elapsed_seconds": elapsed,
                "representatives_per_second": representatives / elapsed,
                "scaled_norm_histogram": counts,
                "phased_harmonic_sums": sums,
            }
        )
        print(
            f"p={p} k={support} reps={representatives} time={elapsed:.6f}s "
            f"rate={representatives / elapsed:.3e}/s",
            flush=True,
        )

    expected_first = 1.0 / (8.0 * (n // 2 + 2))
    if 1 in supports:
        first = aggregate_harmonic[p]
        if np.max(np.abs(first - expected_first)) > 2e-10:
            raise ArithmeticError(
                f"first-shell identity failed: {first} versus {expected_first}"
            )

    aggregate = {
        str(s): {
            "half_representatives": aggregate_counts[s],
            "phased_harmonic_sums": [
                float(value) for value in aggregate_harmonic[s]
            ],
        }
        for s in sorted(aggregate_counts)
    }
    result = {
        "experiment": "r1_sparse_dual_harmonic_gpu",
        "status": "finite_sparse_reconnaissance_not_proof",
        "p": p,
        "n": n,
        "supports": supports,
        "collision_free_union": 2 * max(supports) < p + 1,
        "seed": args.seed,
        "device": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        "compute_capability": cp.cuda.Device().compute_capability,
        "cupy": cp.__version__,
        "kernel_compile_seconds": compile_seconds,
        "admissible_matrices": diagnostics,
        "expected_first_shell_scalar": expected_first,
        "rows": rows,
        "aggregate_by_scaled_norm": aggregate,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
