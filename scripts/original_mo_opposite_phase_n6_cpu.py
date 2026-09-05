#!/usr/bin/env python3
"""Fixed n=6 opposite-phase Gibbs check, independently implemented for NumPy.

This is floating-point numerical evidence, not an optimizer certificate or a
convergence proof.  The discrete signing/state enumeration is exhaustive.
No larger orders, temperature search, or stochastic optimizer is implemented.
"""

import argparse
import datetime as dt
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import time

import numpy as np


N = 6
CS = (0.5, 1.0, 2.0, 4.0, 8.0)
EDGES = tuple(itertools.combinations(range(N), 2))
FREE_EDGES = tuple(itertools.combinations(range(1, N), 2))
TIE_ATOL = 2e-11
CHECK_ATOL = 2e-10
THREAD_KEYS = (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS",
)


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def signing(mask):
    a = np.ones((N, N), dtype=np.int64)
    np.fill_diagonal(a, 0)
    for bit, (i, j) in enumerate(FREE_EDGES):
        a[i, j] = a[j, i] = 1 - 2 * ((mask >> bit) & 1)
    return a


def stable_phase(energy, beta):
    exponent = beta * energy
    shift = float(np.max(exponent))
    weight = np.exp(exponent - shift)
    total = math.fsum(float(w) for w in weight)
    return shift + math.log(total / len(energy)), weight / total


def canonical_mask(a, permutations):
    # All vertex permutations, then the unique gauge with row zero positive.
    permuted = a[permutations[:, :, None], permutations[:, None, :]]
    gauge = permuted[:, 0, :].copy()
    gauge[:, 0] = 1
    normalized = permuted * gauge[:, :, None] * gauge[:, None, :]
    codes = np.zeros(len(permutations), dtype=np.int64)
    for bit, (i, j) in enumerate(FREE_EDGES):
        codes += (normalized[:, i, j] < 0).astype(np.int64) << bit
    # Global sign reversal, followed by gauge normalization, complements bits.
    return int(min(np.min(codes), 1023 - np.max(codes)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(all(os.environ.get(k) == "1" for k in THREAD_KEYS),
            "all six BLAS/OpenMP thread variables must equal 1")
    started = utc_now()
    clock = time.perf_counter()
    print(json.dumps({"event": "start", "utc": started, "pid": os.getpid(),
                      "ppid": os.getppid(), "pgid": os.getpgrp(),
                      "sid": os.getsid(0), "host": platform.node()}), flush=True)

    spins = np.array(list(itertools.product((-1, 1), repeat=N)), dtype=np.int64)
    pair_products = np.array([spins[:, i] * spins[:, j] for i, j in EDGES])
    permutations = np.array(list(itertools.permutations(range(N))), dtype=np.int64)
    records = [[] for _ in CS]
    energy_histograms = {}
    max_trace_error = 0.0
    max_derivative_error = 0.0
    checks = 0

    for mask in range(1024):
        a = signing(mask)
        edge_signs = np.array([a[i, j] for i, j in EDGES], dtype=np.int64)
        energy = edge_signs @ pair_products
        # These identities are exact integer checks, not floating comparisons.
        require(np.array_equal(energy, np.sum((spins @ a) * spins, axis=1) // 2),
                f"quadratic normalization failed at {mask}")
        require(int(np.sum(energy)) == 0 and int(energy @ energy) == 64 * 15,
                f"uniform first/second moments failed at {mask}")
        checks += 2
        values, counts = np.unique(energy, return_counts=True)
        energy_histograms[mask] = [[int(v), int(k)] for v, k in zip(values, counts)]
        cross_square = (spins @ a @ spins.T) ** 2

        for ci, c in enumerate(CS):
            beta = c / math.sqrt(N)
            log_plus, plus = stable_phase(energy, beta)
            log_minus, minus = stable_phase(energy, -beta)
            u = spins.T @ (plus[:, None] * spins)
            v = spins.T @ (minus[:, None] * spins)
            half_product = (log_plus + log_minus) / 2.0
            derivative = (float(plus @ energy) - float(minus @ energy)) / 2.0
            trace = float(np.trace(a @ u @ a @ v))
            trace_direct = float(plus @ cross_square @ minus)
            derivative_trace = float(np.trace(a @ (u - v))) / 4.0
            trace_error = abs(trace - trace_direct)
            derivative_error = abs(derivative - derivative_trace)
            max_trace_error = max(max_trace_error, trace_error)
            max_derivative_error = max(max_derivative_error, derivative_error)
            require(trace_error < CHECK_ATOL,
                    f"mixed-replica trace identity failed {mask},{c}: {trace_error}")
            require(derivative_error < CHECK_ATOL,
                    f"derivative identity failed {mask},{c}: {derivative_error}")
            require(np.max(abs(np.diag(u) - 1)) < CHECK_ATOL
                    and np.max(abs(np.diag(v) - 1)) < CHECK_ATOL,
                    f"covariance normalization failed {mask},{c}")
            require(trace > -CHECK_ATOL, f"negative squared trace {mask},{c}")
            checks += 4
            records[ci].append({
                "mask": mask, "a": half_product, "a_prime": derivative,
                "T": trace, "T_over_n_squared": trace / N**2,
                "gap_over_n": (beta * trace / 2.0 - derivative) / N,
            })

    canonical_cache = {}
    summaries = []
    for c, rows in zip(CS, records):
        minimum = min(row["a"] for row in rows)
        tied = [r for r in rows if abs(r["a"] - minimum) <= TIE_ATOL]
        other = [r["a"] for r in rows if r["a"] > minimum + TIE_ATOL]
        types = {}
        for row in tied:
            mask = row["mask"]
            if mask not in canonical_cache:
                canonical_cache[mask] = canonical_mask(signing(mask), permutations)
            key = canonical_cache[mask]
            types.setdefault(key, []).append(row)
        type_records = []
        for key, members in sorted(types.items()):
            representative = min(members, key=lambda r: r["mask"])
            type_records.append({
                "canonical_switch_permutation_sign_mask": key,
                "multiplicity_in_normalized_enumeration": len(members),
                "masks": [r["mask"] for r in members],
                "representative": representative,
                "representative_matrix": signing(representative["mask"]).tolist(),
                "exact_energy_histogram": energy_histograms[representative["mask"]],
                "T_over_n_squared_range": [min(r["T_over_n_squared"] for r in members),
                                           max(r["T_over_n_squared"] for r in members)],
                "gap_over_n_range": [min(r["gap_over_n"] for r in members),
                                      max(r["gap_over_n"] for r in members)],
            })
        summaries.append({
            "c": c, "beta": c / math.sqrt(N), "minimum_a": minimum,
            "next_non_tied_a_gap": min(other) - minimum if other else None,
            "tied_count": len(tied), "type_count": len(type_records),
            "candidate_minimizing_types": type_records,
        })

    finished = utc_now()
    result = {
        "classification": "exhaustive discrete enumeration with float64 evaluation; numerical candidate minima only",
        "not_a_convergence_proof": True,
        "started_utc": started, "finished_utc": finished,
        "elapsed_seconds": time.perf_counter() - clock,
        "host": platform.node(), "pid": os.getpid(), "ppid": os.getppid(),
        "pgid": os.getpgrp(), "sid": os.getsid(0),
        "python": platform.python_version(), "numpy": np.__version__,
        "thread_environment": {k: os.environ[k] for k in THREAD_KEYS},
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "n": N, "cs": CS, "signing_count": 1024, "spin_count": 64,
        "Q_definition": "sum_(i<j) A_ij x_i x_j = x^T A x / 2",
        "mask_definition": "A_0j=+1; bit k=1 means -1 on kth lexicographic edge of vertices 1,...,5",
        "free_edges": FREE_EDGES,
        "type_equivalence": "diagonal switching, all vertex permutations, and global sign reversal",
        "tie_absolute_tolerance": TIE_ATOL,
        "checks": checks, "max_trace_double_replica_error": max_trace_error,
        "max_derivative_trace_error": max_derivative_error,
        "summaries": summaries,
        "all_values_by_c": [{"c": c, "rows": rows} for c, rows in zip(CS, records)],
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"event": "complete", "utc": finished,
                      "elapsed_seconds": result["elapsed_seconds"], "checks": checks,
                      "output": str(args.output), "summaries": summaries}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
