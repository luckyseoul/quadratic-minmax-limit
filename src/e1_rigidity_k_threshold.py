#!/usr/bin/env python3
"""
Probe: for n=10 Paley C, among matrices with best_k(A,C) >= K, what is min Phi?
If min Phi >= Phi(C) for all K > k_star of true optima, supports
"undercutters must be Hamming-close".
"""
from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from minmax_quadratic import (  # noqa: E402
    edge_hamming,
    paley_conference_prime_power,
    phi,
    random_sym_pm1,
)

C = None
N = 10
PHI_C = 15.0


def _init():
    global C
    C = paley_conference_prime_power(3)


def best_k_exact(A: np.ndarray) -> int:
    best = edge_hamming(A, C)
    for bits in product([-1.0, 1.0], repeat=N - 1):
        eps = np.ones(N)
        eps[1:] = bits
        As = eps[:, None] * A * eps[None, :]
        best = min(best, edge_hamming(As, C))
    return best


def worker(args):
    K, seed, steps = args
    _init()
    rng = np.random.default_rng(seed)
    A = random_sym_pm1(N, rng)

    def energy(M):
        ph = phi(M)
        k = best_k_exact(M)
        pen = 0.0 if k >= K else 500.0 * (K - k) + (ph * 0.0)
        return ph + pen, ph, k

    cur_en, ph, k = energy(A)
    best_ph, best_k = (ph, k) if k >= K else (1e9, -1)
    for t in range(steps):
        i, j = rng.integers(0, N, size=2)
        if i == j:
            continue
        A[i, j] *= -1
        A[j, i] *= -1
        en, ph, k = energy(A)
        T = 4.0 * (0.08 / 4.0) ** (t / max(steps - 1, 1))
        if en <= cur_en or rng.random() < np.exp((cur_en - en) / max(T, 1e-9)):
            cur_en = en
            if k >= K and ph < best_ph:
                best_ph, best_k = ph, k
        else:
            A[i, j] *= -1
            A[j, i] *= -1
    return {"K": K, "seed": seed, "phi": best_ph, "k": best_k}


def main():
    W = max(1, (os.cpu_count() or 2) - 2)
    jobs = []
    for K in (5, 6, 7, 8, 10, 12, 15):
        for s in range(32):
            jobs.append((K, s + 17 * K, 6000))
    print(f"workers={W} jobs={len(jobs)}", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=W) as ex:
        for r in ex.map(worker, jobs, chunksize=1):
            rows.append(r)
            if r["k"] >= r["K"] and r["phi"] < 1e8:
                print(
                    f"K>={r['K']} seed={r['seed']} phi={r['phi']} k={r['k']}",
                    flush=True,
                )
    print("--- summary min Phi with best_k >= K ---", flush=True)
    for K in (5, 6, 7, 8, 10, 12, 15):
        phis = [r["phi"] for r in rows if r["K"] == K and r["k"] >= K and r["phi"] < 1e8]
        if phis:
            print(
                f"K>={K}: hits={len(phis)} min_phi={min(phis)} "
                f"(Phi_C={PHI_C}, undercut? {min(phis) < PHI_C})",
                flush=True,
            )
        else:
            print(f"K>={K}: no feasible hits", flush=True)


if __name__ == "__main__":
    main()
