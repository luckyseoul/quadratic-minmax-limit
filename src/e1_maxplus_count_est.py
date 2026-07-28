#!/usr/bin/env python3
"""Estimate |Max+| at given p via free-coordinate hit rate (ProcessPool)."""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CACHE = Path("/tmp/e1_count_est")
CACHE.mkdir(exist_ok=True)


def prep(p: int):
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    M = C - p * np.eye(n)
    _, s, vt = np.linalg.svd(M)
    rank = int((s > 1e-8).sum())
    nul = n - rank
    nb = vt[rank:].T
    rng = np.random.default_rng(0)
    coords = None
    for _ in range(50000):
        idx = rng.choice(n, size=nul, replace=False)
        if np.linalg.matrix_rank(nb[idx], tol=1e-6) == nul:
            coords = idx
            break
    if coords is None:
        raise RuntimeError("no free coords")
    A = nb @ np.linalg.inv(nb[coords])
    np.save(CACHE / f"A_p{p}.npy", A)
    return n, nul


def worker(args):
    p, seed, trials = args
    A = np.load(CACHE / f"A_p{p}.npy")
    nul = A.shape[1]
    rng = np.random.default_rng(seed)
    found = 0
    B = 5000
    left = trials
    while left > 0:
        b = min(B, left)
        F = rng.choice([-1.0, 1.0], size=(b, nul))
        X = F @ A.T
        found += int(np.sum(np.max(np.abs(np.abs(X) - 1.0), axis=1) < 1e-6))
        left -= b
    return found


def estimate(p: int, trials_per_worker: int = 100000, workers: int = 80):
    n, nul = prep(p)
    print(f"p={p} n={n} nul={nul} 2^nul={2**nul:.4e}", flush=True)
    t0 = time.time()
    tasks = [(p, s, trials_per_worker) for s in range(workers)]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        hits = list(ex.map(worker, tasks, chunksize=1))
    total_hits = sum(hits)
    total_trials = workers * trials_per_worker
    rate = total_hits / total_trials
    N_est = rate * (2**nul)
    out = {
        "p": p,
        "n": n,
        "nul": nul,
        "hits": total_hits,
        "trials": total_trials,
        "rate": rate,
        "N_est": N_est,
        "seconds": time.time() - t0,
    }
    print(json.dumps(out, indent=2), flush=True)
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=int, nargs="+", default=[7, 11])
    ap.add_argument("--trials", type=int, default=100000)
    ap.add_argument("--workers", type=int, default=0)
    args = ap.parse_args()
    w = args.workers or min(80, max(1, (os.cpu_count() or 4) - 2))
    results = []
    for p in args.p:
        results.append(estimate(p, trials_per_worker=args.trials, workers=w))
    path = ROOT / "evidence" / "e1_maxplus_count_est.json"
    path.write_text(json.dumps({"results": results}, indent=2))
    print("wrote", path)


if __name__ == "__main__":
    main()
