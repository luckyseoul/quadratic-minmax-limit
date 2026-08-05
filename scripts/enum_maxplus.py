#!/usr/bin/env python3
"""Enumerate boolean +p eigenvectors (Max+) for Paley conference order n=p²+1."""
from __future__ import annotations

import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

# ensure importable
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def _worker(args: tuple) -> np.ndarray:
    start, count, nul, cpath = args
    nb = np.load(f"{cpath}/nb.npy")
    Binv = np.load(f"{cpath}/Binv.npy")
    found = []
    for bits in range(start, start + count):
        free = np.array(
            [1.0 if (bits >> i) & 1 else -1.0 for i in range(nul)], dtype=float
        )
        x = nb @ (Binv @ free)
        if np.max(np.abs(np.abs(x) - 1.0)) < 1e-6:
            found.append(np.sign(x))
    if not found:
        return np.zeros((0, nb.shape[0]))
    return np.unique(np.round(found).astype(int), axis=0).astype(np.float64)


def enum_maxplus(p: int, workers: int | None = None) -> Path:
    t0 = time.time()
    C = paley_conference_prime_power(p).astype(np.float64)
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
        if np.linalg.matrix_rank(nb[idx], tol=1e-8) == nul:
            coords = idx
            break
    if coords is None:
        raise RuntimeError(f"p={p}: no free coordinate set")
    Binv = np.linalg.inv(nb[coords])
    cache = Path(f"/tmp/e1_p{p}")
    cache.mkdir(exist_ok=True)
    np.save(cache / "nb.npy", nb)
    np.save(cache / "Binv.npy", Binv)
    np.save(cache / "coords.npy", coords)

    total = 1 << nul
    if workers is None:
        workers = min(86, max(1, (os.cpu_count() or 4) - 2))
    print(f"p={p} n={n} nul={nul} total={total} workers={workers}", flush=True)

    if nul <= 16:
        found = []
        for bits in range(total):
            free = np.array(
                [1.0 if (bits >> i) & 1 else -1.0 for i in range(nul)], dtype=float
            )
            x = nb @ (Binv @ free)
            if np.max(np.abs(np.abs(x) - 1.0)) < 1e-6:
                found.append(np.sign(x))
        Mp = np.unique(np.round(found).astype(int), axis=0).astype(np.float64)
    else:
        shard = (total + workers - 1) // workers
        tasks = []
        for w in range(workers):
            start = w * shard
            if start >= total:
                break
            tasks.append((start, min(shard, total - start), nul, str(cache)))
        chunks: list[np.ndarray] = []
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for i, arr in enumerate(ex.map(_worker, tasks, chunksize=1)):
                if len(arr):
                    chunks.append(arr)
                if (i + 1) % 10 == 0 or i + 1 == len(tasks):
                    print(f"  shard {i+1}/{len(tasks)}", flush=True)
        Mp = np.unique(np.vstack(chunks), axis=0)

    out = Path(f"/tmp/maxplus_p{p}.npy")
    np.save(out, Mp)
    if p == 7:
        np.save(cache / "maxplus.npy", Mp)
    print(f"|Max+|={len(Mp)} wrote {out} in {time.time()-t0:.1f}s", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("p", type=int, nargs="+")
    ap.add_argument("-j", type=int, default=None)
    args = ap.parse_args()
    for p in args.p:
        enum_maxplus(p, workers=args.j)


if __name__ == "__main__":
    main()
