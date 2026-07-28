#!/usr/bin/env python3
"""
Full boolean +p evec enumeration at p=7 (2^25 free coords, ProcessPool),
then exact gmin over all disjoint edge pairs (or a complete type census).

Writes evidence/e1_gmin_p7.json
"""
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

CACHE = Path("/tmp/e1_p7")
CACHE.mkdir(exist_ok=True)


def _prep():
    p = 7
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    M = C - p * np.eye(n)
    _, s, vt = np.linalg.svd(M)
    rank = int((s > 1e-8).sum())
    nul = n - rank
    nb = vt[rank:].T
    rng = np.random.default_rng(0)
    coords = None
    for _ in range(20000):
        idx = rng.choice(n, size=nul, replace=False)
        if np.linalg.matrix_rank(nb[idx], tol=1e-8) == nul:
            coords = idx
            break
    if coords is None:
        raise RuntimeError("no free coords")
    Binv = np.linalg.inv(nb[coords])
    np.save(CACHE / "nb.npy", nb)
    np.save(CACHE / "Binv.npy", Binv)
    np.save(CACHE / "coords.npy", coords)
    np.save(CACHE / "C.npy", C)
    return p, n, nul


def _worker(args):
    start, count, nul = args
    nb = np.load(CACHE / "nb.npy")
    Binv = np.load(CACHE / "Binv.npy")
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


def main() -> None:
    t0 = time.time()
    p, n, nul = _prep()
    total = 1 << nul
    workers = min(80, max(1, (os.cpu_count() or 4) - 2))
    # shard bit ranges
    shard = (total + workers - 1) // workers
    tasks = []
    for w in range(workers):
        start = w * shard
        if start >= total:
            break
        count = min(shard, total - start)
        tasks.append((start, count, nul))

    print(f"p={p} n={n} nul={nul} total={total} workers={workers} shards={len(tasks)}", flush=True)
    chunks = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for arr in ex.map(_worker, tasks, chunksize=1):
            if len(arr):
                chunks.append(arr)
            print(f"  chunk size {len(arr)} cumulative_chunks {len(chunks)}", flush=True)

    if not chunks:
        raise RuntimeError("no evecs")
    Mp = np.unique(np.vstack(chunks), axis=0)
    # include both signs? boolean evecs come in ± pairs from free signs already
    print(f"Max+ count {len(Mp)} (expected 11452)", flush=True)
    np.save(CACHE / "maxplus.npy", Mp)
    np.save(ROOT / "evidence" / "maxplus_p7.npy", Mp)  # large - maybe skip
    # actually evidence npy may be big; keep in /tmp only
    Path(ROOT / "evidence" / "maxplus_p7.npy").unlink(missing_ok=True)

    C = np.load(CACHE / "C.npy")
    N = len(Mp)
    # Exact gmin: use all 4-sets is binom(50,4)=230300 - fine
    # For each 4-set, 3 matchings
    from itertools import combinations

    gmin = 1.0
    gmax = -1.0
    # also track histogram of rounded G
    from collections import Counter

    hist = Counter()
    n_pairs = 0
    for quad in combinations(range(n), 4):
        i, j, k, l = quad
        # three perfect matchings
        for pairs in (((i, j), (k, l)), ((i, k), (j, l)), ((i, l), (j, k))):
            (a, b), (c, d) = pairs
            m4 = float(np.mean(Mp[:, a] * Mp[:, b] * Mp[:, c] * Mp[:, d]))
            g = float(C[a, b] * C[c, d] * m4)
            gmin = min(gmin, g)
            gmax = max(gmax, g)
            hist[round(g, 8)] += 1
            n_pairs += 1

    out = {
        "p": p,
        "n": n,
        "n_Max_plus": int(N),
        "expected_Max_plus": 11452,
        "gmin_disjoint": gmin,
        "gmax_disjoint": gmax,
        "n_disjoint_oriented_pairs": n_pairs,
        "minus_1_over_p": -1.0 / p,
        "minus_1_over_15": -1.0 / 15,
        "star_force_level1": gmin > -1.0 / p + 1e-12,
        "matching_blocked_level2": gmin > -1.0 / 15 + 1e-12,
        "g_histogram_top": hist.most_common(15),
        "workers": workers,
        "seconds": time.time() - t0,
    }
    path = ROOT / "evidence" / "e1_gmin_p7.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print("wrote", path)


if __name__ == "__main__":
    main()
