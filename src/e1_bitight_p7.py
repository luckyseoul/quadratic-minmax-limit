#!/usr/bin/env python3
"""
Bi-tight integral feasibility at p=7 level 2.
Max+ from /tmp/e1_p7/maxplus.npy; Max- enumerated via -C (module-level workers).

Writes evidence/e1_bitight_p7.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CACHE = Path("/tmp/e1_p7")


def _mm_worker(args):
    start, count, nul = args
    nb = np.load(CACHE / "nb_m.npy")
    Binv = np.load(CACHE / "Binv_m.npy")
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


def enum_maxminus(C: np.ndarray, p: float) -> np.ndarray:
    path = CACHE / "maxminus.npy"
    if path.is_file() and np.load(path).shape[0] > 1000:
        return np.load(path).astype(float)
    M = -C - p * np.eye(C.shape[0])
    _, s, vt = np.linalg.svd(M)
    rank = int((s > 1e-8).sum())
    nul = C.shape[0] - rank
    nb = vt[rank:].T
    rng = np.random.default_rng(1)
    coords = None
    for _ in range(20000):
        idx = rng.choice(C.shape[0], size=nul, replace=False)
        if np.linalg.matrix_rank(nb[idx], tol=1e-8) == nul:
            coords = idx
            break
    Binv = np.linalg.inv(nb[coords])
    np.save(CACHE / "nb_m.npy", nb)
    np.save(CACHE / "Binv_m.npy", Binv)
    total = 1 << nul
    workers = min(80, max(1, (os.cpu_count() or 4) - 2))
    shard = (total + workers - 1) // workers
    tasks = []
    for w in range(workers):
        start = w * shard
        if start >= total:
            break
        tasks.append((start, min(shard, total - start), nul))
    chunks = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for arr in ex.map(_mm_worker, tasks, chunksize=1):
            if len(arr):
                chunks.append(arr)
    Mm = np.unique(np.vstack(chunks), axis=0)
    np.save(path, Mm)
    return Mm


def main() -> None:
    t0 = time.time()
    p = 7
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    Mp = np.load(CACHE / "maxplus.npy").astype(float)
    print("Max+", len(Mp), flush=True)
    Mm = enum_maxminus(C, float(p))
    print("Max-", len(Mm), flush=True)

    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    print("building Chip/Chim...", flush=True)
    Chip = np.array(
        [[C[a, b] * y[a] * y[b] for a, b in edges] for y in Mp], dtype=np.float64
    )
    Chim = np.array(
        [[C[a, b] * z[a] * z[b] for a, b in edges] for z in Mm], dtype=np.float64
    )

    s_lev = 2
    k = s_lev * p
    Aeq = np.vstack([Chip, Chim, np.ones((1, E))])
    beq = np.concatenate(
        [s_lev * np.ones(len(Mp)), -s_lev * np.ones(len(Mm)), [float(k)]]
    )
    print(f"frac s={s_lev} k={k}...", flush=True)
    resf = linprog(np.zeros(E), A_eq=Aeq, b_eq=beq, bounds=(0, 1), method="highs")
    print("frac", resf.success, flush=True)
    print("integral time_limit=900...", flush=True)
    resi = milp(
        c=np.zeros(E),
        constraints=LinearConstraint(Aeq, beq, beq),
        integrality=np.ones(E),
        bounds=Bounds(0, 1),
        options={"time_limit": 900},
    )
    msg = str(resi.message)
    out = {
        "p": p,
        "n": n,
        "n_Max_plus": int(len(Mp)),
        "n_Max_minus": int(len(Mm)),
        "level_s": s_lev,
        "k": k,
        "fractional_feasible": bool(resf.success),
        "integral_feasible": bool(resi.success),
        "integral_infeasible": (not resi.success) and ("infeas" in msg.lower()),
        "integral_timeout": "time limit" in msg.lower(),
        "message": msg[:200],
        "seconds": time.time() - t0,
    }
    path = ROOT / "evidence" / "e1_bitight_p7.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print("wrote", path)


if __name__ == "__main__":
    main()
