#!/usr/bin/env python3
"""
SA maximizing min_{y in Max+} S_M(y) over perfect matchings of Paley n=p^2+1.

For odd |M|=n/2 (always true for odd p), S_M(y) is always odd, so either
min S >= 1 (Max-cover, candidate undercutter) or min S <= -1 (Phi jumps to
Phi+2 on that Max+ vector alone).

Usage:
  python3 src/e1_matching_mins_sa.py
Writes evidence/e1_matching_minS_sa.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def boolean_plus(C: np.ndarray, p: float) -> np.ndarray:
    n = C.shape[0]
    M = C - p * np.eye(n)
    _, s, vt = np.linalg.svd(M)
    rank = int((s > 1e-8).sum())
    nb = vt[rank:].T
    nul = n - rank
    rng = np.random.default_rng(0)
    coords = None
    for _ in range(8000):
        idx = rng.choice(n, size=nul, replace=False)
        if np.linalg.matrix_rank(nb[idx], tol=1e-8) == nul:
            coords = idx
            break
    if coords is None:
        raise RuntimeError("no free coords for boolean evecs")
    Binv = np.linalg.inv(nb[coords])
    found = []
    for bits in range(1 << nul):
        free = np.array([1.0 if (bits >> i) & 1 else -1.0 for i in range(nul)])
        x = nb @ (Binv @ free)
        if np.max(np.abs(np.abs(x) - 1.0)) < 1e-6:
            found.append(np.sign(x))
    return np.unique(np.round(found).astype(int), axis=0).astype(float)


def pm_min_S_worker(args: tuple) -> dict:
    p, seed, iters = args
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    Max = boolean_plus(C, float(p))
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    M = [(int(min(perm[2 * i], perm[2 * i + 1])), int(max(perm[2 * i], perm[2 * i + 1]))) for i in range(n // 2)]

    def min_S(matching: list) -> tuple[float, float, float]:
        s = np.zeros(len(Max))
        for i, j in matching:
            s += C[i, j] * Max[:, i] * Max[:, j]
        return float(s.min()), float(s.mean()), float(s.max())

    best = min_S(M)
    cur = best
    for t in range(iters):
        a, b = rng.choice(n // 2, size=2, replace=False)
        (u, v), (w, x) = M[a], M[b]
        for pairs in (
            [(min(u, w), max(u, w)), (min(v, x), max(v, x))],
            [(min(u, x), max(u, x)), (min(v, w), max(v, w))],
        ):
            if pairs[0][0] == pairs[0][1] or pairs[1][0] == pairs[1][1]:
                continue
            if len({pairs[0], pairs[1]}) < 2:
                continue
            verts = [pairs[0][0], pairs[0][1], pairs[1][0], pairs[1][1]]
            if len(set(verts)) < 4:
                continue
            newM = M[:]
            newM[a], newM[b] = pairs[0], pairs[1]
            sc = min_S(newM)
            temp = 3.0 + 10.0 * t / max(iters, 1)
            if sc[0] > cur[0] or (abs(sc[0] - cur[0]) < 1e-12 and sc[1] > cur[1]) or rng.random() < np.exp(
                (sc[0] - cur[0]) * temp
            ):
                M = newM
                cur = sc
                if cur[0] > best[0]:
                    best = cur
    return {
        "seed": int(seed),
        "best_min_S": best[0],
        "best_mean_S": best[1],
        "best_max_S": best[2],
    }


def run_p(p: int, n_seeds: int, iters: int, workers: int) -> dict:
    C = paley_conference_prime_power(p)
    Max = boolean_plus(C, float(p))
    n = C.shape[0]
    print(
        f"p={p} n={n} |Max+|={len(Max)} n/2={n // 2} odd={bool(n // 2 % 2)}",
        flush=True,
    )
    t0 = time.time()
    args = [(p, 1000 + s, iters) for s in range(n_seeds)]
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(pm_min_S_worker, a) for a in args]
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            if r["best_min_S"] >= 1 - 1e-9:
                print("COVER FOUND", r, flush=True)
    mins = [r["best_min_S"] for r in rows]
    summary = {
        "p": p,
        "n": int(n),
        "n_Max": int(len(Max)),
        "matching_size": int(n // 2),
        "matching_size_odd": bool((n // 2) % 2),
        "n_seeds": n_seeds,
        "iters": iters,
        "global_best_min_S": float(max(mins)),
        "global_worst_seed_min_S": float(min(mins)),
        "mean_of_seed_bests": float(np.mean(mins)),
        "all_min_S_le_neg1": all(m <= -1 + 1e-9 for m in mins),
        "any_cover": any(m >= 1 - 1e-9 for m in mins),
        "best_minS_counter": {str(k): int(v) for k, v in sorted(Counter(mins).items())},
        "seconds": time.time() - t0,
    }
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main() -> None:
    workers = max(1, (os.cpu_count() or 2) - 2)
    s3 = run_p(3, n_seeds=24, iters=10000, workers=workers)
    s5 = run_p(5, n_seeds=48, iters=25000, workers=workers)
    out = {"p3": s3, "p5": s5, "note": "OPEN: SA evidence, not a general proof of no matching Max-cover"}
    path = ROOT / "evidence" / "e1_matching_minS_sa.json"
    path.write_text(json.dumps(out, indent=2))
    print("wrote", path, flush=True)
    # Sanity: p=3 must find a cover (144 undercutting matchings exist)
    assert s3["any_cover"], "p=3 SA failed to recover known matching covers"
    # p=5: covers can exist (see E1_MATCHING_COVER_SPIKE.md); do not assert absence


if __name__ == "__main__":
    main()
