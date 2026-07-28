#!/usr/bin/env python3
"""
Parallel undercut hunt at Paley n=26: random edge sets by k and degree cap.
Uses phi_mitm. ProcessPool with many workers.
Writes evidence/e1_n26_parallel_undercut_hunt.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power, phi_mitm  # noqa: E402


def worker(job):
    seed, n_samples, k, delta_cap = job
    C = paley_conference_prime_power(5)
    n = C.shape[0]
    Phi = 0.5 * n * 5
    rng = np.random.default_rng(seed)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    best = 10**9
    under = []
    for _t in range(n_samples):
        deg = np.zeros(n, dtype=int)
        chosen = []
        order = rng.permutation(E)
        for ei in order:
            if len(chosen) >= k:
                break
            i, j = edges[int(ei)]
            if deg[i] >= delta_cap or deg[j] >= delta_cap:
                continue
            chosen.append((i, j))
            deg[i] += 1
            deg[j] += 1
        if len(chosen) < k:
            have = {frozenset(e) for e in chosen}
            for ei in order:
                if len(chosen) >= k:
                    break
                i, j = edges[int(ei)]
                fs = frozenset((i, j))
                if fs in have:
                    continue
                chosen.append((i, j))
                have.add(fs)
                deg[i] += 1
                deg[j] += 1
        A = C.copy()
        for i, j in chosen:
            A[i, j] *= -1
            A[j, i] *= -1
        ph = float(phi_mitm(A))
        if ph < best:
            best = ph
        if ph < Phi - 0.5:
            under.append(
                {
                    "k": k,
                    "delta_cap": delta_cap,
                    "phi": ph,
                    "Delta": int(deg.max()),
                }
            )
    return {
        "seed": seed,
        "k": k,
        "delta_cap": delta_cap,
        "n": n_samples,
        "best_phi": best,
        "undercuts": under,
    }


def main():
    W = min(80, (os.cpu_count() or 4) - 2)
    jobs = []
    seed = 0
    for k in [5, 10, 13, 15, 20, 26, 40, 60]:
        for dcap in [1, 2, 3, 5, 10, 25]:
            for _rep in range(max(1, W // 48)):
                jobs.append((seed, 30, k, dcap))
                seed += 1
    while len(jobs) < W:
        jobs.append((seed, 30, 13, 2))
        seed += 1
    print(f"jobs={len(jobs)} workers={W}", flush=True)
    results = []
    under_all = []
    best_global = 10**9
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = [ex.submit(worker, j) for j in jobs]
        for i, f in enumerate(as_completed(futs)):
            r = f.result()
            results.append(r)
            best_global = min(best_global, r["best_phi"])
            under_all.extend(r["undercuts"])
            if (i + 1) % 20 == 0:
                print(
                    f"  done {i+1}/{len(jobs)} best_so_far={best_global}",
                    flush=True,
                )
    out = {
        "n_jobs": len(jobs),
        "workers": W,
        "total_samples_est": int(sum(r["n"] for r in results)),
        "best_phi_seen": float(best_global),
        "Phi_C": 65.0,
        "n_undercuts": len(under_all),
        "undercuts": under_all[:20],
        "status": (
            "OPEN — parallel random F hunt; 0 undercuts is evidence not a proof. "
            "Existence of lim alpha_n remains OPEN."
        ),
    }
    path = ROOT / "evidence" / "e1_n26_parallel_undercut_hunt.json"
    path.write_text(json.dumps(out, indent=2))
    print(
        "wrote",
        path,
        "best",
        best_global,
        "undercuts",
        len(under_all),
        flush=True,
    )


if __name__ == "__main__":
    main()
