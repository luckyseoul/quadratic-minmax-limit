#!/usr/bin/env python3
"""
E(1) exact-Phi SA at n=26 using phi_mitm rescoring.

ProcessPool requires a script file (stdin workers fail).
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

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "1"

from minmax_quadratic import paley_conference_prime_power, phi_local, phi_mitm


def worker(seed: int) -> dict:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    C = paley_conference_prime_power(5)
    n = 26
    rng = np.random.default_rng(seed)
    A = C.copy()
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    k0 = int(rng.integers(0, 50))
    for t in rng.choice(len(edges), size=min(k0, len(edges)), replace=False):
        i, j = edges[t]
        A[i, j] *= -1
        A[j, i] *= -1
    cur = phi_local(A, restarts=6, rng=rng)
    best_local = cur
    bestA = A.copy()
    steps = 2500
    for t in range(steps):
        i, j = rng.integers(0, n, size=2)
        if i == j:
            continue
        A[i, j] *= -1
        A[j, i] *= -1
        sc = phi_local(A, restarts=2, rng=rng)
        T = 5.0 * (0.08 / 5) ** (t / steps)
        if sc <= cur or rng.random() < np.exp((cur - sc) / max(T, 1e-9)):
            cur = sc
            if sc < best_local:
                best_local = sc
                bestA = A.copy()
        else:
            A[i, j] *= -1
            A[j, i] *= -1
    exact = phi_mitm(bestA)
    return {"seed": seed, "phi_local": best_local, "phi_exact": exact}


def main() -> None:
    nseeds = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    nworkers = max(1, (os.cpu_count() or 4) - 2)
    Phi = 65.0
    print(f"n26 MITM-SA seeds={nseeds} workers={nworkers}", flush=True)
    results = []
    with ProcessPoolExecutor(max_workers=nworkers) as ex:
        futs = [ex.submit(worker, s) for s in range(nseeds)]
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            if r["phi_exact"] < Phi - 0.5:
                print("UNDERCUT", r, flush=True)
            if len(results) % 10 == 0:
                mx = min(x["phi_exact"] for x in results)
                print(f"  done {len(results)} min_exact={mx}", flush=True)
    exacts = [r["phi_exact"] for r in results]
    out = {
        "n": 26,
        "Phi_C": Phi,
        "nseeds": nseeds,
        "min_exact": min(exacts),
        "max_exact": max(exacts),
        "n_undercut": sum(1 for e in exacts if e < Phi - 0.5),
        "results": results,
    }
    path = ROOT / "evidence" / "e1_n26_mitm_sa.json"
    path.write_text(json.dumps(out, indent=2))
    print("SUMMARY", {k: out[k] for k in out if k != "results"})
    print("wrote", path)


if __name__ == "__main__":
    main()
