#!/usr/bin/env python3
"""
n=26 SA minimising phi_local, then exact-Φ rescore of finalists (Gray/numba).

Produces a *certified upper bound* on m_26: min over candidates of exact Φ(A).
Does not claim m_26 equality.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import (
    is_conference_matrix,
    paley_conference_prime_power,
    phi_local,
    halfspace_boolean_vector,
    rho_from_boolean_evec,
)
from n26_matching_probe import exact_phi_numba, HAS_NUMBA

for _k in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_k, "1")


def _sa_worker(args):
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    n, iters, seed, start_paley, C_list = args
    rng = np.random.default_rng(seed)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    ne = len(edges)
    if start_paley and C_list is not None:
        A = np.array(C_list, dtype=np.float64)
        # scramble a few edges so we are not stuck at Paley
        for _ in range(int(rng.integers(3, 12))):
            i, j = edges[int(rng.integers(0, ne))]
            A[i, j] *= -1
            A[j, i] *= -1
    else:
        A = np.zeros((n, n), dtype=np.float64)
        for i, j in edges:
            s = float(rng.choice([-1.0, 1.0]))
            A[i, j] = A[j, i] = s
    cur = phi_local(A, restarts=24, rng=rng)
    best = cur
    best_A = A.copy()
    T0 = max(3.0, cur * 0.12)
    for t in range(iters):
        T = T0 * (1.0 - t / iters) + 1e-9
        e = int(rng.integers(0, ne))
        i, j = edges[e]
        A[i, j] *= -1
        A[j, i] *= -1
        # cheap local
        new = phi_local(A, restarts=6, rng=rng)
        d = new - cur
        if d <= 0 or rng.random() < math.exp(-d / T):
            cur = new
            if cur < best:
                best = cur
                best_A = A.copy()
        else:
            A[i, j] *= -1
            A[j, i] *= -1
    # heavier local restarts on best
    best = min(best, phi_local(best_A, restarts=80, rng=rng))
    return {"phi_local": float(best), "A": best_A.tolist(), "seed": seed}


def _rescore(args):
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    A = np.array(args["A"], dtype=np.float64)
    t0 = time.time()
    ph = float(exact_phi_numba(A))
    return {
        "seed": args["seed"],
        "phi_local": args["phi_local"],
        "phi_exact": ph,
        "time_s": time.time() - t0,
    }


def main():
    workers = int(os.environ.get("WORKERS", max(1, (os.cpu_count() or 4) - 2)))
    iters = int(os.environ.get("SA_ITERS", "4000"))
    n_sa = int(os.environ.get("N_SA", str(workers)))
    out_dir = Path(os.environ.get("OUT", str(ROOT / "evidence")))
    out_dir.mkdir(parents=True, exist_ok=True)

    p, n = 5, 26
    C = paley_conference_prime_power(p)
    assert is_conference_matrix(C)
    assert abs(rho_from_boolean_evec(C, halfspace_boolean_vector(p)) - 1) < 1e-8
    phi_paley = 0.5 * n * p  # 65

    print(f"SA wave: {n_sa} workers, iters={iters}, n={n}", flush=True)
    t0 = time.time()
    jobs = []
    C_list = C.tolist()
    for w in range(n_sa):
        start_paley = w % 3 == 0
        jobs.append((n, iters, 70000 + w, start_paley, C_list if start_paley else None))

    finals = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_sa_worker, j) for j in jobs]
        for fut in as_completed(futs):
            finals.append(fut.result())
            print(
                f"  SA done seed={finals[-1]['seed']} phi_local={finals[-1]['phi_local']:.1f}",
                flush=True,
            )
    print(f"SA wall {time.time()-t0:.1f}s; rescoring all with exact Φ ...", flush=True)

    # also score pure Paley
    finals.append({"phi_local": phi_paley, "A": C_list, "seed": -1})

    rescored = []
    t1 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_rescore, f) for f in finals]
        for fut in as_completed(futs):
            r = fut.result()
            rescored.append(r)
            print(
                f"  exact Φ={r['phi_exact']:.1f} (local was {r['phi_local']:.1f}) "
                f"seed={r['seed']} t={r['time_s']:.2f}s",
                flush=True,
            )

    exacts = [r["phi_exact"] for r in rescored]
    report = {
        "n": n,
        "Phi_Paley": phi_paley,
        "has_numba": HAS_NUMBA,
        "n_sa": n_sa,
        "sa_iters": iters,
        "workers": workers,
        "certified_UB_m26": min(exacts),
        "paley_exact_phi": next(r["phi_exact"] for r in rescored if r["seed"] == -1),
        "n_undercut_paley": sum(1 for e in exacts if e < phi_paley - 1e-9),
        "sa_time_s": t1 - t0,
        "rescore_time_s": time.time() - t1,
        "rescored": sorted(rescored, key=lambda z: z["phi_exact"]),
        "note": (
            "Certified upper bound m_26 ≤ min exact Φ among SA finalists + Paley. "
            "phi_local underestimates Φ; only exact column is valid for UB."
        ),
    }
    out = out_dir / "n26_sa_exact_rescore.json"
    out.write_text(json.dumps(report, indent=2))
    print(f"wrote {out}", flush=True)
    print(
        f"CERTIFIED m_26 ≤ {report['certified_UB_m26']}  "
        f"(Paley Φ={phi_paley}, undercuts={report['n_undercut_paley']})",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
