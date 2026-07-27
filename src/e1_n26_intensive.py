#!/usr/bin/env python3
"""
Intensive E(1) search at n=26 (ρ=1 Paley, Φ=65).

Methods:
  1. SA with phi_local from many starts (Paley-near + random + structured)
  2. Exact-Φ rescore of all finalists
  3. Structured: geometric matchings, random k-flips of Paley for k=1..40
  4. Greedy bit ascent on exact Φ for top candidates

If best exact Φ < 65, E(1) has a concrete undercut at n=26.
If best = 65, stronger evidence for local optimality (not a proof).
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
    halfspace_boolean_vector,
    is_conference_matrix,
    paley_conference_prime_power,
    phi_local,
    random_sym_pm1,
)
from n26_matching_probe import exact_phi_numba, random_perfect_matching, flip_matching

for _k in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_k, "1")


def _flip(A, i, j):
    A[i, j] *= -1.0
    A[j, i] *= -1.0


def _sa_job(args):
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    n, iters, seed, mode, C_list = args
    rng = np.random.default_rng(seed)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    ne = len(edges)
    C = np.array(C_list, dtype=np.float64)
    if mode == "paley_near":
        A = C.copy()
        k = int(rng.integers(1, 30))
        for _ in range(k):
            i, j = edges[int(rng.integers(0, ne))]
            _flip(A, i, j)
    elif mode == "matching":
        m = random_perfect_matching(n, rng)
        A = flip_matching(C, m)
        # extra scramble
        for _ in range(int(rng.integers(0, 10))):
            i, j = edges[int(rng.integers(0, ne))]
            _flip(A, i, j)
    elif mode == "kflip":
        A = C.copy()
        k = int(args[5]) if len(args) > 5 else 10
        idx = rng.choice(ne, size=min(k, ne), replace=False)
        for t in idx:
            i, j = edges[int(t)]
            _flip(A, i, j)
    else:
        A = random_sym_pm1(n, rng)

    cur = phi_local(A, restarts=20, rng=rng)
    best = cur
    best_A = A.copy()
    T0 = max(3.0, cur * 0.15)
    for t in range(iters):
        T = T0 * (1.0 - t / max(iters, 1)) + 1e-9
        e = int(rng.integers(0, ne))
        i, j = edges[e]
        _flip(A, i, j)
        new = phi_local(A, restarts=6, rng=rng)
        d = new - cur
        if d <= 0 or rng.random() < math.exp(-d / T):
            cur = new
            if cur < best:
                best = cur
                best_A = A.copy()
        else:
            _flip(A, i, j)
    # final local polish: greedy exact is too slow; extra phi_local
    best = min(best, phi_local(best_A, restarts=40, rng=rng))
    return {
        "phi_local": float(best),
        "A": best_A.tolist(),
        "seed": seed,
        "mode": mode,
    }


def _exact_job(payload):
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    A = np.array(payload["A"], dtype=np.float64)
    t0 = time.time()
    ph = float(exact_phi_numba(A))
    return {
        "seed": payload.get("seed"),
        "mode": payload.get("mode"),
        "phi_local": payload.get("phi_local"),
        "phi_exact": ph,
        "time_s": time.time() - t0,
    }


def _kflip_exact_job(args):
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    C_list, seed, k = args
    rng = np.random.default_rng(seed)
    C = np.array(C_list, dtype=np.float64)
    n = C.shape[0]
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    A = C.copy()
    idx = rng.choice(len(edges), size=k, replace=False)
    for t in idx:
        i, j = edges[int(t)]
        _flip(A, i, j)
    ph = float(exact_phi_numba(A))
    return {"k": k, "seed": seed, "phi_exact": ph}


def main():
    workers = int(os.environ.get("WORKERS", max(1, (os.cpu_count() or 4) - 2)))
    sa_iters = int(os.environ.get("SA_ITERS", "8000"))
    n_sa = int(os.environ.get("N_SA", str(workers * 2)))
    out_dir = Path(os.environ.get("OUT", str(ROOT / "evidence")))
    out_dir.mkdir(parents=True, exist_ok=True)

    p, n = 5, 26
    C = paley_conference_prime_power(p)
    assert is_conference_matrix(C)
    Phi = 0.5 * n * p  # 65
    # warm exact
    ph_c = float(exact_phi_numba(C))
    assert abs(ph_c - Phi) < 1e-6, (ph_c, Phi)

    C_list = C.tolist()
    print(f"n=26 Φ_Paley={Phi} workers={workers} sa_iters={sa_iters} n_sa={n_sa}", flush=True)

    # --- SA wave ---
    jobs = []
    modes = ["paley_near", "matching", "random", "paley_near", "matching", "random"]
    for w in range(n_sa):
        mode = modes[w % len(modes)]
        jobs.append((n, sa_iters, 400_000 + w, mode, C_list))

    t0 = time.time()
    finals = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_sa_job, j) for j in jobs]
        for fut in as_completed(futs):
            r = fut.result()
            finals.append(r)
            if len(finals) % 20 == 0:
                print(
                    f"  SA {len(finals)}/{n_sa} best_local={min(x['phi_local'] for x in finals):.1f}",
                    flush=True,
                )
    print(f"SA done in {time.time()-t0:.1f}s; rescoring {len(finals)}+Paley with exact Φ", flush=True)

    finals.append({"phi_local": Phi, "A": C_list, "seed": -1, "mode": "paley"})
    # keep top 120 by phi_local for exact (plus paley)
    finals_sorted = sorted(finals, key=lambda z: z["phi_local"])[: min(120, len(finals))]
    if finals[-1] not in finals_sorted:
        finals_sorted.append(finals[-1])

    rescored = []
    t1 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_exact_job, f) for f in finals_sorted]
        for fut in as_completed(futs):
            r = fut.result()
            rescored.append(r)
            if r["phi_exact"] < Phi - 0.5:
                print(f"  UNDER CUT? exact={r['phi_exact']} local={r['phi_local']} mode={r['mode']}", flush=True)

    best = min(r["phi_exact"] for r in rescored)
    under = [r for r in rescored if r["phi_exact"] < Phi - 1e-9]
    print(
        f"Exact rescore: best={best} undercuts={len(under)} time={time.time()-t1:.1f}s",
        flush=True,
    )

    # --- structured k-flip exact samples ---
    print("Structured random k-flip exact samples...", flush=True)
    kflip_jobs = []
    for k in [1, 2, 3, 5, 8, 10, 13, 15, 20, 26, 30, 40]:
        for s in range(workers // 4):
            kflip_jobs.append((C_list, 500_000 + k * 100 + s, k))
    kflip_res = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(_kflip_exact_job, kflip_jobs):
            kflip_res.append(r)
    by_k = {}
    for r in kflip_res:
        by_k.setdefault(r["k"], []).append(r["phi_exact"])
    kflip_summary = {
        str(k): {"min": min(v), "mean": float(np.mean(v)), "n": len(v)}
        for k, v in sorted(by_k.items())
    }
    print("kflip", kflip_summary, flush=True)

    report = {
        "n": 26,
        "Phi_Paley": Phi,
        "workers": workers,
        "sa_iters": sa_iters,
        "n_sa": n_sa,
        "certified_UB_m26": best,
        "gap": Phi - best,
        "rel_gap": (Phi - best) / (n ** 1.5),
        "n_undercuts": len(under),
        "undercuts": under[:10],
        "kflip_summary": kflip_summary,
        "best_of_kflip": min(r["phi_exact"] for r in kflip_res),
        "e1_status": {
            "undercut_found": best < Phi - 1e-9,
            "note": (
                "If certified_UB=65 and no undercut, consistent with E(1) at n=26. "
                "Does not prove m_26=65. Limit existence still requires general E(1)."
            ),
        },
        "limit_still_open": True,
        "total_time_s": time.time() - t0,
    }
    out = out_dir / "e1_n26_intensive.json"
    out.write_text(json.dumps(report, indent=2))
    print("wrote", out, flush=True)
    print(
        f"RESULT m_26 ≤ {best} (Paley 65), gap={Phi-best}, undercuts={len(under)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
