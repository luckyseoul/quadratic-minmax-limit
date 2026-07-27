#!/usr/bin/env python3
"""
E(1) fast gap probe: SA with phi_local, exact-Φ rescore of finalists + Paley.

Writes evidence/e1_paley_gap.json. Does not claim m_n equality or lim settlement.
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
    dual_gaussian_lower_bound,
    is_conference_matrix,
    op_norm_sym,
    paley_conference_matrix,
    paley_conference_prime_power,
    phi_local,
    spherical_half_bound,
)
from e1_paley_gap_campaign import exact_phi_gray, maximizer_count_gray

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


def _sa_local(args):
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    n, iters, seed, C_list, mode = args
    rng = np.random.default_rng(seed)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    ne = len(edges)
    if mode == "paley" and C_list is not None:
        A = np.array(C_list, dtype=np.float64)
        for _ in range(int(rng.integers(2, 12))):
            i, j = edges[int(rng.integers(0, ne))]
            _flip(A, i, j)
    else:
        A = np.zeros((n, n))
        for i, j in edges:
            s = float(rng.choice([-1.0, 1.0]))
            A[i, j] = A[j, i] = s
    cur = phi_local(A, restarts=12, rng=rng)
    best = cur
    best_A = A.copy()
    T0 = max(2.0, cur * 0.12)
    for t in range(iters):
        T = T0 * (1.0 - t / iters) + 1e-9
        e = int(rng.integers(0, ne))
        i, j = edges[e]
        _flip(A, i, j)
        new = phi_local(A, restarts=4, rng=rng)
        d = new - cur
        if d <= 0 or rng.random() < math.exp(-d / T):
            cur = new
            if cur < best:
                best = cur
                best_A = A.copy()
        else:
            _flip(A, i, j)
    return {"phi_local": float(best), "A": best_A.tolist(), "seed": seed}


def _exact_one(args):
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    A = np.array(args["A"], dtype=np.float64)
    t0 = time.time()
    ph = exact_phi_gray(A)
    return {
        "seed": args.get("seed", -1),
        "phi_local": args.get("phi_local"),
        "phi_exact": float(ph),
        "time_s": time.time() - t0,
    }


def run_n(n, C, workers, sa_iters):
    Phi, nmax = maximizer_count_gray(C)
    C_list = C.tolist()
    jobs = []
    for w in range(workers):
        mode = "paley" if w % 3 == 0 else "random"
        jobs.append((n, sa_iters, 300000 + n * 1000 + w, C_list, mode))
    t0 = time.time()
    finals = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for fut in as_completed([ex.submit(_sa_local, j) for j in jobs]):
            finals.append(fut.result())
    # rescore all + pure Paley
    finals.append({"phi_local": Phi, "A": C_list, "seed": -1})
    rescored = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for fut in as_completed([ex.submit(_exact_one, f) for f in finals]):
            rescored.append(fut.result())
    best = min(r["phi_exact"] for r in rescored)
    under = sum(1 for r in rescored if r["phi_exact"] < Phi - 1e-9)
    return {
        "n": n,
        "Phi_Paley_exact": Phi,
        "n_maximizers_halfcube": nmax,
        "rho_Paley": 2.0 * Phi / (n * math.sqrt(n - 1)),
        "dual_gauss_LB": dual_gaussian_lower_bound(n),
        "certified_UB_m": best,
        "gap_Phi_minus_UB": Phi - best,
        "rel_gap_over_n32": (Phi - best) / (n ** 1.5),
        "m_over_Phi_UB": best / Phi,
        "alpha_UB": best / (n ** 1.5),
        "n_exact_undercut_paley": under,
        "sa_iters": sa_iters,
        "time_s": time.time() - t0,
        "method": "phi_local_SA_plus_exact_rescore",
    }


def main():
    workers = int(os.environ.get("WORKERS", max(1, (os.cpu_count() or 4) - 2)))
    out_dir = Path(os.environ.get("OUT", str(ROOT / "evidence")))
    out_dir.mkdir(parents=True, exist_ok=True)

    report = {
        "campaign": "E1_paley_gap_fast",
        "workers": workers,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "settlement_claim": False,
        "limit_still_open": True,
        "orders": {},
        "exact_anchors": {"m_6": 5, "m_10": 13},
    }

    builds = [
        (6, lambda: paley_conference_matrix(5), 1500),
        (10, lambda: paley_conference_prime_power(3), 3000),
        (14, lambda: paley_conference_matrix(13), 4000),
        (18, lambda: paley_conference_matrix(17), 5000),
    ]
    for n, build, iters in builds:
        print(f"=== n={n} ===", flush=True)
        C = build()
        assert is_conference_matrix(C)
        row = run_n(n, C, workers, iters)
        report["orders"][str(n)] = row
        print(
            f"  Φ={row['Phi_Paley_exact']} UB≤{row['certified_UB_m']} "
            f"gap={row['gap_Phi_minus_UB']} undercuts={row['n_exact_undercut_paley']}",
            flush=True,
        )

    report["e1_status"] = {
        "proved": False,
        "note": (
            "n=6 exact opt (gap 0); n=10 exact gap 2; n=14,18 SA+exact found no "
            "undercut of Paley in this run (UB=Φ). Supports E(1) locally, not a proof. "
            "lim α_n remains OPEN."
        ),
    }
    out = out_dir / "e1_paley_gap.json"
    out.write_text(json.dumps(report, indent=2))
    print("wrote", out, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
