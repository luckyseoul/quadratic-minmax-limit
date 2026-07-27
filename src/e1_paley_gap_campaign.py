#!/usr/bin/env python3
"""
E(1) campaign: asymptotic optimality of m_n vs Paley Φ along conference orders.

Goal relevance: E(1) along the dense ρ=1 family (or all Paley) + Prop 6.2
would settle L = lim α_n. This script does NOT claim settlement; it produces
certified evidence about gaps Φ(C) - UB(m_n) at small Paley orders and
maximizer-set sizes.

Certified upper bounds on m_n: min exact Φ over candidate matrices
(SA finalists + Paley + structured flips). Exact Φ via Gray code (n≤26 feasible).
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import (
    dual_gaussian_lower_bound,
    form_Q,
    halfspace_boolean_vector,
    is_conference_matrix,
    op_norm_sym,
    paley_conference_matrix,
    paley_conference_prime_power,
    phi,
    phi_local,
    random_sym_pm1,
    spherical_half_bound,
)

for _k in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_k, "1")


def exact_phi_gray(A: np.ndarray) -> float:
    """Exact Φ via Gray code over half-cube. O(n 2^{n-1})."""
    n = A.shape[0]
    x = np.ones(n, dtype=np.float64)
    Ax = A @ x
    q2 = float(x @ Ax)
    best = abs(q2) * 0.5
    gray = 0
    npat = 1 << (n - 1)
    for i in range(1, npat):
        ng = i ^ (i >> 1)
        bit = (gray ^ ng).bit_length() - 1
        gray = ng
        j = bit + 1
        xj = x[j]
        q2 += -4.0 * xj * Ax[j]
        Ax -= 2.0 * xj * A[:, j]
        x[j] = -xj
        aq = abs(q2) * 0.5
        if aq > best:
            best = aq
    return best


def maximizer_count_gray(C: np.ndarray) -> Tuple[float, int]:
    Phi = exact_phi_gray(C)
    n = C.shape[0]
    x = np.ones(n, dtype=np.float64)
    Ax = C @ x
    q2 = float(x @ Ax)
    cnt = 1 if abs(abs(q2) * 0.5 - Phi) < 1e-9 else 0
    gray = 0
    npat = 1 << (n - 1)
    for i in range(1, npat):
        ng = i ^ (i >> 1)
        bit = (gray ^ ng).bit_length() - 1
        gray = ng
        j = bit + 1
        xj = x[j]
        q2 += -4.0 * xj * Ax[j]
        Ax -= 2.0 * xj * C[:, j]
        x[j] = -xj
        if abs(abs(q2) * 0.5 - Phi) < 1e-9:
            cnt += 1
    return Phi, cnt


def _flip(A: np.ndarray, i: int, j: int) -> None:
    A[i, j] *= -1.0
    A[j, i] *= -1.0


def _sa_worker(args) -> Dict:
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    n, iters, seed, C_list, start_mode = args
    rng = np.random.default_rng(seed)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    ne = len(edges)
    if start_mode == "paley" and C_list is not None:
        A = np.array(C_list, dtype=np.float64)
        for _ in range(int(rng.integers(2, 10))):
            i, j = edges[int(rng.integers(0, ne))]
            _flip(A, i, j)
    else:
        A = np.zeros((n, n), dtype=np.float64)
        for i, j in edges:
            s = float(rng.choice([-1.0, 1.0]))
            A[i, j] = A[j, i] = s
    # SA on exact Φ (n≤18: 2^{17}=~1e5, fine)
    cur = exact_phi_gray(A)
    best = cur
    best_A = A.copy()
    T0 = max(2.0, cur * 0.12)
    for t in range(iters):
        T = T0 * (1.0 - t / iters) + 1e-9
        e = int(rng.integers(0, ne))
        i, j = edges[e]
        _flip(A, i, j)
        new = exact_phi_gray(A)
        d = new - cur
        if d <= 0 or rng.random() < math.exp(-d / T):
            cur = new
            if cur < best:
                best = cur
                best_A = A.copy()
        else:
            _flip(A, i, j)
    return {"phi": float(best), "seed": seed, "start": start_mode, "A": best_A.tolist()}


def run_order(
    n: int,
    C: np.ndarray,
    workers: int,
    sa_iters: int,
    n_sa: int,
) -> Dict:
    Phi_C = exact_phi_gray(C)
    dual = dual_gaussian_lower_bound(n)
    half = spherical_half_bound(n)
    rho = 2.0 * Phi_C / (n * math.sqrt(n - 1))
    _, nmax = maximizer_count_gray(C)

    C_list = C.tolist()
    jobs = []
    for w in range(n_sa):
        mode = "paley" if w % 3 == 0 else "random"
        jobs.append((n, sa_iters, 200_000 + n * 1000 + w, C_list, mode))

    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_sa_worker, j) for j in jobs]
        for fut in as_completed(futs):
            results.append(fut.result())

    # also score pure Paley
    paley_phi = Phi_C
    best_ub = min([r["phi"] for r in results] + [paley_phi])
    undercuts = [r for r in results if r["phi"] < Phi_C - 1e-9]

    # fingerprint best
    best_r = min(results, key=lambda z: z["phi"])
    A_best = np.array(best_r["A"], dtype=np.float64)
    op = op_norm_sym(A_best)
    r_prod = 2.0 * best_r["phi"] / (n * math.sqrt(n - 1))

    return {
        "n": n,
        "Phi_Paley_exact": Phi_C,
        "rho_Paley": rho,
        "spherical_half": half,
        "dual_gauss_LB": dual,
        "n_maximizers_halfcube": nmax,
        "certified_UB_m": best_ub,
        "gap_Phi_minus_UB": Phi_C - best_ub,
        "rel_gap_over_n32": (Phi_C - best_ub) / (n ** 1.5),
        "m_over_Phi_UB": best_ub / Phi_C,
        "alpha_UB": best_ub / (n ** 1.5),
        "n_sa_undercut_paley": len(undercuts),
        "sa_iters": sa_iters,
        "n_sa": n_sa,
        "time_s": time.time() - t0,
        "best_fingerprint": {
            "phi": best_r["phi"],
            "op": op,
            "op_over_sqrt_n1": op / math.sqrt(n - 1),
            "r": r_prod,
            "seed": best_r["seed"],
            "start": best_r["start"],
        },
        "note": (
            "certified_UB_m is a valid upper bound on m_n (min exact Φ over candidates). "
            "Does not prove m_n equals this value unless exhaustive."
        ),
    }


def main():
    workers = int(os.environ.get("WORKERS", max(1, (os.cpu_count() or 4) - 2)))
    sa_iters = int(os.environ.get("SA_ITERS", "3000"))
    out_dir = Path(os.environ.get("OUT", str(ROOT / "evidence")))
    out_dir.mkdir(parents=True, exist_ok=True)

    report: Dict = {
        "campaign": "E1_paley_gap",
        "workers": workers,
        "sa_iters": sa_iters,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "settlement_claim": False,
        "orders": {},
    }

    # n=6 Paley q=5, n=10 p=3, n=14 q=13, n=18 q=17
    specs = [
        ("6", lambda: paley_conference_matrix(5), 2000, workers),
        ("10", lambda: paley_conference_prime_power(3), 4000, workers),
        ("14", lambda: paley_conference_matrix(13), sa_iters, workers),
        ("18", lambda: paley_conference_matrix(17), sa_iters, workers),
    ]

    for name, build, iters, w in specs:
        print(f"=== order n={name} ===", flush=True)
        C = build()
        assert is_conference_matrix(C)
        n = C.shape[0]
        row = run_order(n, C, workers=w, sa_iters=iters, n_sa=w)
        report["orders"][name] = row
        print(
            f"  Φ={row['Phi_Paley_exact']} UB_m≤{row['certified_UB_m']} "
            f"gap={row['gap_Phi_minus_UB']} rel={row['rel_gap_over_n32']:.6f} "
            f"undercuts={row['n_sa_undercut_paley']} maxs={row['n_maximizers_halfcube']}",
            flush=True,
        )

    # Known exact anchors
    report["exact_anchors"] = {
        "m_6": 5,
        "m_10": 13,
        "note": "m_6=Φ; m_10=13<15; SA UB should meet these when search is good",
    }
    report["e1_status"] = {
        "proved": False,
        "numeric_support": (
            "At n=14,18 if certified_UB_m stays at Phi_Paley with zero undercuts "
            "under heavy SA, consistent with local E(1) but not a proof. "
            "At n=10 exact gap 2 is o(n^{3/2})."
        ),
        "limit_still_open": True,
    }

    out = out_dir / "e1_paley_gap.json"
    out.write_text(json.dumps(report, indent=2))
    print(f"wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
