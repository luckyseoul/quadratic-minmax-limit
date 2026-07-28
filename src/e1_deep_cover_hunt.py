#!/usr/bin/env python3
"""
Hunt integral deep two-sided Max-covers at p=5:
  min_{Max+} S >= 2 and max_{Max-} S <= -2, various k.

Also minimize max Max- S subject to min Max+ S >= 2 (master-style).

Writes evidence/e1_deep_cover_hunt.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from e1_bitight_infeas import boolean_evecs  # noqa: E402

CACHE = Path("/tmp/e1_deep_p5")


def _prepare(p: int = 5) -> dict:
    CACHE.mkdir(exist_ok=True)
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    mp_path = CACHE / "maxplus.npy"
    mm_path = CACHE / "maxminus.npy"
    if mp_path.is_file() and mm_path.is_file():
        Mp = np.load(mp_path).astype(float)
        Mm = np.load(mm_path).astype(float)
    else:
        Mp = boolean_evecs(C, float(p), seed=0)
        Mm = boolean_evecs(-C, float(p), seed=1)
        np.save(mp_path, Mp)
        np.save(mm_path, Mm)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    Chip = np.array(
        [[C[a, b] * y[a] * y[b] for a, b in edges] for y in Mp], dtype=np.float64
    )
    Chim = np.array(
        [[C[a, b] * z[a] * z[b] for a, b in edges] for z in Mm], dtype=np.float64
    )
    np.save(CACHE / "chip.npy", Chip)
    np.save(CACHE / "chim.npy", Chim)
    return {
        "p": p,
        "n": n,
        "n_max": int(len(Mp)),
        "E": int(len(edges)),
        "Chip": Chip,
        "Chim": Chim,
    }


def _deep_feas_worker(k: int, time_limit: float = 90.0) -> dict:
    Chip = np.load(CACHE / "chip.npy")
    Chim = np.load(CACHE / "chim.npy")
    N, E = Chip.shape
    # Chip x >= 2, Chim x <= -2, sum x = k
    A_ub = np.vstack([-Chip, Chim])
    b_ub = np.concatenate([-2.0 * np.ones(N), -2.0 * np.ones(N)])
    cons = [
        LinearConstraint(A_ub, -np.inf * np.ones(2 * N), b_ub),
        LinearConstraint(np.ones((1, E)), [float(k)], [float(k)]),
    ]
    t0 = time.time()
    res = milp(
        c=np.zeros(E),
        constraints=cons,
        integrality=np.ones(E),
        bounds=Bounds(0, 1),
        options={"time_limit": time_limit},
    )
    return {
        "k": k,
        "integral_deep_feasible": bool(res.success),
        "message": str(res.message)[:120],
        "seconds": time.time() - t0,
    }


def _min_max_sm_worker(k: int, time_limit: float = 90.0) -> dict:
    """min t s.t. Chim x <= t, Chip x >= 2, sum x = k, x binary."""
    Chip = np.load(CACHE / "chip.npy")
    Chim = np.load(CACHE / "chim.npy")
    N, E = Chip.shape
    c = np.zeros(E + 1)
    c[-1] = 1.0
    A1 = np.hstack([-Chip, np.zeros((N, 1))])
    A2 = np.hstack([Chim, -np.ones((N, 1))])
    A3 = np.hstack([np.ones((1, E)), [[0.0]]])
    cons = [
        LinearConstraint(A1, -np.inf * np.ones(N), -2.0 * np.ones(N)),
        LinearConstraint(A2, -np.inf * np.ones(N), np.zeros(N)),
        LinearConstraint(A3, [float(k)], [float(k)]),
    ]
    integ = np.ones(E + 1)
    integ[-1] = 0
    t0 = time.time()
    res = milp(
        c=c,
        constraints=cons,
        integrality=integ,
        bounds=Bounds([0.0] * E + [-40.0], [1.0] * E + [40.0]),
        options={"time_limit": time_limit},
    )
    tval = float(res.x[-1]) if res.success and res.x is not None else None
    return {
        "k": k,
        "success": bool(res.success),
        "min_max_S_minus": tval,
        "message": str(res.message)[:120],
        "seconds": time.time() - t0,
    }


def main() -> None:
    t0 = time.time()
    meta = _prepare(5)
    Chip, Chim = meta["Chip"], meta["Chim"]
    N, E = Chip.shape

    # fractional deep feasibility for a few k
    frac = []
    for k in range(10, 41, 5):
        A_ub = np.vstack([-Chip, Chim])
        b_ub = np.concatenate([-2.0 * np.ones(N), -2.0 * np.ones(N)])
        res = linprog(
            np.zeros(E),
            A_ub=A_ub,
            b_ub=b_ub,
            A_eq=np.ones((1, E)),
            b_eq=[float(k)],
            bounds=(0, 1),
            method="highs",
        )
        frac.append({"k": k, "fractional_deep": bool(res.success)})

    workers = min(20, max(1, (os.cpu_count() or 4) - 2))
    ks_deep = list(range(10, 31, 2))
    ks_min = [10, 12, 14, 15, 16, 20, 25, 30]

    deep_results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_deep_feas_worker, k, 75.0): k for k in ks_deep}
        for f in as_completed(futs):
            deep_results.append(f.result())
            print("deep", deep_results[-1], flush=True)
    deep_results.sort(key=lambda r: r["k"])

    minmax_results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_min_max_sm_worker, k, 75.0): k for k in ks_min}
        for f in as_completed(futs):
            minmax_results.append(f.result())
            print("minmax", minmax_results[-1], flush=True)
    minmax_results.sort(key=lambda r: r["k"])

    any_deep = any(r["integral_deep_feasible"] for r in deep_results)
    min_of_minmax = [
        r["min_max_S_minus"]
        for r in minmax_results
        if r["success"] and r["min_max_S_minus"] is not None
    ]

    out = {
        "p": 5,
        "n": meta["n"],
        "n_Max": meta["n_max"],
        "fractional_deep": frac,
        "integral_deep_feas": deep_results,
        "any_integral_deep_found": any_deep,
        "min_max_S_minus": minmax_results,
        "all_min_max_S_minus_ge_0": bool(min_of_minmax) and all(t >= -1e-9 for t in min_of_minmax),
        "workers": workers,
        "seconds": time.time() - t0,
        "status": (
            "Integral deep two-sided cover search at p=5; "
            "if none and min max S- >= 0 for s+>=2 covers, deep undercutters impossible at p=5."
        ),
    }
    path = ROOT / "evidence" / "e1_deep_cover_hunt.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps({k: out[k] for k in out if k not in ()}, indent=2))
    # Non-fatal: search may time out without proving infeasibility
    print("wrote", path)


if __name__ == "__main__":
    main()
