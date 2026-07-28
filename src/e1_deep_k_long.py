#!/usr/bin/env python3
"""
Longer MILP for deep two-sided feasibility at p=5, k in {14,15,16,18,20}.
Also min max S- for k=12,14,15 with longer time.
File-based ProcessPool. Writes evidence/e1_deep_k_long.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

CACHE = Path("/tmp/e1_deep_p5")


def _deep_worker(k: int, tlim: float = 300.0) -> dict:
    Chip = np.load(CACHE / "chip.npy")
    Chim = np.load(CACHE / "chim.npy")
    N, E = Chip.shape
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
        options={"time_limit": tlim},
    )
    msg = str(res.message)
    return {
        "k": k,
        "ok": bool(res.success),
        "infeasible": (not res.success) and ("infeas" in msg.lower()),
        "timeout": "time limit" in msg.lower(),
        "message": msg[:140],
        "seconds": time.time() - t0,
    }


def _minmax_worker(k: int, tlim: float = 300.0) -> dict:
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
        options={"time_limit": tlim},
    )
    msg = str(res.message)
    return {
        "k": k,
        "success": bool(res.success),
        "min_max_S_minus": float(res.x[-1]) if res.success and res.x is not None else None,
        "timeout": "time limit" in msg.lower(),
        "message": msg[:140],
        "seconds": time.time() - t0,
    }


def main() -> None:
    assert (CACHE / "chip.npy").is_file(), "run e1_deep_cover_hunt.py first"
    workers = min(12, max(1, (os.cpu_count() or 4) - 2))
    t0 = time.time()
    deep_ks = [14, 15, 16, 18, 20]
    min_ks = [12, 14, 15]

    deep = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_deep_worker, k, 240.0): k for k in deep_ks}
        for f in as_completed(futs):
            deep.append(f.result())
            print("deep", deep[-1], flush=True)
    deep.sort(key=lambda r: r["k"])

    mm = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_minmax_worker, k, 240.0): k for k in min_ks}
        for f in as_completed(futs):
            mm.append(f.result())
            print("minmax", mm[-1], flush=True)
    mm.sort(key=lambda r: r["k"])

    out = {
        "p": 5,
        "deep": deep,
        "min_max_S_minus": mm,
        "any_deep_found": any(r["ok"] for r in deep),
        "workers": workers,
        "seconds": time.time() - t0,
    }
    path = ROOT / "evidence" / "e1_deep_k_long.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print("wrote", path)


if __name__ == "__main__":
    main()
