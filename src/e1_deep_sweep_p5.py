#!/usr/bin/env python3
"""Sweep deep two-sided integral feasibility for k=10..40 at p=5."""
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
CACHE = Path("/tmp/e1_deep_p5")


def worker(k: int) -> dict:
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
        options={"time_limit": 90},
    )
    msg = str(res.message)
    return {
        "k": int(k),
        "ok": bool(res.success),
        "infeas": (not res.success) and ("infeas" in msg.lower()),
        "timeout": "time limit" in msg.lower(),
        "sec": time.time() - t0,
        "msg": msg[:100],
    }


def main() -> None:
    ks = list(range(10, 41))
    workers = min(30, max(1, (os.cpu_count() or 4) - 2))
    results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(worker, k): k for k in ks}
        for f in as_completed(futs):
            r = f.result()
            results.append(r)
            print(r, flush=True)
    results.sort(key=lambda r: r["k"])
    out = {
        "p": 5,
        "results": results,
        "infeasible_ks": [r["k"] for r in results if r["infeas"]],
        "feasible_ks": [r["k"] for r in results if r["ok"]],
        "timeout_ks": [r["k"] for r in results if r["timeout"]],
    }
    path = ROOT / "evidence" / "e1_deep_sweep_p5.json"
    path.write_text(json.dumps(out, indent=2))
    print("INFEAS", out["infeasible_ks"])
    print("FEAS", out["feasible_ks"])
    print("TIMEOUT", out["timeout_ks"])
    print("wrote", path)


if __name__ == "__main__":
    main()
