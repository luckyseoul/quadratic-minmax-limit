#!/usr/bin/env python3
"""2-opt local search for Phi of the n=50 plus-I lift. ProcessPool 86."""
from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BPATH = ROOT / "evidence" / "coherent_BI_lift_20260919" / "best_n25.json"


def lift(b: np.ndarray) -> np.ndarray:
    d = b + np.eye(b.shape[0])
    return np.block([[b, d], [d, -b]])


def climb(k: np.ndarray, rng: np.random.Generator, pair_rounds: int = 8) -> float:
    n = k.shape[0]
    x = rng.choice([-1.0, 1.0], size=n)
    kx = k @ x
    q = float(x @ kx)
    # 1-opt to local peak of |Q|
    improved = True
    while improved:
        improved = False
        for i in range(n):
            delta = -4.0 * x[i] * kx[i]
            if abs(q + delta) > abs(q) + 1e-12:
                xi = x[i]
                x[i] = -xi
                kx -= 2.0 * xi * k[:, i]
                q += delta
                improved = True
    # a few random-accept 2-opt probes
    for _ in range(pair_rounds):
        i, j = int(rng.integers(0, n)), int(rng.integers(0, n - 1))
        j = j if j < i else j + 1
        # flip i then j
        d1 = -4.0 * x[i] * kx[i]
        xi = x[i]
        x[i] = -xi
        kx -= 2.0 * xi * k[:, i]
        q1 = q + d1
        d2 = -4.0 * x[j] * kx[j]
        if abs(q1 + d2) > abs(q) + 1e-12:
            xj = x[j]
            x[j] = -xj
            kx -= 2.0 * xj * k[:, j]
            q = q1 + d2
            # re-1-opt
            improved = True
            while improved:
                improved = False
                for t in range(n):
                    delta = -4.0 * x[t] * kx[t]
                    if abs(q + delta) > abs(q) + 1e-12:
                        xt = x[t]
                        x[t] = -xt
                        kx -= 2.0 * xt * k[:, t]
                        q += delta
                        improved = True
        else:
            # revert i
            x[i] = xi
            kx += 2.0 * xi * k[:, i]
    return abs(q) / 2.0


def one(seed: int) -> float:
    b = np.array(json.loads(BPATH.read_text())["A"], dtype=np.float64)
    k = lift(b)
    rng = np.random.default_rng(seed)
    best = 0.0
    for _ in range(30):
        best = max(best, climb(k, rng))
    return best


def main() -> None:
    best = 0.0
    with ProcessPoolExecutor(max_workers=86) as ex:
        futs = [ex.submit(one, s) for s in range(86)]
        for fut in as_completed(futs):
            v = float(fut.result())
            if v > best:
                best = v
                print("best", best, flush=True)
    print("FINAL", best, "Paley50", 175.0, flush=True)


if __name__ == "__main__":
    main()
