#!/usr/bin/env python3
"""Local-search Phi on K=[[B,B+I],[B+I,-B]] at n=50 (B the n=25 record)."""
from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from minmax_quadratic import phi_local

ROOT = Path(__file__).resolve().parents[1]
# extracted campaign winner; fall back if already copied
CANDIDATES = [
    Path("/tmp/original-mo-cpu-structured.Z49T3i/final_winners/best_n25.json"),
    ROOT / "evidence" / "coherent_BI_lift_20260919" / "best_n25.json",
]


def load_B25() -> np.ndarray:
    for p in CANDIDATES:
        if p.exists():
            return np.array(json.loads(p.read_text())["A"], dtype=np.float64)
    raise FileNotFoundError("best_n25.json")


def lift(B: np.ndarray) -> np.ndarray:
    d = B + np.eye(B.shape[0])
    return np.block([[B, d], [d, -B]])


def one(payload: tuple) -> float:
    seed, klist = payload
    k = np.asarray(klist, dtype=np.float64)
    return phi_local(k, restarts=50, rng=np.random.default_rng(seed))


def main() -> None:
    b = load_B25()
    k = lift(b)
    klist = k.tolist()
    best = 0.0
    with ProcessPoolExecutor(max_workers=86) as ex:
        futs = [ex.submit(one, (s, klist)) for s in range(86)]
        for fut in as_completed(futs):
            v = float(fut.result())
            if v > best:
                best = v
                print("best", best, flush=True)
    print("FINAL", best, "Paley50", 175.0, flush=True)


if __name__ == "__main__":
    main()
