#!/usr/bin/env python3
"""Exact max of |Q_B(x)-Q_B(y)+x^T(B+I)y| over high-energy x and all y.

n=25 cube is 2^24. Implication: a witness |Q_K|>=175 shows this lift does
not beat Paley C_50; a certified max <175 on a complete slice is not a
global Phi bound.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BPATH = ROOT / "evidence" / "coherent_BI_lift_20260919" / "best_n25.json"


def cube(m: int) -> np.ndarray:
    nst = 1 << (m - 1)
    idx = np.arange(nst, dtype=np.int64)
    x = np.empty((nst, m), dtype=np.int8)
    x[:, 0] = 1
    for i in range(1, m):
        x[:, i] = (1 - 2 * ((idx >> (i - 1)) & 1)).astype(np.int8)
    return x


def main() -> None:
    b = np.array(json.loads(BPATH.read_text())["A"], dtype=np.float64)
    m = b.shape[0]
    x = cube(m).astype(np.float64)
    q = 0.5 * np.sum((x @ b) * x, axis=1)
    absq = np.abs(q)
    print("B Phi", float(absq.max()), "nstates", len(q), flush=True)
    thr = 60.0
    idx = np.where(absq >= thr - 1e-6)[0]
    print("n high |Q_B|>=", thr, ":", len(idx), flush=True)
    # precompute all y = x (same cube), q_y = q
    best = 0.0
    wmat = x @ (b + np.eye(m))  # (nst, m)
    # scan high x against all y in chunks
    chunk = 4096
    for t, i in enumerate(idx):
        c = q[i]
        w = wmat[i]
        # |c - q[j] + w·x[j]|
        dots = x @ w
        val = np.max(np.abs(c - q + dots))
        if val > best:
            best = float(val)
            print("best", best, "after", t + 1, "high-x", flush=True)
        if best >= 175:
            print("WITNESS >= Paley50", flush=True)
            break
    print("FINAL best over high-x x all-y", best, flush=True)


if __name__ == "__main__":
    main()
