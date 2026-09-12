#!/usr/bin/env python3
"""Exact order-four screen for the conditional-cross spectral bridge.

This is deliberately a finite diagnostic, not an all-orders claim.  It
enumerates every 4-by-4 cross signing for one order-four optimal internal
source and records the operator norms of *all* conditional minimizers.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np


def signs(n: int) -> np.ndarray:
    return np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.int8)


def main() -> None:
    n = 4
    # This is the order-four Phi=4 source used by diagonal_majorizer_gap_scan.
    a = np.array(
        [[0, -1, -1, -1], [-1, 0, -1, -1], [-1, -1, 0, 1], [-1, -1, 1, 0]],
        dtype=np.int8,
    )
    x = signs(n)
    q = np.einsum("pi,ij,pj->p", x, a, x, optimize=True) // 2
    internal = np.abs(q[:, None] - q[None, :])
    records: list[dict[str, object]] = []
    best = None
    hist: dict[str, int] = {}
    for start in range(0, 1 << (n * n), 512):
        values = np.arange(start, min(start + 512, 1 << (n * n)), dtype=np.uint32)
        bits = ((values[:, None] >> np.arange(n * n)) & 1).astype(np.int8)
        b = 2 * bits.reshape(-1, n, n) - 1
        cross = np.einsum("pi,bij,qj->bpq", x, b, x, optimize=True)
        scores = np.max(internal[None, :, :] + np.abs(cross), axis=(1, 2))
        local = int(scores.min())
        if best is None or local < best:
            best = local
            records = []
            hist = {}
        if local == best:
            for idx in np.flatnonzero(scores == best):
                bb = b[idx].astype(float)
                op = float(np.linalg.svd(bb, compute_uv=False)[0])
                beta = int(np.max(np.abs(cross[idx])))
                key = f"{op:.12g}"
                hist[key] = hist.get(key, 0) + 1
                records.append({"mask": int(values[idx]), "op_norm": op, "beta": beta})
    assert best is not None
    out = {
        "classification": "exhaustive finite diagnostic only; not an all-orders spectral bridge",
        "n": n,
        "source_phi": int(np.max(np.abs(q))),
        "cross_signings_examined": 1 << (n * n),
        "conditional_value": best,
        "conditional_minimizer_count": len(records),
        "operator_norm_histogram": hist,
        "minimizers": records,
    }
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
