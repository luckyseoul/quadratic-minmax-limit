#!/usr/bin/env python3
"""Exact small-order conditional-cross minimax probe for the spectral bridge.

One conditional minimizer is obtained with CP-SAT.  It is finite evidence
only and does not classify all minimizers or prove an asymptotic bridge.
"""
from __future__ import annotations

import argparse
import itertools
import json

import numpy as np
from ortools.sat.python import cp_model


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=(5, 6, 7), default=5)
    parser.add_argument("--workers", type=int, default=88)
    parser.add_argument("--seconds", type=float, default=600)
    args = parser.parse_args()
    n = args.n
    if n == 5:
        a = np.ones((n, n), dtype=np.int8)
        np.fill_diagonal(a, 0)
        for i in range(n):
            a[i, (i + 1) % n] = a[(i + 1) % n, i] = -1
    elif n == 6:
        a = np.array(
            [[0, 1, 1, 1, 1, 1], [1, 0, 1, -1, -1, 1],
             [1, 1, 0, 1, -1, -1], [1, -1, 1, 0, 1, -1],
             [1, -1, -1, 1, 0, 1], [1, 1, -1, -1, 1, 0]],
            dtype=np.int8,
        )
    else:
        a = np.array(
            [[0, 1, 1, 1, 1, 1, 1], [1, 0, -1, 1, 1, 1, 1],
             [1, -1, 0, 1, 1, -1, -1], [1, 1, 1, 0, -1, -1, 1],
             [1, 1, 1, -1, 0, 1, -1], [1, 1, -1, -1, 1, 0, -1],
             [1, 1, -1, 1, -1, -1, 0]],
            dtype=np.int8,
        )
    x = np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.int8)
    q = np.einsum("pi,ij,pj->p", x, a, x, optimize=True) // 2
    model = cp_model.CpModel()
    bits = [model.new_bool_var(f"b_{i}_{j}") for i in range(n) for j in range(n)]
    score = model.new_int_var(0, 100, "score")
    for xp, qp in zip(x, q, strict=True):
        for yq, qq in zip(x, q, strict=True):
            internal = abs(int(qp - qq))
            coeff = [int(xp[i] * yq[j]) for i in range(n) for j in range(n)]
            linear = sum(c * (2 * bit - 1) for c, bit in zip(coeff, bits, strict=True))
            model.add(score >= internal + linear)
            model.add(score >= internal - linear)
    model.minimize(score)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = args.workers
    solver.parameters.max_time_in_seconds = args.seconds
    status = solver.solve(model)
    if status != cp_model.OPTIMAL:
        raise RuntimeError(f"expected OPTIMAL, got {solver.status_name(status)}")
    b = np.array([2 * solver.value(bit) - 1 for bit in bits], dtype=np.int8).reshape(n, n)
    cross = np.einsum("pi,ij,qj->pq", x, b, x, optimize=True)
    direct = int(np.max(np.abs(q[:, None] - q[None, :]) + np.abs(cross)))
    out = {
        "classification": "exact finite n=5 conditional-minimax solve; one minimizer only",
        "n": n,
        "source_phi": int(np.max(np.abs(q))),
        "status": solver.status_name(status),
        "conditional_value": int(solver.value(score)),
        "direct_replay_value": direct,
        "op_norm": float(np.linalg.svd(b.astype(float), compute_uv=False)[0]),
        "beta": int(np.max(np.abs(cross))),
        "matrix": b.tolist(),
        "workers": args.workers,
    }
    assert out["conditional_value"] == out["direct_replay_value"]
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
