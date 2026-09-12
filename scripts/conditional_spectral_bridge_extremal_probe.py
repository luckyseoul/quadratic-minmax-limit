#!/usr/bin/env python3
"""Exact finite stress test for the conditional-cross spectral bridge.

For a fixed internal source A and its known conditional optimum T, this
solves, for every Boolean right vector y, the exact problem

  maximize ||B y||_2^2 subject to B in {-1,1}^{n x n} and
  max_{x,z} |q_A(x)-q_A(z)| + |x^T B z| <= T.

The maximum divided by n is a certified lower bound for ||B||_op^2 among
conditional minimizers.  This is a finite obstruction screen only: it
neither bounds every singular direction nor proves an all-orders bridge.
"""
from __future__ import annotations

import argparse
import itertools
import json

import numpy as np
from ortools.sat.python import cp_model


def source(n: int) -> np.ndarray:
    if n == 5:
        a = np.ones((n, n), dtype=np.int8)
        np.fill_diagonal(a, 0)
        for i in range(n):
            a[i, (i + 1) % n] = a[(i + 1) % n, i] = -1
        return a
    if n == 6:
        return np.array(
            [[0, 1, 1, 1, 1, 1], [1, 0, 1, -1, -1, 1],
             [1, 1, 0, 1, -1, -1], [1, -1, 1, 0, 1, -1],
             [1, -1, -1, 1, 0, 1], [1, 1, -1, -1, 1, 0]],
            dtype=np.int8,
        )
    raise ValueError(n)


def solve_direction(a: np.ndarray, x: np.ndarray, q: np.ndarray, target: int,
                    y: np.ndarray, workers: int, seconds: float) -> dict[str, object]:
    n = a.shape[0]
    model = cp_model.CpModel()
    bits = [model.new_bool_var(f"b_{i}_{j}") for i in range(n) for j in range(n)]
    for xp, qp in zip(x, q, strict=True):
        for z, qz in zip(x, q, strict=True):
            linear = sum(int(xp[i] * z[j]) * (2 * bit - 1)
                         for (i, j), bit in zip(itertools.product(range(n), repeat=2), bits,
                                                  strict=True))
            slack = target - abs(int(qp - qz))
            model.add(linear <= slack)
            model.add(-linear <= slack)
    squares: list[cp_model.IntVar] = []
    for i in range(n):
        row = model.new_int_var(-n, n, f"row_{i}")
        model.add(row == sum(int(y[j]) * (2 * bits[i * n + j] - 1) for j in range(n)))
        index = model.new_int_var(0, n, f"index_{i}")
        model.add(row + n == 2 * index)
        square = model.new_int_var(0, n * n, f"square_{i}")
        model.add_element(index, [(2 * k - n) ** 2 for k in range(n + 1)], square)
        squares.append(square)
    energy = model.new_int_var(0, n**3, "energy")
    model.add(energy == sum(squares))
    model.maximize(energy)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.max_time_in_seconds = seconds
    status = solver.solve(model)
    if status != cp_model.OPTIMAL:
        raise RuntimeError(f"expected OPTIMAL, got {solver.status_name(status)}")
    b = np.array([2 * solver.value(bit) - 1 for bit in bits], dtype=np.int8).reshape(n, n)
    cross = np.einsum("pi,ij,qj->pq", x, b, x, optimize=True)
    replay = int(np.max(np.abs(q[:, None] - q[None, :]) + np.abs(cross)))
    if replay != target:
        raise AssertionError((replay, target))
    return {"energy": int(solver.value(energy)), "right_vector": y.tolist(),
            "matrix": b.tolist(), "beta": int(np.max(np.abs(cross))),
            "op_norm": float(np.linalg.svd(b.astype(float), compute_uv=False)[0])}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=(5, 6), default=5)
    parser.add_argument("--workers", type=int, default=88)
    parser.add_argument("--seconds", type=float, default=600)
    args = parser.parse_args()
    a = source(args.n)
    x = np.array(list(itertools.product((-1, 1), repeat=args.n)), dtype=np.int8)
    q = np.einsum("pi,ij,pj->p", x, a, x, optimize=True) // 2
    target = {5: 13, 6: 18}[args.n]
    records = [solve_direction(a, x, q, target, y, args.workers, args.seconds)
               for y in x]
    best = max(records, key=lambda row: int(row["energy"]))
    out = {"classification": "exact finite lower-bound screen; not an all-orders spectral bridge",
           "n": args.n, "source_phi": int(np.max(np.abs(q))),
           "conditional_value": target, "directions_examined": len(records),
           "max_boolean_direction_energy": best["energy"],
           "op_norm_squared_lower_bound": best["energy"] / args.n,
           "maximizing_direction_witness": best, "workers": args.workers}
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
