#!/usr/bin/env python3
"""Exact Peisert(121) conference matrix and Boolean quadratic optimization."""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
from ortools.sat.python import cp_model


P = 11


def add(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def neg(x: tuple[int, int]) -> tuple[int, int]:
    return ((-x[0]) % P, (-x[1]) % P)


def mul(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    # alpha^2 = 2, irreducible because 2 is nonsquare modulo 11.
    return ((x[0] * y[0] + 2 * x[1] * y[1]) % P,
            (x[0] * y[1] + x[1] * y[0]) % P)


def pow_ff(x: tuple[int, int], e: int) -> tuple[int, int]:
    y = (1, 0)
    while e:
        if e & 1:
            y = mul(y, x)
        x = mul(x, x)
        e //= 2
    return y


def order(x: tuple[int, int]) -> int:
    if x == (0, 0):
        return 0
    y = (1, 0)
    for k in range(1, 121):
        y = mul(y, x)
        if y == (1, 0):
            return k
    raise AssertionError


def construct() -> tuple[np.ndarray, dict]:
    pts = [(a, b) for a in range(P) for b in range(P)]
    idx = {x: i for i, x in enumerate(pts)}
    primitive = next(x for x in pts[1:] if order(x) == 120)
    powers = [pow_ff(primitive, j) for j in range(120)]
    connection = {powers[j] for j in range(120) if j % 4 in (0, 1)}
    assert len(connection) == 60
    assert {neg(x) for x in connection} == connection

    A = np.zeros((121, 121), dtype=np.int64)
    for i, x in enumerate(pts):
        for d in connection:
            A[i, idx[add(x, d)]] = 1
    assert np.array_equal(A, A.T)
    assert np.all(np.diag(A) == 0)
    degrees = A.sum(axis=1)
    common = A @ A
    off_adj = common[(A == 1)]
    off_non = common[(A == 0) & (~np.eye(121, dtype=bool))]
    assert np.all(degrees == 60)
    assert np.all(off_adj == 29)
    assert np.all(off_non == 30)

    S = np.ones((121, 121), dtype=np.int64) - np.eye(121, dtype=np.int64) - 2 * A
    C = np.zeros((122, 122), dtype=np.int64)
    C[0, 1:] = 1
    C[1:, 0] = 1
    C[1:, 1:] = S
    assert np.array_equal(C, C.T)
    assert np.all(np.diag(C) == 0)
    assert np.array_equal(C @ C, 121 * np.eye(122, dtype=np.int64))
    meta = {
        "field_polynomial": "alpha^2-2 over F_11",
        "primitive": list(primitive),
        "connection_exponents_mod_4": [0, 1],
        "srg": [121, 60, 29, 30],
        "conference_order": 122,
        "conference_square_verified": True,
    }
    return C, meta


def score(C: np.ndarray, bits: list[int]) -> int:
    x = 1 - 2 * np.asarray(bits, dtype=np.int64)
    return int(x @ C @ x // 2)


def solve(C: np.ndarray, mode: str, seconds: float, workers: int) -> dict:
    n = len(C)
    model = cp_model.CpModel()
    b = [model.NewBoolVar(f"b{i}") for i in range(n)]
    # Global x -> -x symmetry fixes the normalized conference coordinate.
    # Do not fix a second coordinate: switched candidates need not retain
    # any translation automorphism.
    model.Add(b[0] == 0)

    if mode in ("eigplus", "eigminus"):
        lam = 11 if mode == "eigplus" else -11
        # C(1-2b)=lam(1-2b), exact integer linear equations.
        rowsums = C.sum(axis=1)
        for i in range(n):
            model.Add(sum(int(C[i, j]) * b[j] for j in range(n))
                      == (int(rowsums[i]) - lam) // 2 + lam * b[i])
    else:
        terms = []
        for i in range(n):
            for j in range(i + 1, n):
                z = model.NewBoolVar(f"z{i}_{j}")
                model.Add(z >= b[i] - b[j])
                model.Add(z >= b[j] - b[i])
                model.Add(z <= b[i] + b[j])
                model.Add(z <= 2 - b[i] - b[j])
                terms.append(int(C[i, j]) * z)
        total_upper = int(np.triu(C, 1).sum())
        q = total_upper - 2 * sum(terms)
        if mode == "q669":
            model.Add(q >= 669)
        elif mode == "opt":
            model.Maximize(q)
        else:
            raise ValueError(mode)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.log_search_progress = False
    solver.parameters.symmetry_level = 3
    solver.parameters.cp_model_presolve = True
    t0 = time.time()
    status = solver.Solve(model)
    elapsed = time.time() - t0
    out = {
        "mode": mode,
        "status": solver.StatusName(status),
        "elapsed_seconds": elapsed,
        "response_stats": solver.ResponseStats(),
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        bits = [solver.Value(v) for v in b]
        out["bits"] = bits
        out["score"] = score(C, bits)
        x = 1 - 2 * np.asarray(bits, dtype=np.int64)
        out["conference_residual_plus_maxabs"] = int(np.max(np.abs(C @ x - 11*x)))
        out["conference_residual_minus_maxabs"] = int(np.max(np.abs(C @ x + 11*x)))
        if mode == "opt":
            out["objective"] = solver.ObjectiveValue()
            out["best_bound"] = solver.BestObjectiveBound()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["construct", "eigplus", "eigminus", "q669", "opt"])
    ap.add_argument("--seconds", type=float, default=3600)
    ap.add_argument("--workers", type=int, default=44)
    args = ap.parse_args()
    C, meta = construct()
    out = {"construction": meta}
    if args.mode != "construct":
        out["solve"] = solve(C, args.mode, args.seconds, args.workers)
    print("RESULT_JSON=" + json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()
