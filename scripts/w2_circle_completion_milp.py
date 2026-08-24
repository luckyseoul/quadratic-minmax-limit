#!/usr/bin/env python3
"""Exact MILP search for the Frobenius-pair circle completion needed by W2.

For the standard nonsquare circle S through infinity, find a sign vector y
such that

    C y = -p y,
    y|_S is the sparse signed -p eigenvector on S, and
    C_ij y_i y_j = -1

for one Frobenius-conjugate pair {i,j} outside S.  The last orbit is the
only outside-pair orbit not supplied uniformly by the named halfspace
completion when p == 3 (mod 4).

Writing y=1-2b makes this a binary linear feasibility problem.  Redundant
real equations are removed by pivoted QR before SciPy/HiGHS is called.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15613 import _finv  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from w2_circle_completion_probe import signed_sparse_eigenvector  # noqa: E402


def conjugate_pair(p: int, circle: np.ndarray):
    """First pair outside S conjugate under S's intrinsic F_p-structure."""
    q, mul, _add, chi, frob, _norm, _ia, _ib = field_ctx(p)
    sigma = next(x for x in range(1, q) if chi(x) == -1)
    sinv = _finv(mul, q, sigma)
    outside = set(int(i) for i in np.flatnonzero(circle == 0))
    for i in sorted(outside):
        if i == 0:
            continue
        u = mul(sinv, i - 1)
        j = 1 + mul(sigma, frob(u))
        if j in outside and i < j:
            return i, int(j)
    raise RuntimeError("no outside Frobenius pair")


def independent_rows(A: np.ndarray):
    from scipy.linalg import qr

    _q, r, piv = qr(A.T, mode="economic", pivoting=True)
    diag = np.abs(np.diag(r))
    tol = max(A.shape) * np.finfo(float).eps * (diag.max() if len(diag) else 0.0)
    rank = int(np.count_nonzero(diag > tol))
    return np.asarray(piv[:rank], dtype=np.int64), rank


def solve_orientation(
    p: int,
    C: np.ndarray,
    v: np.ndarray,
    pair: tuple[int, int],
    first_pair_sign: int,
    time_limit: float | None,
):
    from scipy.optimize import Bounds, LinearConstraint, milp

    n = len(v)
    A = C.astype(np.float64) + p * np.eye(n)
    rows, rank = independent_rows(A)
    A = A[rows]
    rhs = (A @ np.ones(n)) / 2.0

    fixed_y = {int(k): int(v[k]) for k in np.flatnonzero(v)}
    i, j = pair
    fixed_y[i] = first_pair_sign
    fixed_y[j] = -int(C[i, j]) * first_pair_sign
    fixed_b = {k: (1 - y) // 2 for k, y in fixed_y.items()}
    fixed = np.asarray(sorted(fixed_b), dtype=np.int64)
    free = np.asarray([k for k in range(n) if k not in fixed_b], dtype=np.int64)
    bfix = np.asarray([fixed_b[int(k)] for k in fixed], dtype=np.float64)
    reduced_rhs = rhs - A[:, fixed] @ bfix
    reduced_A = A[:, free]

    options = {"presolve": True}
    if time_limit is not None:
        options["time_limit"] = float(time_limit)
    t0 = time.time()
    result = milp(
        np.zeros(len(free)),
        integrality=np.ones(len(free), dtype=np.int8),
        bounds=Bounds(np.zeros(len(free)), np.ones(len(free))),
        constraints=LinearConstraint(reduced_A, reduced_rhs, reduced_rhs),
        options=options,
    )
    elapsed = time.time() - t0
    node_count = getattr(result, "mip_node_count", None)
    mip_gap = getattr(result, "mip_gap", None)
    rec = {
        "backend": "highs",
        "first_pair_sign": first_pair_sign,
        "status": int(result.status),
        "success": bool(result.success),
        "message": str(result.message),
        "equation_rank": rank,
        "n_fixed": int(len(fixed)),
        "n_free": int(len(free)),
        "mip_node_count": None if node_count is None else int(node_count),
        "mip_gap": None if mip_gap is None else float(mip_gap),
        "seconds": round(elapsed, 3),
    }
    if not result.success:
        return rec, None
    b = np.zeros(n, dtype=np.int8)
    b[fixed] = bfix.astype(np.int8)
    b[free] = np.rint(result.x).astype(np.int8)
    y = (1 - 2 * b).astype(np.int8)
    exact = np.array_equal(C.astype(np.int64) @ y.astype(np.int64), -p * y)
    circle_ok = np.array_equal(y[np.flatnonzero(v)], v[np.flatnonzero(v)])
    pair_ok = int(C[i, j]) * int(y[i]) * int(y[j]) == -1
    rec.update(
        {
            "exact_eigenvector": bool(exact),
            "circle_signs_fixed": bool(circle_ok),
            "pair_in_U": bool(pair_ok),
            "weight": int(((1 - y) // 2).sum()),
            "sign_sum": int(y.sum()),
        }
    )
    if not (exact and circle_ok and pair_ok):
        raise RuntimeError("MILP returned a solution that failed exact validation")
    return rec, y


def solve_orientation_cpsat(
    p: int,
    C: np.ndarray,
    v: np.ndarray,
    pair: tuple[int, int],
    first_pair_sign: int,
    time_limit: float | None,
    workers: int,
):
    from ortools.sat.python import cp_model

    n = len(v)
    A_full = C.astype(np.int64) + p * np.eye(n, dtype=np.int64)
    rows, rank = independent_rows(A_full.astype(np.float64))
    A = A_full[rows]
    rhs = (A @ np.ones(n, dtype=np.int64)) // 2

    fixed_y = {int(k): int(v[k]) for k in np.flatnonzero(v)}
    i, j = pair
    fixed_y[i] = first_pair_sign
    fixed_y[j] = -int(C[i, j]) * first_pair_sign
    fixed_b = {k: (1 - y) // 2 for k, y in fixed_y.items()}
    fixed = np.asarray(sorted(fixed_b), dtype=np.int64)
    free = np.asarray([k for k in range(n) if k not in fixed_b], dtype=np.int64)
    bfix = np.asarray([fixed_b[int(k)] for k in fixed], dtype=np.int64)
    reduced_rhs = rhs - A[:, fixed] @ bfix
    reduced_A = A[:, free]

    model = cp_model.CpModel()
    variables = [model.new_bool_var(f"b_{int(k)}") for k in free]
    for row, value in zip(reduced_A, reduced_rhs):
        terms = [int(a) * x for a, x in zip(row, variables) if a]
        model.add(sum(terms) == int(value))
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    if time_limit is not None:
        solver.parameters.max_time_in_seconds = float(time_limit)
    t0 = time.time()
    status = solver.solve(model)
    elapsed = time.time() - t0
    success = status in (cp_model.OPTIMAL, cp_model.FEASIBLE)
    rec = {
        "backend": "cp-sat",
        "first_pair_sign": first_pair_sign,
        "status": int(status),
        "status_name": solver.status_name(status),
        "success": bool(success),
        "equation_rank": rank,
        "n_fixed": int(len(fixed)),
        "n_free": int(len(free)),
        "workers": workers,
        "conflicts": int(solver.num_conflicts),
        "branches": int(solver.num_branches),
        "wall_time": float(solver.wall_time),
        "seconds": round(elapsed, 3),
    }
    if not success:
        return rec, None
    b = np.zeros(n, dtype=np.int8)
    b[fixed] = bfix.astype(np.int8)
    b[free] = np.asarray([solver.value(x) for x in variables], dtype=np.int8)
    y = (1 - 2 * b).astype(np.int8)
    exact = np.array_equal(C.astype(np.int64) @ y.astype(np.int64), -p * y)
    circle_ok = np.array_equal(y[np.flatnonzero(v)], v[np.flatnonzero(v)])
    pair_ok = int(C[i, j]) * int(y[i]) * int(y[j]) == -1
    rec.update(
        {
            "exact_eigenvector": bool(exact),
            "circle_signs_fixed": bool(circle_ok),
            "pair_in_U": bool(pair_ok),
            "weight": int(((1 - y) // 2).sum()),
            "sign_sum": int(y.sum()),
        }
    )
    if not (exact and circle_ok and pair_ok):
        raise RuntimeError("CP-SAT returned a solution that failed exact validation")
    return rec, y


def run(
    p: int,
    time_limit: float | None = None,
    pair_sign: int | None = None,
    backend: str = "highs",
    workers: int = 1,
) -> dict:
    t0 = time.time()
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    circle, v = signed_sparse_eigenvector(p, C)
    pair = conjugate_pair(p, circle)
    attempts = []
    witness = None
    signs = (pair_sign,) if pair_sign is not None else (1, -1)
    for sign in signs:
        if backend == "cp-sat":
            rec, y = solve_orientation_cpsat(
                p, C, v, pair, sign, time_limit, workers
            )
        else:
            rec, y = solve_orientation(p, C, v, pair, sign, time_limit)
        attempts.append(rec)
        print(json.dumps(rec, indent=2), flush=True)
        if y is not None:
            witness = y
            break
    proven_infeasible = pair_sign is None and len(attempts) == 2 and all(
        (a.get("status_name") == "INFEASIBLE")
        or (a.get("status") == 2 and a.get("backend") != "cp-sat")
        for a in attempts
    )
    out = {
        "p": p,
        "n": p * p + 1,
        "circle_support": np.flatnonzero(circle).tolist(),
        "frobenius_pair": list(pair),
        "attempts": attempts,
        "feasible": witness is not None,
        "proven_infeasible": proven_infeasible,
        "solve_state": (
            "feasible"
            if witness is not None
            else "infeasible"
            if proven_infeasible
            else "unknown"
        ),
        "seconds": round(time.time() - t0, 3),
    }
    if witness is not None:
        out["witness_bits_hex"] = np.packbits(
            ((1 - witness) // 2).astype(np.uint8), bitorder="little"
        ).tobytes().hex()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("p", type=int)
    ap.add_argument("--time-limit", type=float)
    ap.add_argument("--pair-sign", type=int, choices=(-1, 1))
    ap.add_argument("--backend", choices=("highs", "cp-sat"), default="highs")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = run(
        args.p,
        args.time_limit,
        args.pair_sign,
        args.backend,
        args.workers,
    )
    print(json.dumps(out, indent=2), flush=True)
    if args.output:
        args.output.write_text(json.dumps(out, indent=2) + "\n")
        print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
