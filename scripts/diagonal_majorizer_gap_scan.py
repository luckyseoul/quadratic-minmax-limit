"""Finite canonical-gap diagnostic for a fixed complete source/cross family.

The SDP is numerical.  This script records primal residuals and treats its
output only as a lead for the all-orders small-canonical-gap premise.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import time

import cvxpy as cp
import numpy as np


def source_matrix(n: int) -> np.ndarray:
    # At n=3 every signing is switching-equivalent to this one; the fixed
    # source is retained explicitly rather than inferred from that fact.
    if n == 3:
        a = -np.ones((n, n), dtype=float)
        np.fill_diagonal(a, 0.0)
        return a
    if n == 4:
        return np.array([[0, -1, -1, -1], [-1, 0, -1, -1],
                         [-1, -1, 0, 1], [-1, -1, 1, 0]], dtype=float)
    raise ValueError("only the fixed exact source orders 3 and 4 are supported")


def solve_one(task: tuple[int, int, float]) -> dict:
    n, mask, tol = task
    a = source_matrix(n)
    b = np.ones((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            e = i * n + j
            b[i, j] = -1.0 if (mask >> e) & 1 else 1.0
    k = np.block([[a, b], [b.T, -a]])
    N = 2 * n
    d = cp.Variable(N, nonneg=True)
    D = cp.diag(d)
    problem = cp.Problem(cp.Minimize(cp.sum(d)), [D - k >> 0, D + k >> 0])
    problem.solve(solver="CLARABEL", tol_gap_abs=tol, tol_feas=tol, tol_gap_rel=tol,
                  max_iter=500)
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError(f"mask={mask}: {problem.status}")
    dv = np.asarray(d.value, dtype=float)
    eigminus = float(np.linalg.eigvalsh(np.diag(dv) - k).min())
    eigplus = float(np.linalg.eigvalsh(np.diag(dv) + k).min())
    eig = np.linalg.eigvalsh(k)
    s3 = float(np.abs(eig) ** 3 @ np.ones(N))
    S = float(dv.sum())
    gap = S - s3 / (N - 1)
    return {"mask": mask, "tau": S, "s3_over_q": s3 / (N - 1),
            "gap": gap, "relative_gap": gap / S,
            "min_eig_D_minus_K": eigminus, "min_eig_D_plus_K": eigplus,
            "diagonal": dv.tolist(), "status": problem.status}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=3)
    p.add_argument("--workers", type=int, required=True)
    p.add_argument("--tolerance", type=float, default=1e-8)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise ValueError("Refusing to overwrite output")
    if args.n not in (3, 4):
        raise ValueError("This diagnostic is deliberately limited to n=3,4")
    if not 1 <= args.workers <= len(os.sched_getaffinity(0)):
        raise ValueError("Invalid worker budget")
    total = 1 << (args.n * args.n)
    start = time.monotonic()
    with ProcessPoolExecutor(max_workers=min(args.workers, total)) as pool:
        rows = list(pool.map(solve_one, [(args.n, m, args.tolerance) for m in range(total)]))
    rows.sort(key=lambda x: (x["relative_gap"], x["mask"]))
    # All numerical PSD residuals must be no worse than the solver tolerance
    # by a small eigenvalue-reconstruction allowance.
    assert min(min(r["min_eig_D_minus_K"], r["min_eig_D_plus_K"]) for r in rows) >= -5e-7
    source = source_matrix(args.n)
    payload = {"classification": "numerical finite SDP diagnostic; not a proof or asymptotic estimate",
               "n": args.n, "source": source.astype(int).tolist(), "cross_blocks": total,
               "workers": min(args.workers, total), "tolerance": args.tolerance,
               "seconds": time.monotonic() - start,
               "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "rows": rows}
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"cross_blocks": total, "seconds": payload["seconds"],
                      "best": rows[0], "worst_relative_gap": rows[-1]["relative_gap"]}, indent=2))


if __name__ == "__main__":
    main()
