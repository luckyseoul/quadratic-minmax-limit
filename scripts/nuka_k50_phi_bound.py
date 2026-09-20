#!/usr/bin/env python3
"""Nuka (5700X3D vcache) CP-SAT bound on Phi of the n=50 plus-I lift.

Q is already sum_{i<j} K_ij X_i X_j (the CORE quadratic, not x^T K x).
Paley C_50 has Phi=175. A proved max |Q| <= 173 is a certified undercut.

Workers default to nproc-2 (14 on nuka). Do not use soulkiller's 86 here.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

ROOT = Path(__file__).resolve().parents[1]
BPATH = ROOT / "evidence" / "coherent_BI_lift_20260919" / "best_n25.json"


def lift(b: np.ndarray) -> np.ndarray:
    d = b + np.eye(b.shape[0])
    return np.block([[b, d], [d, -b]])


def build_q(k: np.ndarray):
    n = k.shape[0]
    model = cp_model.CpModel()
    x = [model.new_bool_var(f"x{i}") for i in range(n)]
    model.add(x[0] == 1)
    qterms = []
    for i in range(n):
        for j in range(i + 1, n):
            p = model.new_bool_var(f"p{i}_{j}")
            model.add_multiplication_equality(p, [x[i], x[j]])
            kij = int(k[i, j])
            qterms.append(kij * (4 * p - 2 * x[i] - 2 * x[j] + 1))
    q = sum(qterms)
    return model, x, q


def solve(model: cp_model.CpModel, workers: int, tlim: float, tag: str) -> dict:
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.max_time_in_seconds = float(tlim)
    solver.parameters.linearization_level = 2
    t0 = time.time()
    print(f"BEGIN {tag} workers={workers} tlim={tlim}", flush=True)
    status = solver.solve(model)
    rec = {
        "tag": tag,
        "status": solver.status_name(status),
        "wall_s": round(time.time() - t0, 3),
        "workers": workers,
        "tlim": tlim,
        "objective": None,
        "bound": None,
        "conflicts": int(solver.num_conflicts),
        "branches": int(solver.num_branches),
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        rec["objective"] = int(solver.objective_value)
        rec["bound"] = int(solver.best_objective_bound)
    elif status == cp_model.UNKNOWN:
        rec["bound"] = int(solver.best_objective_bound)
    print(
        f"END {tag} {rec['status']} obj={rec['objective']} "
        f"bound={rec['bound']} wall={rec['wall_s']} "
        f"conflicts={rec['conflicts']}",
        flush=True,
    )
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 2))
    ap.add_argument("--feas-time", type=float, default=900.0)
    ap.add_argument("--opt-time", type=float, default=2400.0)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    b = np.array(json.loads(BPATH.read_text())["A"], dtype=np.int64)
    k = np.rint(lift(b.astype(float))).astype(np.int64)
    out_dir = ROOT / "evidence" / "k50_nuka_20260919"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out or (out_dir / "phi_bound.json")

    receipt = {
        "host": os.uname().nodename,
        "nproc": os.cpu_count(),
        "workers": args.workers,
        "n": 50,
        "paley": 175,
        "phases": [],
        "undercut_certified": False,
        "note": "Q is CORE Phi (sum i<j). Do not divide by 2.",
    }

    for target in (175, 173, 171):
        model, _x, q = build_q(k)
        model.add(q >= target)
        rec = solve(model, args.workers, args.feas_time, f"feas_Q>={target}")
        rec["target"] = target
        receipt["phases"].append(rec)
        out.write_text(json.dumps(receipt, indent=2))
        if rec["status"] == "INFEASIBLE":
            # Q <= target-2 (Q odd: binom(50,2)=1225 odd, Q ≡ 1225 (mod 2))
            receipt["undercut_certified"] = target <= 175
            receipt["phi_upper"] = target - 2
            print(
                f"INFEASIBLE Q>={target} => Phi <= {target-2} vs Paley 175",
                flush=True,
            )
            break
        if rec["status"] == "FEASIBLE" and target == 175:
            receipt["undercut_certified"] = False
            receipt["phi_lower"] = 175
            print("FEASIBLE Q>=175: this lift does not undercut Paley", flush=True)
            break

    model, _x, q = build_q(k)
    model.maximize(q)
    rec = solve(model, args.workers, args.opt_time, "maximize_Q")
    receipt["phases"].append(rec)
    if rec["objective"] is not None:
        receipt["phi_lower"] = max(receipt.get("phi_lower") or 0, rec["objective"])
    if rec["status"] == "OPTIMAL" and rec["objective"] is not None:
        receipt["phi_exact"] = rec["objective"]
        receipt["undercut_certified"] = rec["objective"] < 175
    out.write_text(json.dumps(receipt, indent=2))
    print("WROTE", out, flush=True)
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()
