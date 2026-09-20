#!/usr/bin/env python3
"""CP-SAT exact max of Q_K on the n=50 plus-I lift. 50 binary variables."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

ROOT = Path(__file__).resolve().parents[1]
BPATH = ROOT / "evidence" / "coherent_BI_lift_20260919" / "best_n25.json"


def lift(b: np.ndarray) -> np.ndarray:
    d = b + np.eye(b.shape[0])
    return np.block([[b, d], [d, -b]])


def main() -> None:
    b = np.array(json.loads(BPATH.read_text())["A"], dtype=np.int64)
    k = np.rint(lift(b.astype(float))).astype(np.int64)
    n = k.shape[0]
    model = cp_model.CpModel()
    x = [model.new_bool_var(f"x{i}") for i in range(n)]
    model.add(x[0] == 1)
    # Q = sum_{i<j} K_ij X_i X_j with X=2x-1
    # X_i X_j = 4 x_i x_j - 2x_i - 2x_j + 1
    prod = {}
    qterms = []
    for i in range(n):
        for j in range(i + 1, n):
            p = model.new_bool_var(f"p{i}_{j}")
            model.add_multiplication_equality(p, [x[i], x[j]])
            prod[i, j] = p
            kij = int(k[i, j])
            # contribution kij * X_i X_j = kij * (4p - 2x_i - 2x_j + 1)
            qterms.append(kij * (4 * p - 2 * x[i] - 2 * x[j] + 1))
    q = sum(qterms)
    model.maximize(q)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 86
    solver.parameters.max_time_in_seconds = 300.0
    print("CP-SAT max Q, n=50, workers=86, 300s", flush=True)
    status = solver.solve(model)
    print("status", solver.status_name(status), flush=True)
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        qv = solver.objective_value
        print("max Q", qv, "Phi>=", abs(qv) / 2, "bound", solver.best_objective_bound / 2)
        print("Paley50", 175.0)


if __name__ == "__main__":
    main()
