#!/usr/bin/env python3
"""Feasibility B<=34 at n=9. INFEASIBLE 32 already proved min B in {34,36}."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from original_mo_two_half_geometry import (  # noqa: E402
    boolean_states,
    edge_pairs,
    quadratic_values,
)

OUT = ROOT / "evidence" / "two_half_n9_n10_20260919"
A_PATH = OUT / "A9_exact.npy"
HINT = "001000000101101111010111000100001100"
CAP = 34


def main() -> None:
    a = np.load(A_PATH)
    n = a.shape[0]
    edges = edge_pairs(n)
    states = boolean_states(n, fix_first=True)
    q = quadratic_values(a, states)
    model = cp_model.CpModel()
    bits = [model.new_bool_var(f"r_{i}_{j}") for i, j in edges]
    model.add(bits[0] == 0)
    for xi, x in enumerate(states):
        for yi, y in enumerate(states):
            u = abs(int(q[xi] + q[yi]))
            coeffs = [int(x[i] * y[j] - x[j] * y[i]) for i, j in edges]
            expr = sum(c * (2 * bit - 1) for c, bit in zip(coeffs, bits, strict=True))
            model.add(u + expr <= CAP)
            model.add(u - expr <= CAP)
    for bit, value in zip(bits, HINT, strict=True):
        model.add_hint(bit, int(value))
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 86
    solver.parameters.max_time_in_seconds = 600.0
    t0 = time.time()
    print(f"feasibility B<={CAP}, workers=86", flush=True)
    status = solver.solve(model)
    rec = {
        "n": 9,
        "M": 12,
        "cap": CAP,
        "status": solver.status_name(status),
        "seconds": round(time.time() - t0, 3),
        "num_conflicts": int(solver.num_conflicts),
        "num_branches": int(solver.num_branches),
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        rec["witness_bits"] = "".join(str(solver.value(b)) for b in bits)
        rec["min_B"] = CAP
    elif status == cp_model.INFEASIBLE:
        rec["min_B"] = 36
    (OUT / "n9_feas34.json").write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=2), flush=True)


if __name__ == "__main__":
    main()
