#!/usr/bin/env python3
"""Search for a p=13,k=7 scalar-7 solution outside one signed-PSL orbit.

The packed forbidden rows are independently extracted from the complete orbit
by ``k7_p13_extract_orbit_representatives.py``.  INFEASIBLE certifies that the
known orbit exhausts this nonzero quintic scalar class; FEASIBLE returns a seed
for a new orbit.  UNKNOWN is only a timed diagnostic.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "evidence"))

from k5_p23_coefficient_sieve import quartic_kernel  # noqa: E402
from k7_p13_cpsat_probe import P, SolutionProbe, build_model  # noqa: E402


def unpack_finite(packed: np.ndarray) -> np.ndarray:
    q = P * P
    indices = np.arange(1, q + 1, dtype=np.int64)
    return (
        (packed[:, indices // 64] >> (indices % 64).astype(np.uint64))
        & np.uint64(1)
    ).astype(np.int8)


def main() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("representatives", type=Path)
    parser.add_argument("--top-scalar", type=int, default=7)
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    packed = np.load(args.representatives)
    forbidden = unpack_finite(packed)
    if len(np.unique(forbidden, axis=0)) != len(forbidden):
        raise RuntimeError("duplicate forbidden representatives")

    started = time.monotonic()
    model, negative, metadata = build_model(args.top_scalar)
    model.add_forbidden_assignments(negative, forbidden.tolist())
    kernel_real, kernel_imag = quartic_kernel(P)
    callback = SolutionProbe(
        negative,
        kernel_real,
        kernel_imag,
        max_solutions=1,
        stop_at_limit=True,
    )
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.time_limit
    solver.parameters.num_search_workers = args.workers
    solver.parameters.random_seed = args.seed
    solver.parameters.randomize_search = True
    status = solver.solve(model, callback)
    report = {
        **metadata,
        "algorithm": "CP-SAT feasibility after forbidding a complete signed-PSL orbit slice",
        "forbidden_representatives": len(forbidden),
        "representatives_path": str(args.representatives),
        "status": solver.status_name(status),
        "solve_seconds": solver.wall_time,
        "elapsed_seconds": time.monotonic() - started,
        "workers": args.workers,
        "random_seed": args.seed,
        "known_orbit_exhausts_scalar_class": status == cp_model.INFEASIBLE,
        "outside_orbit_solution_found": callback.count == 1,
        "outside_orbit_Zpsi": {
            "real": callback.best_real,
            "imag": callback.best_imag,
            "abs_sq": callback.min_abs_z_sq,
            "negative_indices": callback.best_negative_indices,
        },
        "finite_diagnostic_only": status != cp_model.INFEASIBLE,
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == "__main__":
    main()
