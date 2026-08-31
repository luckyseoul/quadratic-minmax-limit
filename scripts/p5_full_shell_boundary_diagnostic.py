#!/usr/bin/env python3
"""Cheap complete-shell diagnostics for p=5 boundary-orbit sources.

No finite solver is invoked.  The program reconstructs the two parity
vectors, their fixed lift masses, immediate parity-mass contradictions, and
zero-lift modular contradictions for every orbit.  It also counts distinct
parity-vector pairs so later exact solves can be deduplicated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p5_size_four_full_shell_mod5_batch as shell  # noqa: E402
from residual_boundary_four_lift_cpsat import atomic_write  # noqa: E402


def zero_lift_infeasible(eps: int, parity: np.ndarray) -> bool:
    dependencies = shell._SHELL_DATA[eps]["dependencies"]
    base = (
        dependencies[:, :2] @ np.asarray([21, 1], dtype=np.int64)
        + dependencies[:, 2:] @ (9 - parity.astype(np.int64))
    ) % 5
    return bool(np.any(base))


def scan(source_path: Path) -> dict:
    source = json.loads(source_path.read_text())
    if int(source["p"]) != 5:
        raise ValueError("source must have p=5")
    shell.build_linear_data()
    c_h = int(source["c_H"])
    mass_pairs = Counter()
    lift_pairs = Counter()
    status_counts = Counter()
    parity_pair_keys = set()
    rows = []
    for index, orbit in enumerate(source["orbits"]):
        boundary = tuple(int(value) for value in orbit["representative_vertices"])
        parities = {
            eps: shell.parity_vector(eps, c_h, boundary) for eps in (-1, 1)
        }
        masses = {eps: int(parities[eps].sum()) for eps in (-1, 1)}
        lifts = {
            eps: None if masses[eps] > 78 else (78 - masses[eps]) // 2
            for eps in (-1, 1)
        }
        immediate = []
        for eps in (-1, 1):
            if masses[eps] > 78:
                immediate.append([eps, "PARITY_MASS_INFEASIBLE"])
            elif lifts[eps] == 0 and zero_lift_infeasible(eps, parities[eps]):
                immediate.append([eps, "ZERO_LIFT_MOD5_INFEASIBLE"])
        status = "IMMEDIATE_EXCLUSION" if immediate else "NEEDS_BOUNDED_SOLVE"
        status_counts[status] += 1
        mass_pair = (masses[-1], masses[1])
        lift_pair = (lifts[-1], lifts[1])
        mass_pairs[mass_pair] += 1
        lift_pairs[lift_pair] += 1
        parity_key = hashlib.sha256(
            parities[-1].tobytes() + parities[1].tobytes()
        ).hexdigest()
        parity_pair_keys.add(parity_key)
        rows.append(
            {
                "orbit_index": index,
                "orbit_size": int(orbit["size"]),
                "boundary": list(boundary),
                "parity_masses": {str(eps): masses[eps] for eps in (-1, 1)},
                "lift_masses": {str(eps): lifts[eps] for eps in (-1, 1)},
                "status": status,
                "immediate_reasons": immediate,
                "parity_pair_sha256": parity_key,
            }
        )
    return {
        "experiment": "p5_full_shell_boundary_diagnostic",
        "status": "complete_no_solver_parity_and_zero_lift_scan",
        "source": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "p": 5,
        "c_H": c_h,
        "boundary_size": int(source.get("boundary_size", 4)),
        "infinity_value": int(source["infinity_value"]),
        "orbit_count": len(rows),
        "boundary_count": sum(int(row["orbit_size"]) for row in rows),
        "status_counts": dict(sorted(status_counts.items())),
        "distinct_parity_vector_pairs": len(parity_pair_keys),
        "parity_mass_pair_histogram": [
            {"masses": list(pair), "orbits": count}
            for pair, count in sorted(mass_pairs.items())
        ],
        "lift_mass_pair_histogram": [
            {"lifts": list(pair), "orbits": count}
            for pair, count in sorted(lift_pairs.items(), key=lambda item: str(item[0]))
        ],
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = scan(args.source)
    atomic_write(args.output, result)
    print(json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    main()
