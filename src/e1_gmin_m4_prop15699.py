#!/usr/bin/env python3
"""Prop. 15.699 -- close the p=19 second all-finite endpoint.

After Proposition 15.698, three parity profiles of slack 24, 28, and 32
remain.  Their exact affine-Radon inverse/profile models are UNSAT before
any edge-lift variables or floor assumptions are introduced.  Completed
runs on nuka, jellyfin, and soulkiller ECC close all three profiles.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15688 import p19_residue_zero_profiles


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "evidence" / "p19_endpoint_boundary_unsat"
RAW_HASHES = {
    "slack24_nuka.json": "154b94df1c192fa7770c715f5bbc4c0c99f529e13c7bc9faf54385eef8bea16e",
    "slack24_soulkiller_ecc.json": "0adc5d01c140cedf0c3c5bc01ad9bf2134945cd711c6b61281cbe39daafc784c",
    "slack28_soulkiller.json": "5dd0a19d449d74d158c625870a716f19dc2b235a5e6795176c950944c6d7b5d7",
    "slack32_jellyfin.json": "b633d52c6edbdabbd7cc765de0f5bd81487621a0f34796014e7640f619936dc5",
    "slack32_soulkiller_ecc.json": "b64d63404dcb4406d67077f8a5d02071b89b3b0e38ee32e6b541d1696a74d403",
}
EXPECTED = {
    0: (32, {"0": {"0": 5, "16": 5}, "1": {"2": 9, "14": 1}}),
    1: (28, {"0": {"0": 5, "16": 5}, "1": {"2": 9, "10": 1}}),
    2: (24, {"0": {"0": 5, "16": 5}, "1": {"2": 9, "6": 1}}),
}


def _atomic_write(path: Path, payload: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def p19_endpoint_boundary_unsat_certificate() -> dict[str, object]:
    profiles = [
        row for row in p19_residue_zero_profiles()["profiles"]
        if int(row["pair_slack"]) in (24, 28, 32)
        and row["phase_profiles_b"]["0"] == {0: 5, 16: 5}
    ]
    if Counter(int(row["pair_slack"]) for row in profiles) != Counter({24: 1, 28: 1, 32: 1}):
        raise ArithmeticError("final p=19 profile ledger changed")

    runs = []
    counts: Counter[int] = Counter()
    for filename, expected_hash in RAW_HASHES.items():
        path = ARCHIVE / filename
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected_hash:
            raise ArithmeticError(f"raw endpoint hash changed: {filename}")
        row = json.loads(path.read_text())
        index = int(row["profile_index"])
        slack, phase_profile = EXPECTED[index]
        if not (
            row["solver_status"] == "UNSATISFIABLE"
            and row["finite_infeasibility_only"] is True
            and row["feasible_boundary_profile"] is False
            and row["repair_five_constraints"] is False
            and row["pair_slack"] == slack
            and row["phase_profiles_b"] == phase_profile
            and row["native_xor_constraints"] == 741
            and row["cardinality_constraints"] == 45
            and row["normalization"]["mode"] == "phase-zero-b0-pair"
        ):
            raise ArithmeticError(f"raw endpoint UNSAT record changed: {filename}")
        counts[slack] += 1
        runs.append({
            "file": filename,
            "sha256": expected_hash,
            "pair_slack": slack,
            "threads": row["threads"],
            "clauses": row["clauses"],
            "solve_seconds": row["solve_seconds"],
            "solver_status": row["solver_status"],
        })
    if counts != Counter({24: 2, 28: 1, 32: 2}):
        raise ArithmeticError("endpoint replay coverage changed")
    return {
        "profiles": sorted(profiles, key=lambda row: int(row["pair_slack"])),
        "model": {
            "boundary_size": 16,
            "point_variables": 361,
            "radon_equations": "r=A*x and x=A^T*r over F_2",
            "uses_edge_lift_variables": False,
            "uses_floor_relaxation": False,
            "normalization": (
                "a pair on a phase-zero b=0 line is sent to 0,1 by a "
                "square affine similarity"
            ),
            "runs": runs,
            "all_completed_unsat": True,
        },
        "sign_transfer": {
            "nonsquare_dilation_flips_eps_and_c_H_together": True,
            "both_c_H_signs_excluded": True,
        },
        "excluded_slacks": [24, 28, 32],
        "proved_computationally": True,
    }


def p19_second_endpoint_exclusion() -> dict[str, object]:
    certificate = p19_endpoint_boundary_unsat_certificate()
    return {
        "proposition": "15.699",
        "certificate": certificate,
        "p19_profiles_before": 3,
        "p19_profiles_after": 0,
        "remaining_slack_histogram": {},
        "p19_second_all_finite_endpoint_closed": True,
        "closes_residual_ii": False,
        "closes_R1": False,
        "closes_type_I": False,
        "L_status": "OPEN",
        "proved": True,
    }


def main() -> None:
    theorem = p19_second_endpoint_exclusion()
    _atomic_write(ROOT / "evidence" / "e1_gmin_m4_prop15699.json", theorem)
    print("Prop. 15.699: p=19 second all-finite endpoint closed")


if __name__ == "__main__":
    main()
