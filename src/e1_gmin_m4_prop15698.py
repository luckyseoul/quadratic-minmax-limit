#!/usr/bin/env python3
"""Prop. 15.698 -- exclude the p=19 all-b2 slack-twenty profile.

Proposition 15.694 forces every slack-twenty boundary to split into an
11-arc core and five deleted points, each on exactly one core secant.  The
exact native-XOR model imposes that repair structure together with the
complete affine Radon equations and the profile ``{0:5,16:5}/{2:10}``.
Two completed CryptoMiniSat runs return UNSAT, so the boundary itself does
not exist and no edge lift can realize the profile.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path

from e1_gmin_m4_prop15688 import p19_residue_zero_profiles
from e1_gmin_m4_prop15694 import p19_slack_twenty_equality_normal_form
from e1_gmin_m4_prop15697 import p19_allb2_infinity_degree_reduction


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "evidence" / "p19_allb2_boundary_unsat"
RAW_HASHES = {
    "nuka.json": "cb38e482d3610e214c23364711937fdce8e969caec9242ba80bb38cbcac5e4c9",
    "soulkiller_ecc.json": "a580391837db13d890a46bf7cd0f55c677fd9e4f5e653dd120d6b9349fb21935",
}


def _atomic_write(path: Path, payload: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def p19_allb2_boundary_unsat_certificate() -> dict[str, object]:
    profile = [
        row for row in p19_residue_zero_profiles()["profiles"]
        if int(row["pair_slack"]) == 20
        and row["phase_profiles_b"] == {"0": {0: 5, 16: 5}, "1": {2: 10}}
    ]
    if len(profile) != 1:
        raise ArithmeticError("all-b2 slack-twenty profile changed")
    repair = p19_slack_twenty_equality_normal_form()["repair_normal_form"]
    if repair != {
        "core_size": 11,
        "deleted_size": 5,
        "core_is_arc": True,
        "deleted_is_arc": True,
        "deleted_core_secant_multiplicities": [1] * 5,
        "global_slack_equality": "slack(S)=4*sum_{x in D} mu_A(x)=20",
    }:
        raise ArithmeticError("slack-twenty repair normal form changed")

    runs = []
    for filename, expected_hash in RAW_HASHES.items():
        path = ARCHIVE / filename
        if _sha256(path) != expected_hash:
            raise ArithmeticError(f"raw solver hash changed: {filename}")
        row = json.loads(path.read_text())
        if not (
            row["solver_status"] == "UNSATISFIABLE"
            and row["finite_infeasibility_only"] is True
            and row["feasible_boundary_profile"] is False
            and row["repair_five_constraints"] is True
            and row["profile_index"] == 3
            and row["pair_slack"] == 20
            and row["phase_profiles_b"] == {"0": {"0": 5, "16": 5}, "1": {"2": 10}}
            and row["native_xor_constraints"] == 741
            and row["clauses"] == 1184892
            and row["normalization"]["mode"] == "phase-zero-b0-pair"
        ):
            raise ArithmeticError(f"raw UNSAT record changed: {filename}")
        runs.append(
            {
                "file": filename,
                "sha256": expected_hash,
                "threads": row["threads"],
                "solve_seconds": row["solve_seconds"],
                "solver_status": row["solver_status"],
            }
        )

    degrees = p19_allb2_infinity_degree_reduction()["remaining_infinity_degrees"]
    if degrees != [0, 20, 38] or any(value & 1 for value in degrees):
        raise ArithmeticError("all-b2 infinity parity changed")
    return {
        "profile": profile[0],
        "repair_normal_form": repair,
        "model": {
            "point_variables": 361,
            "affine_line_parity_variables": 380,
            "radon_equations": "r=A*x and x=A^T*r over F_2",
            "normalization": (
                "choose a retained core point and its partner on any b=0 "
                "phase-zero line, then apply a square affine similarity to 0,1"
            ),
            "raw_runs": runs,
            "all_completed_unsat": True,
        },
        "sign_transfer": {
            "admissible_infinity_degrees": degrees,
            "finite_finite_edge_count_is_odd": True,
            "nonsquare_dilation_flips_c_H_and_direction_type": True,
            "both_c_H_signs_excluded": True,
        },
        "profile_excluded": True,
        "proved_computationally": True,
    }


def p19_allb2_profile_exclusion() -> dict[str, object]:
    certificate = p19_allb2_boundary_unsat_certificate()
    return {
        "proposition": "15.698",
        "certificate": certificate,
        "p19_profiles_before": 4,
        "p19_profiles_after": 3,
        "remaining_slack_histogram": {24: 1, 28: 1, 32: 1},
        "all_p19_slack_twenty_profiles_closed": True,
        "p19_second_all_finite_endpoint_closed": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "closes_type_I": False,
        "L_status": "OPEN",
        "proved": True,
    }


def main() -> None:
    theorem = p19_allb2_profile_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15698.json"
    _atomic_write(target, theorem)
    print("Prop. 15.698: p=19 all-b2 slack-twenty boundary is UNSAT; three profiles remain")


if __name__ == "__main__":
    main()
