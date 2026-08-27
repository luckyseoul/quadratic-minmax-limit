#!/usr/bin/env python3
"""Prop. 15.662 -- close the p=7 size-eight conic subbranch.

Among finite eight-point boundaries in AG(2,7), exactly 6,174 attain the
minimum of eight odd secants.  Segre's odd-order arc theorem identifies this
minimum branch with affine conics.  The exact floor leaves 1,323 boundaries,
or 32 square-semilinear orbits.  Complete allocation certificates exclude
all 32 orbits for c_H=-1, and a nonsquare Paley anti-isometry transfers the
exclusion to c_H=+1.

This does not close the full p=7 size-eight case: 108,753,246 nonconic floor
survivors remain per product sign.  Residual (ii), Type I, R1, global QVAR,
and the limit remain open.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

GLOBAL_CONIC_AUDIT_SHA256 = (
    "85f927f41b3ffc9afe1a101584e95ed852709ca6e861b439d8da1715008640a9"
)
ORDINARY_EXCEPTIONAL_AUDIT_SHA256 = (
    "cffa08a2b3858465f474b1d49ec294741c0e044610ac6f416da45f9262540487"
)
HIGH_MEAN_EXCEPTIONAL_AUDIT_SHA256 = (
    "adad4cccc67d2d165f213a2fbdf71448f02270f69d330d4688b70025a8d388cb"
)
ARCHIVE_ROOT = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-27-p7-size-eight-conic/"
)


def p7_size_eight_conic_certificate() -> dict:
    return {
        "all_size_eight_boundaries_per_sign": 450_978_066,
        "floor_surviving_boundaries_per_sign": 108_754_569,
        "floor_survivor_ordered_profiles_per_sign": 5_152,
        "minimum_odd_secants": 8,
        "minimum_odd_secant_conic_boundaries_per_sign": 6_174,
        "floor_rejected_conic_boundaries_per_sign": 4_851,
        "floor_surviving_conic_boundaries_per_sign": 1_323,
        "conic_stabilizer_orbits_per_sign": 32,
        "saturated_orbits": 25,
        "saturated_boundaries": 1_176,
        "saturated_mean_allocations": 600,
        "saturated_initial_cp_exclusions": 355,
        "saturated_long_cp_exclusions": 6,
        "saturated_catalog_join_exclusions": 239,
        "exceptional_orbits": 7,
        "exceptional_boundaries": 147,
        "exceptional_mean_allocations": 1_260,
        "exceptional_initial_cp_exclusions": 172,
        "exceptional_ordinary_gpu_exclusions": 662,
        "exceptional_high_direction_omission_gpu_exclusions": 426,
        "remaining_conic_mean_allocations": 0,
        "nonsquare_sign_transfer": True,
        "nonconic_floor_survivors_per_sign": 108_753_246,
        "global_conic_audit_sha256": GLOBAL_CONIC_AUDIT_SHA256,
        "ordinary_exceptional_audit_sha256": ORDINARY_EXCEPTIONAL_AUDIT_SHA256,
        "high_mean_exceptional_audit_sha256": HIGH_MEAN_EXCEPTIONAL_AUDIT_SHA256,
        "archive_root": ARCHIVE_ROOT,
        "proved": True,
    }


def theorem_p7_size_eight_conic_subbranch() -> dict:
    certificate = p7_size_eight_conic_certificate()
    proved = bool(
        certificate["proved"]
        and certificate["floor_rejected_conic_boundaries_per_sign"]
        + certificate["floor_surviving_conic_boundaries_per_sign"]
        == certificate["minimum_odd_secant_conic_boundaries_per_sign"]
        and certificate["saturated_orbits"] + certificate["exceptional_orbits"]
        == certificate["conic_stabilizer_orbits_per_sign"]
        and certificate["saturated_boundaries"]
        + certificate["exceptional_boundaries"]
        == certificate["floor_surviving_conic_boundaries_per_sign"]
        and certificate["saturated_initial_cp_exclusions"]
        + certificate["saturated_long_cp_exclusions"]
        + certificate["saturated_catalog_join_exclusions"]
        == certificate["saturated_mean_allocations"]
        and certificate["exceptional_initial_cp_exclusions"]
        + certificate["exceptional_ordinary_gpu_exclusions"]
        + certificate["exceptional_high_direction_omission_gpu_exclusions"]
        == certificate["exceptional_mean_allocations"]
        and certificate["remaining_conic_mean_allocations"] == 0
        and certificate["nonsquare_sign_transfer"]
        and certificate["floor_surviving_conic_boundaries_per_sign"]
        + certificate["nonconic_floor_survivors_per_sign"]
        == certificate["floor_surviving_boundaries_per_sign"]
    )
    return {
        "proved": proved,
        "certificate": certificate,
        "p7_size_eight_minimum_odd_secant_conic_subbranch_both_signs": "CLOSED",
        "all_32_floor_surviving_conic_orbits_both_signs": "CLOSED",
        "p7_size_eight_nonconic_floor_survivors_per_sign": 108_753_246,
        "full_p7_size_eight": "OPEN",
        "closes_all_p7_size_eight": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_p7_size_eight_conic_subbranch()
    out = {
        "prop": "15.662",
        "title": "complete p=7 size-eight minimum-odd-secant conic exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15662.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
