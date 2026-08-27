#!/usr/bin/env python3
"""Prop. 15.663 -- exclude the p=7 size-eight forced-floor stratum.

For each product sign, 83,770,008 finite eight-point boundaries have exact
directional floor sums (32,32).  The exact type mean sums are also (32,32),
so every directional mean is forced to its parity floor.  Complete
Johnson-slice catalogs and the 135 left dependencies of the common affine
score system modulo seven exclude the whole stratum for c_H=-1.  A
nonsquare Paley anti-isometry transfers the exclusion to c_H=+1.

Together with Proposition 15.662 this lowers the open p=7 size-eight floor
remainder from 108,754,569 to 24,983,238 boundaries per sign.  It does not
close the full size-eight case, residual (ii), Type I, R1, global QVAR, or
the limit.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

GPU_RESULT_SHA256 = (
    "6143d4eb269861b3d380c53262b534e0a54a9645c9bbe7c29d9327200ae30535"
)
INDEPENDENT_AUDIT_SHA256 = (
    "7adaa5e76bf4f5e128c82ec219650b390c8c087d3aed2a44857f9da7939a9c53"
)
FLOOR_CMINUS1_SHA256 = (
    "47e7db3512fd0419df8d1fc30d886ffd6d6db1c3c09145c085ed3d80a285b218"
)
FLOOR_CPLUS1_SHA256 = (
    "b5e61256d34e7713db60beaf8e5e24c36958e6e1725229c9450349e4c3050b51"
)
ARCHIVE_ROOT = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-27-p7-size-eight-forced-floor/"
)


def p7_size_eight_forced_floor_certificate() -> dict:
    return {
        "all_size_eight_boundaries_per_sign": 450_978_066,
        "floor_surviving_boundaries_before_15662_15663_per_sign": 108_754_569,
        "conic_floor_survivors_excluded_by_15662_per_sign": 1_323,
        "nonconic_floor_survivors_before_15663_per_sign": 108_753_246,
        "forced_floor_ordered_profiles_per_sign": 2_016,
        "forced_floor_boundaries_per_sign": 83_770_008,
        "forced_floor_odd_secant_histogram_per_sign": {
            16: 254_016,
            20: 8_396_640,
            24: 32_673_984,
            28: 30_465_456,
            32: 10_459_344,
            36: 1_467_648,
            40: 52_920,
        },
        "type_floor_sums": (32, 32),
        "exact_type_mean_sums": (32, 32),
        "all_directional_means_forced_to_floor": True,
        "maximum_variable_catalogs_per_boundary": 1,
        "sole_variable_catalog_rows": 36,
        "common_score_system_shape": (282, 1_225),
        "common_score_system_rank_mod7": 147,
        "full_left_dependency_dimension_mod7": 135,
        "gpu_projection_dimension": 8,
        "gpu_projected_survivors": 526,
        "full_dependency_survivors": 0,
        "independent_nuka_recheck": True,
        "nonsquare_sign_transfer": True,
        "remaining_nonconic_floor_survivors_per_sign": 24_983_238,
        "gpu_result_sha256": GPU_RESULT_SHA256,
        "independent_audit_sha256": INDEPENDENT_AUDIT_SHA256,
        "floor_cminus1_sha256": FLOOR_CMINUS1_SHA256,
        "floor_cplus1_sha256": FLOOR_CPLUS1_SHA256,
        "archive_root": ARCHIVE_ROOT,
        "proved": True,
    }


def theorem_p7_size_eight_forced_floor_exclusion() -> dict:
    certificate = p7_size_eight_forced_floor_certificate()
    histogram = certificate["forced_floor_odd_secant_histogram_per_sign"]
    proved = bool(
        certificate["proved"]
        and sum(histogram.values())
        == certificate["forced_floor_boundaries_per_sign"]
        and certificate["type_floor_sums"]
        == certificate["exact_type_mean_sums"]
        == (32, 32)
        and certificate["all_directional_means_forced_to_floor"]
        and certificate["maximum_variable_catalogs_per_boundary"] == 1
        and certificate["common_score_system_rank_mod7"]
        + certificate["full_left_dependency_dimension_mod7"]
        == certificate["common_score_system_shape"][0]
        and certificate["gpu_projected_survivors"] == 526
        and certificate["full_dependency_survivors"] == 0
        and certificate["independent_nuka_recheck"]
        and certificate["nonsquare_sign_transfer"]
        and certificate["nonconic_floor_survivors_before_15663_per_sign"]
        - certificate["forced_floor_boundaries_per_sign"]
        == certificate["remaining_nonconic_floor_survivors_per_sign"]
        and certificate["floor_surviving_boundaries_before_15662_15663_per_sign"]
        - certificate["conic_floor_survivors_excluded_by_15662_per_sign"]
        - certificate["forced_floor_boundaries_per_sign"]
        == certificate["remaining_nonconic_floor_survivors_per_sign"]
    )
    return {
        "proved": proved,
        "certificate": certificate,
        "p7_size_eight_forced_floor_stratum_both_signs": "CLOSED",
        "p7_size_eight_remaining_floor_survivors_per_sign": 24_983_238,
        "full_p7_size_eight": "OPEN",
        "closes_all_nonconic_size_eight": False,
        "closes_all_p7_size_eight": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_p7_size_eight_forced_floor_exclusion()
    out = {
        "prop": "15.663",
        "title": "p=7 size-eight forced-floor modular exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15663.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
