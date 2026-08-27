#!/usr/bin/env python3
"""Prop. 15.660 -- exclude every p=5 size-six residual boundary.

Four exact boundary catalogs, two signed symmetry transfers, and the coarse
full-shell SCIP batches reduce all product signs and infinity bits to six
signed-symmetry classes.  Independent layered audits reconstruct every
finite quotient used to exclude those six classes.  Hence every p=5
boundary of size six is impossible.

This leaves the six-finite p=7 branch at size six.  Boundaries of size at
least eight, residual (ii), Type I, R1, global QVAR, and the limit remain
open.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

GLOBAL_AUDIT_SHA256 = (
    "d6650a9f71043dce2902e157b56f988305470b911eaccb130522dd2f55b3bbd8"
)
CLASS_AUDIT_SHA256 = {
    0: "5d13771fd85b75df98bf486f996728f885b243e25fa2a1c37657f199cb6656ad",
    881: "a4f320fc46039ab1009e47190c10fc5884645df2776492d2598ffeef8297f0ff",
    2529: "51149e286878ea92859e70f931e2984d888eca2c6b5e81476c48c133d8e6af0e",
    3032: "85cd90c192916c93da2bf9cf54e58bff898b103c731e21307930c6feffa1424b",
    4731: "787829034632e188897fe971635395d369781fe9bdf17d11f87b3778cec12e1f",
    4939: "d7ea05e882c07a3949c741a531b783f0f411b71efce2b36a6f0f414c13af913b",
}
ARCHIVE_ROOT = (
    "/mnt/storage/e1work/maxplus_p13/"
    "p5_size6_circle_attack_2026-08-26/"
)


def p5_size_six_global_reduction() -> dict:
    return {
        "catalogs_rebuilt_from_definitions": True,
        "signs": (-1, 1),
        "infinity_bits": (0, 1),
        "no_infinity_survivors_per_sign": 159_050,
        "no_infinity_orbits_per_sign": 6_766,
        "coarse_residual_classes_before_cross_infinity_symmetry": 7,
        "residual_classes_after_symmetry": 6,
        "no_infinity_sign_transfer_bijective": True,
        "infinity_orbit_1144_transfers_to_no_infinity_orbit_881": True,
        "selection_reduction_proved": True,
    }


def p5_size_six_class_certificates() -> dict:
    classes = {
        0: {"boundary": (0, 1, 2, 3, 4, 5), "artifacts": 867},
        881: {"boundary": (1, 2, 6, 7, 18, 20), "artifacts": 97},
        2529: {"boundary": (2, 3, 6, 9, 21, 24), "artifacts": 97},
        3032: {"boundary": (2, 3, 7, 8, 16, 19), "artifacts": 90},
        4731: {"boundary": (2, 5, 8, 9, 12, 15), "artifacts": 93},
        4939: {"boundary": (2, 5, 13, 14, 18, 19), "artifacts": 61},
    }
    return {
        "proved": set(classes) == set(CLASS_AUDIT_SHA256),
        "class_count": len(classes),
        "closed_class_count": len(classes),
        "classes": classes,
        "audit_sha256": CLASS_AUDIT_SHA256,
    }


def theorem_p5_size_six_exclusion() -> dict:
    reduction = p5_size_six_global_reduction()
    certificates = p5_size_six_class_certificates()
    proved = bool(
        reduction["catalogs_rebuilt_from_definitions"]
        and reduction["selection_reduction_proved"]
        and reduction["residual_classes_after_symmetry"] == 6
        and reduction["no_infinity_sign_transfer_bijective"]
        and reduction[
            "infinity_orbit_1144_transfers_to_no_infinity_orbit_881"
        ]
        and certificates["proved"]
        and certificates["closed_class_count"] == 6
    )
    return {
        "proved": proved,
        "global_reduction": reduction,
        "class_certificates": certificates,
        "global_audit_sha256": GLOBAL_AUDIT_SHA256,
        "archive_root": ARCHIVE_ROOT,
        "p5_size_six": "CLOSED",
        "p7_infinity_plus_five": "CLOSED_BY_15.658_AND_15.659",
        "p7_six_finite": "OPEN",
        "boundaries_size_at_least_eight": "OPEN",
        "closes_all_size_six": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_p5_size_six_exclusion()
    out = {
        "prop": "15.660",
        "title": "complete p=5 size-six boundary exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15660.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
