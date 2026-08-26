#!/usr/bin/env python3
"""Prop. 15.658 -- close p=7 positive infinity-plus-five boundaries.

Directional scaled slacks of one quadratic type are congruent modulo
``p+1``.  At ``p=7``, ``c_H=+1``, and infinity in a six-point boundary,
all eight directions instead already saturate at scaled mean eight.  The
complete rank-21 ``J(7,4)`` classification gives a unique slack polynomial
for each possible odd-fibre size.  The resulting 282 edge equations have
135 left-null dependencies over F_7.  A complete sweep of all C(49,5)
finite boundaries leaves no compatible right side.  Independent V100 and
CPU implementations agree exactly.

This closes one p=7 size-six branch.  It does not close the negative-product
infinity branch, six-finite p=7, any p=5 size-six branch, larger boundaries,
residual (ii), R1, or the limit.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from e1_gmin_m4_prop15632 import scaled_direction_floor

ROOT = Path(__file__).resolve().parents[1]

ARCHIVE = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-26-p7-size-six-positive-infinity/"
)
V100_SHA256 = "7d3a828732a35192188a76958e5b9844c14668695cff455f3c61f8ffc5948b07"
NUKA_SHA256 = "05b622ffcb827e5b54411b27cfac91575f2efd44800775864f6fe76210875ffa"
SWEEP_SCRIPT_SHA256 = (
    "d0e54d2749a1fcd2841674134301fef241acbed2abdc5777aedfc5c36e87330a"
)


def same_type_scaled_slack_congruence(p: int) -> dict:
    """Verify the edgewise coefficient proof modulo p+1.

    For a direction of type ``eps``, an infinity edge contributes one.  A
    finite edge of type ``eps`` contributes ``p`` in its parallel direction
    and ``-1`` in every transverse direction; an opposite-type finite edge
    contributes ``+1`` throughout.  Thus every edge, and hence every scaled
    directional slack, has one common residue within a quadratic type.
    """
    if p < 3 or p % 2 == 0:
        raise ValueError("p must be odd")
    modulus = p + 1
    cases = {
        "infinity_edge": (1 % modulus, 1 % modulus),
        "same_type_finite_edge": (p % modulus, (-1) % modulus),
        "opposite_type_finite_edge": (1 % modulus, 1 % modulus),
    }
    return {
        "p": p,
        "modulus": modulus,
        "edge_coefficient_residues": cases,
        "proved": all(left == right for left, right in cases.values()),
    }


def p7_positive_infinity_slack_classification() -> dict:
    floors = {b: scaled_direction_floor(7, b, 0) for b in (1, 3, 5)}
    return {
        "phase": 0,
        "possible_odd_fibre_sizes": [1, 3, 5],
        "scaled_floors": floors,
        "type_budget": 32,
        "directions_per_type": 4,
        "all_scaled_means": 8,
        "unique_slacks": {
            "b=1,5": "A(X)=|X cap B| mod 2",
            "b=3": "A(X)=(|X cap B|-2)^2",
        },
        "proved": floors == {1: 8, 3: 8, 5: 8},
    }


def p7_positive_infinity_certificate() -> dict:
    return {
        "finite_boundaries": math.comb(49, 5),
        "equations": 282,
        "edge_variables": 1225,
        "rank_mod_7": 147,
        "left_dependency_dimension": 135,
        "v100_checked": 1_906_884,
        "v100_survivors": 0,
        "v100_elapsed_seconds": 2.8290963172912598,
        "nuka_cpu_checked": 1_906_884,
        "nuka_cpu_survivors": 0,
        "nuka_cpu_elapsed_seconds": 4.471272230148315,
        "matching_direction_mask_histogram": {
            "1": 2_923_536,
            "3": 9_507_960,
            "5": 2_823_576,
        },
        "v100_sha256": V100_SHA256,
        "nuka_sha256": NUKA_SHA256,
        "sweep_script_sha256": SWEEP_SCRIPT_SHA256,
        "archive": ARCHIVE,
        "independent_cpu_reproduction": True,
    }


def theorem_p7_positive_infinity_size_six_exclusion() -> dict:
    congruence = same_type_scaled_slack_congruence(7)
    slack = p7_positive_infinity_slack_classification()
    certificate = p7_positive_infinity_certificate()
    proved = bool(
        congruence["proved"]
        and slack["proved"]
        and certificate["finite_boundaries"] == 1_906_884
        and certificate["v100_checked"] == certificate["finite_boundaries"]
        and certificate["nuka_cpu_checked"] == certificate["finite_boundaries"]
        and certificate["v100_survivors"] == 0
        and certificate["nuka_cpu_survivors"] == 0
        and certificate["rank_mod_7"]
        + certificate["left_dependency_dimension"]
        == certificate["equations"]
        and certificate["independent_cpu_reproduction"]
    )
    return {
        "proved": proved,
        "same_type_slack_congruence_mod_p_plus_1": congruence,
        "p7_cH_positive_infinity_plus_five_finite": "CLOSED",
        "slack_classification": slack,
        "certificate": certificate,
        "p7_cH_negative_infinity_plus_five_finite": "OPEN",
        "p7_six_finite": "OPEN",
        "p5_size_six": "OPEN",
        "boundaries_size_at_least_eight": "OPEN",
        "closes_all_p7_size_six": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_p7_positive_infinity_size_six_exclusion()
    out = {
        "prop": "15.658",
        "title": "p=7 positive infinity-plus-five boundary exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15658.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
