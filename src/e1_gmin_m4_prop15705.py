#!/usr/bin/env python3
"""Prop. 15.705 -- exclude the final p=17 slack-sixteen profiles.

Every slack-sixteen boundary can be repaired to a 12-arc by deleting four
points (padding a shallower repair if necessary).  For a repaired core A and
deleted point x, let mu_A(x) be the number of A-secants through x.  Core
secants contribute at least 4 sum_x mu_A(x) to pair slack, so only four-point
extensions with sum mu_A(x) <= 4 need be considered.

Orbiter build 3361 gives 629 PGL classes of 12-arcs in PG(2,17), independently
split here as Sticker's published 553 complete and 76 extendible classes.  The
exact extension census has 97,122 charge-admissible quadruples.  Only 47
extensions (ten distinct point sets) have one of the three required global
line patterns, always four 3-secants.  Enumerating all 6,345 disjoint choices
of line at infinity gives 317 unlabelled Paley-phase profiles.  Neither phase
labelling meets any of the thirteen arithmetic targets.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from e1_gmin_m4_prop15704 import p17_slack_sixteen_free_direction_exclusion


ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "evidence" / "p17_arcs_d2_reps_lvl_12.csv"
CENSUS = ROOT / "evidence" / "p17_slack16_orbiter_extension.json"
CSV_SHA256 = "0a57481731e10d55eb16a24158d57ca738240a9b32d3f66b9a39d85a64f16e24"


def p17_slack_sixteen_orbit_exclusion() -> dict[str, object]:
    """Proposition 15.705, conditional on Orbiter's exhaustive PGL census."""
    previous = p17_slack_sixteen_free_direction_exclusion()
    census = json.loads(CENSUS.read_text())
    if hashlib.sha256(CSV.read_bytes()).hexdigest() != CSV_SHA256:
        raise ArithmeticError("p17 Orbiter representative archive changed")
    expected = {
        "twelve_arc_orbits": 629,
        "complete_twelve_arc_orbits": 553,
        "extendible_twelve_arc_orbits": 76,
        "target_profile_count": 13,
        "raw_four_point_extensions_with_core_secant_charge_at_most_four": 97122,
        "occupancy_valid_extensions": 47,
        "unique_boundary_rank_sets": 10,
        "disjoint_line_affine_chart_cases": 6345,
        "distinct_unlabelled_phase_profiles": 317,
        "phase_labelled_target_hits": 0,
    }
    if any(int(census[key]) != value for key, value in expected.items()):
        raise ArithmeticError("p17 slack-sixteen orbit census changed")
    if census["occupancy_pattern_histogram"] != {"n3=4,n4=0": 47}:
        raise ArithmeticError("p17 slack-sixteen line patterns changed")
    if census["source_sha256"] != CSV_SHA256 or census["hits"]:
        raise ArithmeticError("p17 slack-sixteen source/hit ledger changed")

    before = int(previous["profile_count_after"])
    excluded = int(previous["remaining_slack_sixteen_profiles"])
    after = before - excluded
    histogram = dict(previous["remaining_pair_slack_histogram"])
    del histogram[16]
    if before != 654 or excluded != 13 or after != 641 or sum(histogram.values()) != after:
        raise ArithmeticError("p17 post-orbit accounting changed")
    return {
        "proposition": "15.705",
        "p": 17,
        "boundary_size": 16,
        "profile_count_before": before,
        "profiles_excluded_here": excluded,
        "profile_count_after": after,
        "remaining_pair_slack_histogram": histogram,
        "remaining_slack_zero_profiles": 2,
        "remaining_profiles_of_slack_at_least_twenty": after - 2,
        "remaining_slack_sixteen_profiles": 0,
        "repair_normal_form": (
            "delete four points to a 12-arc; core secants force "
            "sum of deleted-point secant indices at most four"
        ),
        "orbit_extension_census": census,
        "p17_second_all_finite_endpoint_closed": False,
        "top_level_gates_changed": False,
        "proved_conditional_on_orbiter_exhaustive_pgl_orbit_census": True,
    }


def main() -> None:
    theorem = p17_slack_sixteen_orbit_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15705.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.705: p=17 second-boundary profiles "
        f"{theorem['profile_count_before']} -> {theorem['profile_count_after']}"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
