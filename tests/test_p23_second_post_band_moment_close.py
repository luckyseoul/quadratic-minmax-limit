"""Fail-when-wrong tests for the exceptional p=23 second post-band close."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_p23_second_post_band_moment_close import (
    P3_CARRIED_BRANCH,
    P3_NEW_LOCAL_BRANCH,
    p23_carried_sharp_moment_exclusion,
    p23_p_minus_one_local_exclusion,
    p23_second_post_band_moment_close,
    p23_second_post_band_residue_ledger,
    write_evidence,
)
from e1_gmin_m4_prop15734 import BRANCH_B2, BRANCH_P3_LAST


ROOT = Path(__file__).resolve().parents[1]


def test_residue_ledger_leaves_exactly_u9_u10_u11() -> None:
    row = p23_second_post_band_residue_ledger()
    assert row["layer_index_t"] == 10
    assert row["original_k"] == 112
    assert row["H_edge_count"] == 113
    assert row["guaranteed_isolated_vertices"] == 304
    assert row["arithmetic_surviving_residues"] == [9, 10, 11]
    live = {
        entry["u"]: [candidate["classification"] for candidate in entry["live_rows"]]
        for entry in row["rows"]
        if entry["live_rows"]
    }
    assert live == {
        9: ["sharp_p_minus_3", "sharp_p_minus_3"],
        10: ["p_minus_one", "p_minus_one"],
        11: ["exact", "exact"],
    }
    assert row["proved"] is True


def test_mass_p_minus_one_branch_is_boolean_and_catalog_absent() -> None:
    row = p23_p_minus_one_local_exclusion()
    assert row["scaled_mass"] == 22
    assert row["H_at_least_two_scaled_floor"] == 24
    assert row["therefore_height_one_boolean"] is True
    assert row["density"] == "11/46"
    assert Fraction(row["largest_zero_influence_class_complement_bound"]) == Fraction(
        5929, 1058
    ) < 6
    assert row["junta_coordinates_at_most"] == 5
    assert row["cube_active_coordinates_at_most"] == 4
    assert row["target_density_absent"] is True
    assert row["proved"] is True


def test_carried_sharp_branch_has_eleven_roots_and_one_f5_survivor() -> None:
    row = p23_carried_sharp_moment_exclusion()
    ledgers = row["family_ledgers"]
    assert [entry["coefficient_offset"] for entry in ledgers] == [2, 4, 3, 5]
    assert all(entry["low_hard_direction_count"] == 11 for entry in ledgers)
    assert all(entry["unique_high_direction_count"] == 1 for entry in ledgers)
    assert all(
        entry["forced_high_parallel_count"] == entry["coefficient_offset"] + 1
        for entry in ledgers
    )
    assert all(entry["forbidden_scaled_mass"] == 12 for entry in ledgers)
    assert all(entry["forced_low_scaled_mass"] == 36 for entry in ledgers)
    assert all(
        entry["surplus_after_every_opposite_Q_at_least_forced_low"] == 4
        for entry in ledgers
    )
    assert all(entry["directions_at_forced_low_Q_at_least"] == 8 for entry in ledgers)
    assert row["unique_survivor_before_moments"] == {
        "hard_P": 4,
        "opposite_Q": 5,
        "opposite_form": "F5",
    }
    assert row["low_triangle_minus_star_projective_roots"] == 11
    assert row["maximum_common_form_degree"] == 8
    assert row["low_roots_force_G4_and_G8_identically_zero"] is True
    assert row["opposite_K5_simultaneous_zero_count"] == 0
    assert row["excluded"] is True
    assert row["proved"] is True


def test_old_endpoints_and_every_new_branch_are_excluded() -> None:
    row = p23_second_post_band_moment_close()
    branches = row["branch_exclusions"]
    assert set(branches) == {
        BRANCH_B2,
        BRANCH_P3_LAST,
        P3_CARRIED_BRANCH,
        P3_NEW_LOCAL_BRANCH,
    }
    assert branches[BRANCH_B2]["forced_next_scaled_mass"] == 32
    assert branches[BRANCH_P3_LAST]["forced_next_scaled_mass"] == 32
    assert branches[P3_CARRIED_BRANCH]["proved"] is True
    assert branches[P3_NEW_LOCAL_BRANCH]["difference_scaled_mass"] == 22
    assert all(branch["proved"] for branch in branches.values())


def test_package_closes_only_p23_k112_and_replays_evidence(tmp_path: Path) -> None:
    row = p23_second_post_band_moment_close()
    assert row["p23_k112_closed"] is True
    assert row["all_boundary_sizes_excluded"] is True
    assert row["new_graph_or_residual_configuration_census_used"] is False
    assert row["fixed_p23_five_set_coefficient_certificate_reused"] is True
    assert row["later_layers_closed"] is False
    assert row["residual_ii_closed_globally"] is False
    assert row["E1_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["proved"] is True

    expected = json.loads(
        (
            ROOT
            / "evidence"
            / "e1_gmin_m4_p23_second_post_band_moment_close.json"
        ).read_text()
    )
    assert expected == row
    replay = tmp_path / "p23-second.json"
    assert write_evidence(replay) == replay
    assert json.loads(replay.read_text()) == row
    assert not list(tmp_path.glob("*.tmp"))
