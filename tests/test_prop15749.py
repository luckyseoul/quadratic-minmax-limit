"""Focused tests for Proposition 15.749."""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
from e1_gmin_m4_prop15749 import (
    EXPECTED_EVALUATION_ALPHABET_SHA256,
    EXPECTED_INTERSECTION,
    EXPECTED_INTERSECTION_SHA256,
    EXPECTED_Q4_MOMENT_SHA256,
    EXPECTED_Q4_ROW_SHA256,
    p13_t4_u4_close,
    proposition_15749,
    q4_translated_cut_moment_certificate,
    survivor_q4_moment_intersection_certificate,
    translated_cut_coordinate_bound_certificate,
    write_evidence,
)


ROOT = Path(__file__).resolve().parents[1]


def test_translated_cut_duals_force_integral_box() -> None:
    row = translated_cut_coordinate_bound_certificate()
    assert row["rational_lower_bound"] == "-52/9"
    assert row["rational_upper_bound"] == "26/15"
    assert row["integral_coordinate_bounds"] == [-5, 1]
    assert row["cut_catalog_invariant_under_all_six_permutations"] is True
    assert row["distance_action_transitive"] is True
    assert row["proved"] is True


def test_q4_moment_list_is_exact_and_hash_pinned() -> None:
    row = q4_translated_cut_moment_certificate()
    assert row["candidate_box_size"] == 117649
    assert row["admissible_row_count"] == 522
    assert row["admissible_row_sha256"] == EXPECTED_Q4_ROW_SHA256
    assert row["admissible_moment_count"] == 492
    assert row["admissible_moment_sha256"] == EXPECTED_Q4_MOMENT_SHA256
    assert row["proved"] is True


def test_survivor_intersection_forces_fourth_moment_zero() -> None:
    row = survivor_q4_moment_intersection_certificate()
    assert row["common_intersection"] == [list(value) for value in EXPECTED_INTERSECTION]
    assert row["common_intersection_has_N4_zero"] is True
    for sign in ("-1", "1"):
        sign_row = row["sign_rows"][sign]
        assert sign_row["z2_survivor_count"] == 336
        assert sign_row["evaluation_alphabet_size"] == 48
        assert sign_row["evaluation_alphabet_sha256"] == EXPECTED_EVALUATION_ALPHABET_SHA256
        assert sign_row["admissible_intersection_size"] == 12
        assert sign_row["admissible_intersection_sha256"] == EXPECTED_INTERSECTION_SHA256
        assert sign_row["compatible_direction_count_histogram"] == {
            "0": 252,
            "1": 42,
            "2": 42,
        }
        assert sign_row["maximum_compatible_Q4_directions_in_one_survivor"] == 2
        assert sign_row["proved"] is True
    assert row["proved"] is True


def test_prop15749_closes_only_u4_and_matches_evidence(tmp_path: Path) -> None:
    close = p13_t4_u4_close()
    row = proposition_15749()
    assert close["P5_branch_closed"] is True
    assert close["p13_t4_u4_closed"] is True
    assert close["total_forced_projective_M4_roots"] == 7
    assert close["M4_homogeneous_degree"] == 4
    assert close["remaining_p13_t4_residues"] == [6]
    assert row["p13_t4_u4_closed"] is True
    assert row["p13_k_eq_60_closed"] is False
    assert row["remaining_p13_t4_residues"] == [6]
    assert row["residual_ii_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert residual_ii_k_ge_4p_ND_closed() is False
    assert row["proved"] is True

    expected = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15749.json").read_text()
    )
    assert expected == row

    replay = tmp_path / "prop15749.json"
    assert write_evidence(replay) == replay
    assert json.loads(replay.read_text()) == row
    assert not list(tmp_path.glob("*.tmp"))
