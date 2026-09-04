"""Fail-when-wrong tests for Proposition 15.769."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
from e1_gmin_m4_prop15734 import BRANCH_B2, BRANCH_P3_LAST
from e1_gmin_m4_prop15769 import (
    NEW_BRANCH,
    first_uncovered_p3_layer_exclusion,
    first_uncovered_p3_residue_ledger,
    hard_family_catalog,
    p23_local_threshold_witness,
    p_plus_thirteen_local_exclusion,
    proposition_15769,
    sharp_lift_branch_exclusion,
    sharp_p_minus_three_boolean_classification,
    sharp_p_minus_three_four_bit_catalog,
    write_evidence,
)


ROOT = Path(__file__).resolve().parents[1]


def test_sharp_p_minus_three_catalog_is_the_fixed_ten_table_catalog() -> None:
    row = sharp_p_minus_three_four_bit_catalog()
    assert row["source_catalog_valid_tables"] == 222
    assert row["selected_table_count"] == 10
    assert row["selected_pair_table_count"] == 6
    assert row["all_equal_triple_table_count"] == 4
    assert row["target_layer_profiles"] == [
        [0, 0, 1, 2, 1],
        [1, 1, 0, 1, 1],
    ]
    assert row["selected_tables_sha256"] == (
        "a52ce893775839abe2c2d3d5f7371f858591640358172ee77f0ff75a21e0dd34"
    )
    assert row["proved"] is True


def test_sharp_boolean_reduction_is_uniform_from_p31() -> None:
    for p in (31, 43, 47, 59):
        row = sharp_p_minus_three_boolean_classification(p)
        assert Fraction(row["largest_zero_influence_class_complement_bound"]) < 6
        assert row["junta_coordinates_at_most"] == 5
        assert row["five_less_than_both_slice_sides"] is True
        assert row["cube_active_coordinates_at_most"] == 4
        assert sorted(row["matching_four_bit_layer_profiles"]) == [
            [0, 0, 1, 2, 1],
            [1, 1, 0, 1, 1],
        ]
        assert row["original_slice_families"] == [
            "omitted_pair",
            "all_equal_triple",
        ]
        assert row["proved"] is True


def test_two_baselines_times_two_lifts_give_exact_offsets() -> None:
    row = hard_family_catalog(31)
    families = row["families"]
    assert sorted(family["coefficient_offset"] for family in families) == [
        2,
        3,
        4,
        5,
    ]
    assert {family["baseline"] for family in families} == {
        "complement_literal",
        "XNOR",
    }
    assert {family["lift"] for family in families} == {
        "omitted_pair",
        "all_equal_triple",
    }
    assert all(family["total_scaled_mean"] == 58 for family in families)
    assert row["proved"] is True


def test_p_plus_thirteen_local_height_and_boolean_cases_are_excluded() -> None:
    for p in (31, 43, 47, 59):
        row = p_plus_thirteen_local_exclusion(p)
        height = row["height_at_least_two"]
        boolean = row["height_one_boolean"]
        assert Fraction(height["height_lower_bound"]) > 3
        assert Fraction(height["paired_cube_average_upper_bound"]) < Fraction(3, 4)
        assert Fraction(boolean["largest_zero_influence_class_complement_bound"]) < 8
        assert boolean["junta_coordinates_at_most"] == 7
        assert boolean["seven_less_than_both_complementary_slice_sizes"] is True
        assert boolean["cube_coordinates_actually_needed_at_most"] == 4
        assert boolean["target_absent"] is True
        assert row["finite_prime_or_slice_census_used"] is False
        assert row["excluded"] is True


def test_p23_local_witness_only_blocks_the_parameterized_shortcut() -> None:
    row = p23_local_threshold_witness()
    assert row["layer_values"] == [3, 1, 0, 0, 1]
    assert row["mean"] == "9/23"
    assert row["scaled_mass_4p_E_C"] == 36
    assert row["scaled_mass_4p_E_C"] == row["p"] + 13
    assert row["is_only_a_local_quadratic_not_a_residual_graph"] is True


def test_first_uncovered_residue_ledger_has_only_two_residues() -> None:
    row = first_uncovered_p3_residue_ledger(31)
    assert row["q"] == 15
    assert row["m"] == 16
    assert row["layer_index_t"] == 13
    assert row["original_k"] == 150
    assert row["H_edge_count"] == 151
    assert row["surviving_residues"] == [13, 15]
    assert row["possible_branches"] == [BRANCH_B2, BRANCH_P3_LAST, NEW_BRANCH]
    assert row["new_branch_all_quotients_equal_one"] is True
    assert row["new_branch_scaled_mean"] == 58
    assert row["new_branch_sharp_lift_mass"] == 28
    assert row["proved"] is True


def test_new_branch_common_row_forces_four_separate_ledgers() -> None:
    row = sharp_lift_branch_exclusion(31)
    ledgers = sorted(
        row["family_ledgers"], key=lambda item: item["coefficient_offset"]
    )
    assert [ledger["forced_hard_parallel_count"] for ledger in ledgers] == [
        2,
        3,
        4,
        5,
    ]
    assert [ledger["hard_sign_times_global_T"] for ledger in ledgers] == [
        -87,
        -55,
        -23,
        9,
    ]
    assert [ledger["forbidden_mass_twelve_Q"] for ledger in ledgers] == [
        6,
        5,
        4,
        3,
    ]
    assert [ledger["forced_low_Q"] for ledger in ledgers] == [7, 6, 5, 4]
    assert all(
        ledger["surplus_after_every_Q_at_least_forced_low_Q"] == 7
        for ledger in ledgers
    )
    assert all(ledger["directions_at_forced_low_Q_at_least"] == 9 for ledger in ledgers)
    assert all(ledger["forced_low_Q_scaled_mean"] == 44 for ledger in ledgers)
    assert all(ledger["forced_low_Q_is_b_zero"] is True for ledger in ledgers)
    assert row["excluded"] is True


def test_old_endpoint_branches_reuse_p_plus_nine_and_new_branch_closes() -> None:
    row = first_uncovered_p3_layer_exclusion(31)
    branches = row["branch_exclusions"]
    for branch in (BRANCH_B2, BRANCH_P3_LAST):
        assert branches[branch]["forced_next_scaled_mean"] == 40
        assert branches[branch]["surplus_after_raising_every_Q"] == 9
        assert branches[branch]["local_dependency"] == (
            "Proposition 15.752 p+9 exclusion"
        )
        assert branches[branch]["excluded"] is True
    assert branches[NEW_BRANCH]["common_forced_opposite_local_mass"] == "p+13"
    assert branches[NEW_BRANCH]["excluded"] is True
    assert row["original_k"] == 150
    assert row["all_boundary_sizes_excluded"] is True
    assert row["residual_ii_layer_excluded"] is True


@pytest.mark.parametrize("bad_p", [23, 29, 33, 39, 41, 51])
def test_api_rejects_out_of_scope_parameters(bad_p: int) -> None:
    with pytest.raises(
        ValueError, match="prime p>=31 congruent to 3 modulo 4"
    ):
        first_uncovered_p3_layer_exclusion(bad_p)


def test_packaged_theorem_closes_one_family_not_residual_ii(tmp_path: Path) -> None:
    row = proposition_15769()
    assert row["prop"] == "15.769"
    assert row["first_layer_beyond_prop_15752"] is True
    assert row["p23_same_layer_closed"] is True
    assert row["p23_original_k"] == 110
    assert row["p23_exceptional_equality_moment_certificate"]["proved"] is True
    assert (
        row["p23_exceptional_equality_moment_certificate"]["p23_k110_closed"]
        is True
    )
    assert row["later_layers_closed"] is False
    assert "dual-bad/two-level" in row["prop_15274_slope_scope"]
    assert row["residual_ii_closed_globally"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert residual_ii_k_ge_4p_ND_closed() is False

    expected = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15769.json").read_text()
    )
    assert expected == row
    replay = tmp_path / "prop15769.json"
    assert write_evidence(replay) == replay
    assert json.loads(replay.read_text()) == row
    assert not list(tmp_path.glob("*.tmp"))
