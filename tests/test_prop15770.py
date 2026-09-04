"""Tests for Proposition 15.770."""

import json
from fractions import Fraction
from pathlib import Path

import pytest

from e1_gmin_m4_prop15734 import BRANCH_B2, BRANCH_P1_LAST, BRANCH_P3_LAST
from e1_gmin_m4_prop15770 import (
    P1_CARRIED_BRANCH,
    P1_NEW_SHARP_BRANCH,
    P3_CARRIED_BRANCH,
    P3_NEW_LOCAL_BRANCH,
    p1_next_layer_exclusion,
    p1_next_residue_ledger,
    p1_p_plus_thirteen_local_exclusion,
    p1_sharp_family_catalog,
    p3_next_layer_exclusion,
    p3_next_residue_ledger,
    p_minus_one_local_exclusion,
    proposition_15770,
    sharp_p_minus_three_classification_all_odd,
)


def test_p29_p_plus_thirteen_height_endpoint_uses_three_quarter_cube() -> None:
    row = p1_p_plus_thirteen_local_exclusion(29)
    height = row["height_at_least_two"]
    assert height["initial_raw_height_lower_bound"] == "9/2"
    assert height["initial_integral_height_lower_bound"] == 5
    assert height["refined_three_quarter_raw_height_lower_bound"] == "12"
    assert height["stabilizer_height_upper_bound"] == "12"
    assert height["forced_height"] == 12
    assert height["paired_cube_average_at_forced_height"] == "3/4"
    assert height["three_quarter_cube_maximum_upper_bound"] == 6
    assert row["p29_uses_three_quarter_cube_endpoint"] is True
    assert row["proved"] is True


def test_p29_p_plus_thirteen_boolean_case_misses_fixed_catalog() -> None:
    row = p1_p_plus_thirteen_local_exclusion(29)["height_one_boolean"]
    assert row["density"] == "21/58"
    assert Fraction(row["largest_zero_influence_class_complement_bound"]) == Fraction(
        391608, 54665
    )
    assert Fraction(row["largest_zero_influence_class_complement_bound"]) < 8
    assert row["junta_coordinates_at_most"] == 7
    assert row["cube_active_coordinates_at_most"] == 4
    assert row["target_density_absent"] is True


def test_sharp_p_minus_three_classification_extends_to_p1() -> None:
    row = sharp_p_minus_three_classification_all_odd(29)
    assert row["density"] == "13/58"
    assert row["matching_four_bit_layer_profiles"] == [
        [0, 0, 1, 2, 1],
        [1, 1, 0, 1, 1],
    ]
    assert row["selected_pair_table_count"] == 6
    assert row["all_equal_triple_table_count"] == 4
    assert row["original_slice_families"] == ["omitted_pair", "all_equal_triple"]
    assert row["proved"] is True


def test_p1_new_sharp_family_offsets_are_three_and_five() -> None:
    row = p1_sharp_family_catalog(29)
    assert [family["coefficient_offset"] for family in row["families"]] == [3, 5]
    assert row["difference_scaled_mass"] == 26
    assert row["proved"] is True


def test_p1_next_residue_ledger_has_exactly_four_arithmetic_survivors() -> None:
    row = p1_next_residue_ledger(29)
    assert row["t"] == 12
    assert row["k"] == 140
    assert row["H_edge_count"] == 141
    assert row["arithmetic_surviving_residues"] == [0, 11, 12, 14]
    live = {
        entry["u"]: [
            (candidate["b"], candidate["classification"])
            for candidate in entry["live_rows"]
        ]
        for entry in row["rows"]
        if entry["live_rows"]
    }
    assert live == {
        0: [(28, "exact")],
        11: [(26, "exact")],
        12: [(2, "sharp_p_minus_3")],
        14: [(2, "exact")],
    }


def test_p1_next_layer_closes_all_old_carried_and_new_branches() -> None:
    row = p1_next_layer_exclusion(29)
    branches = row["branch_exclusions"]
    assert set(branches) == {
        BRANCH_P1_LAST,
        BRANCH_B2,
        P1_CARRIED_BRANCH,
        P1_NEW_SHARP_BRANCH,
    }
    assert branches[P1_CARRIED_BRANCH]["hard_edge_count"] == 31
    assert branches[P1_CARRIED_BRANCH]["opposite_edge_count"] == 110
    assert branches[P1_CARRIED_BRANCH]["forced_high_direction_parallel_count"] == 3
    assert branches[P1_CARRIED_BRANCH]["forced_next_scaled_mean"] == 44
    sharp_rows = branches[P1_NEW_SHARP_BRANCH]["family_ledgers"]
    assert [entry["coefficient_offset"] for entry in sharp_rows] == [3, 5]
    assert all(entry["forced_next_scaled_mean"] == 42 for entry in sharp_rows)
    assert row["residual_ii_layer_excluded"] is True


def test_p1_complement_triple_excess_two_uses_punctured_gap_not_sharp_lift() -> None:
    row = p1_next_residue_ledger(29)
    residue = next(entry for entry in row["rows"] if entry["u"] == 12)
    candidate = next(entry for entry in residue["candidate_rows"] if entry["b"] == 26)
    assert candidate["excess"] == 2
    assert candidate["classification"] == "excluded_complement_triple_punctured_gap_two"
    dependency = row["complement_triple_punctured_gap_dependency"]
    assert dependency["old_baseline_is_pointwise_parity_minimum_globally"] is False
    assert dependency["excess_two_excluded"] is True
    assert dependency["excess_four_excluded"] is False


def test_p1_ledger_refuses_to_close_if_punctured_gap_dependency_fails(monkeypatch) -> None:
    import e1_gmin_m4_prop15770 as module

    monkeypatch.setattr(module, "complement_triple_gap_certificate", lambda p: {
        "proved": False, "excess_two_excluded": False,
    })
    with pytest.raises(ArithmeticError, match="punctured gap-two bridge"):
        module.p1_next_residue_ledger(29)


def test_p_minus_one_local_branch_is_boolean_and_catalog_absent() -> None:
    row = p_minus_one_local_exclusion(31)
    assert row["scaled_mass"] == 30
    assert row["H_at_least_two_scaled_floor"] == 32
    assert row["therefore_height_one_boolean"] is True
    assert row["density"] == "15/62"
    assert Fraction(row["largest_zero_influence_class_complement_bound"]) < 6
    assert row["junta_coordinates_at_most"] == 5
    assert row["target_density_absent"] is True
    assert row["proved"] is True


def test_p3_next_residue_ledger_separates_new_local_branch() -> None:
    row = p3_next_residue_ledger(31)
    assert row["t"] == 14
    assert row["k"] == 152
    assert row["H_edge_count"] == 153
    assert row["arithmetic_surviving_residues"] == [13, 14, 15]
    live = {
        entry["u"]: [candidate["classification"] for candidate in entry["live_rows"]]
        for entry in row["rows"]
        if entry["live_rows"]
    }
    assert live == {
        13: ["sharp_p_minus_3", "sharp_p_minus_3"],
        14: ["p_minus_one", "p_minus_one"],
        15: ["exact", "exact"],
    }


def test_p3_next_layer_carries_all_four_sharp_families() -> None:
    row = p3_next_layer_exclusion(31)
    branches = row["branch_exclusions"]
    assert set(branches) == {
        BRANCH_B2,
        BRANCH_P3_LAST,
        P3_CARRIED_BRANCH,
        P3_NEW_LOCAL_BRANCH,
    }
    carried = branches[P3_CARRIED_BRANCH]["family_ledgers"]
    assert [entry["coefficient_offset"] for entry in carried] == [2, 4, 3, 5]
    assert all(entry["forced_next_scaled_mean"] == 44 for entry in carried)
    assert all(
        entry["forced_high_direction_parallel_count"]
        == entry["coefficient_offset"] + 1
        for entry in carried
    )
    assert branches[P3_NEW_LOCAL_BRANCH]["difference_scaled_mass"] == 30
    assert row["residual_ii_layer_excluded"] is True


def test_parameter_checks_reject_wrong_ranges_and_congruences() -> None:
    with pytest.raises(ValueError):
        p1_next_layer_exclusion(17)
    with pytest.raises(ValueError):
        p1_next_layer_exclusion(31)
    with pytest.raises(ValueError):
        p3_next_layer_exclusion(23)
    with pytest.raises(ValueError):
        p3_next_layer_exclusion(37)


def test_packaged_scope_is_honest() -> None:
    row = proposition_15770()
    assert row["proved"] is True
    assert row["p23_same_layer_closed"] is True
    assert row["p23_original_k"] == 112
    assert row["p23_exceptional_eleven_root_certificate"]["proved"] is True
    assert row["p23_exceptional_eleven_root_certificate"]["p23_k112_closed"] is True
    assert row["finite_prime_graph_or_slice_census_used"] is False
    assert row["later_layers_closed"] is False
    assert row["residual_ii_closed_globally"] is False
    assert row["quadratic_minmax_limit_closed"] is False


def test_saved_evidence_equals_the_repaired_live_payload() -> None:
    root = Path(__file__).resolve().parents[1]
    saved = json.loads((root / "evidence/e1_gmin_m4_prop15770.json").read_text())
    row = proposition_15770()
    assert saved == json.loads(json.dumps(row))
    assert "punctured-gap theorem" in row["complement_triple_gap_two_justification"]
