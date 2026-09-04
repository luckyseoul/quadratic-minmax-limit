"""Exact local tests; no prime, graph, slice, or new Boolean catalog census."""
import json
from fractions import Fraction
from pathlib import Path

import pytest

from e1_gmin_m4_complement_triple_gap import (
    complement_triple_gap_certificate,
    local_bridge_package,
    p1_p_minus_one_local_exclusion,
    p1_p_plus_eleven_local_exclusion,
)


def test_neighboring_slice_bound_uses_composite_odd_order_without_prime_assumption() -> None:
    row = complement_triple_gap_certificate(29)
    section = row["neighboring_slice_bound"]
    assert section["outside_order"] == 26
    assert section["conditional_weights"] == [14, 12]
    assert section["odd_section_order"] == 25
    assert section["odd_section_need_not_be_prime"] is True
    assert section["identically_zero_coordinate_sections_at_most"] == 2
    assert section["nonzero_odd_section_mean_floor"] == "11/50"
    assert section["nonzero_neighbor_mean_floor"] == "66/325"
    assert section["one_nonzero_outer_contact_contribution_at_least"] == "132/25"
    gap = section["strict_gap_above_four"]
    assert gap["coefficients_in_p_descending"] == [1, -28, 99]
    assert gap["coefficients_in_x_descending"] == [1, 30, 128]
    assert gap["translation_checked_coefficientwise"] is True


def test_punctured_difference_is_not_treated_as_a_global_nonnegative_lift() -> None:
    row = complement_triple_gap_certificate(29)
    assert row["old_baseline_is_pointwise_parity_minimum_globally"] is False
    assert row["difference_L_is_nonnegative_on_contacts_only"] == [1, 2, 3]
    assert row["conditional_means_are_symmetrized_over_small_patterns_of_equal_weight"] is True
    assert row["quadrature_nodes"] == [1, 2, 3]
    assert row["quadrature_weights"] == ["39/58", "3/29", "13/58"]
    assert row["globalization"]["outside_slice_identity"] == "sum y=m-3+sum w"
    assert row["globalization"]["kernel_elimination_sign"] == "PLUS"
    assert row["globalization"]["contacts_in_w"] == [0, 2]


def test_excess_two_excluded_but_excess_four_has_exactly_three_offset_four_forms() -> None:
    row = complement_triple_gap_certificate(29)
    assert row["allowed_excesses_in_zero_to_four"] == [0, 4]
    assert row["excess_two_excluded"] is True
    assert row["excess_four_excluded"] is False
    assert row["excess_zero_coefficient_offset"] == 2
    assert row["excess_four_coefficient_offsets"] == [4]
    forms = row["excess_four_forms"]
    assert [form["pair"] for form in forms] == [[0, 1], [0, 2], [1, 2]]
    assert [form["complement_literal_coordinate"] for form in forms] == [2, 1, 0]
    assert all(form["half_difference_from_old_baseline"][0] == -1 for form in forms)
    assert forms[0]["values_in_binary_mask_order"] == [2, 1, 1, 2, 1, 0, 0, 1]
    assert all(form["coefficient_offset"] == 4 for form in forms)
    assert row["proved"] is True


def test_gap_theorem_itself_accepts_general_odd_not_just_prime_orders() -> None:
    row = complement_triple_gap_certificate(33)
    assert row["proved"] is True
    assert row["excess_two_excluded"] is True
    assert row["new_graph_prime_slice_or_catalog_census_used"] is False


def test_p_minus_one_height_equality_is_excluded_by_maximizing_half_mean_cubes() -> None:
    row = p1_p_minus_one_local_exclusion(29)
    height = row["height_at_least_two"]
    assert row["scaled_mass"] == 28
    assert height["paired_height_lower_bound"] == "8"
    assert height["stabilizer_height_upper_bound"] == "8"
    assert height["forced_height"] == 8
    assert height["average_maximizing_cube_mean"] == "1/2"
    assert height["therefore_every_maximizing_cube_mean"] == "1/2"
    assert height["half_mean_cube_maximum_upper_bound"] == 3
    assert height["excluded"] is True


def test_p_minus_one_boolean_gap_is_a_generic_positive_polynomial() -> None:
    row = p1_p_minus_one_local_exclusion(29)["height_one_boolean"]
    assert row["density"] == "7/29"
    assert row["largest_zero_influence_class_complement_bound"] == "310464/54665"
    assert row["junta_coordinates_at_most"] == 5
    assert row["generic_junta_gap_polynomial"]["coefficients_in_x_descending"] == [
        10, 830, 22912, 210312,
    ]
    assert row["target_strict_bracket"] == ["13/58", "15/58"]
    assert row["target_density_absent"] is True


def test_p_plus_eleven_height_endpoint_has_strict_half_mean_gap() -> None:
    row = p1_p_plus_eleven_local_exclusion(29)
    height = row["height_at_least_two"]
    assert row["scaled_mass"] == 40
    assert height["paired_height_lower_bound"] == "5"
    assert height["stabilizer_height_upper_bound"] == "80/7"
    assert height["average_maximizing_cube_mean_upper_bound"] == "5/7"
    assert height["strict_gap_below_three_quarters"] == "1/28"
    assert height["some_maximizing_cube_has_mean_exactly"] == "1/2"
    assert height["excluded"] is True


def test_p_plus_eleven_needs_the_correct_seven_coordinate_bound_at_p29() -> None:
    row = p1_p_plus_eleven_local_exclusion(29)["height_one_boolean"]
    assert row["density"] == "10/29"
    assert row["largest_zero_influence_class_complement_bound"] == "76608/10933"
    assert 7 < Fraction(row["largest_zero_influence_class_complement_bound"]) < 8
    assert row["junta_coordinates_at_most"] == 7
    assert row["cube_active_coordinates_at_most"] == 4
    assert row["generic_junta_gap_polynomial"]["coefficients_in_x_descending"] == [
        2, 190, 6776, 107936, 651360,
    ]
    assert row["target_strict_bracket"] == ["15/58", "14/29"]
    assert row["target_density_absent"] is True
    assert row["new_catalog_used"] is False


def test_parameterized_local_lemmas_replay_away_from_the_endpoint() -> None:
    for function in (p1_p_minus_one_local_exclusion, p1_p_plus_eleven_local_exclusion):
        row = function(37)
        assert row["proved"] is True
        assert row["excluded"] is True
        assert row["height_one_boolean"]["fixed_four_bit_catalog_reused"] is True


@pytest.mark.parametrize("parameter", [True, 27, 28, 29.0])
def test_gap_parameter_validation(parameter: object) -> None:
    with pytest.raises(ValueError):
        complement_triple_gap_certificate(parameter)


@pytest.mark.parametrize("parameter", [25, 31, 33])
def test_p1_local_parameter_validation(parameter: int) -> None:
    for function in (p1_p_minus_one_local_exclusion, p1_p_plus_eleven_local_exclusion):
        with pytest.raises(ValueError):
            function(parameter)


def test_package_keeps_local_scope_and_does_not_claim_an_endpoint() -> None:
    row = local_bridge_package()
    assert row["proved"] is True
    assert row["all_parameter_inequalities_proved_by_exact_identities_not_endpoint_sampling"] is True
    assert row["complement_triple_excess_four_is_attained_not_excluded"] is True
    assert row["standalone_endpoint_or_global_closure_claimed"] is False


def test_saved_local_evidence_equals_the_live_payload() -> None:
    root = Path(__file__).resolve().parents[1]
    saved = json.loads((root / "evidence/e1_gmin_m4_complement_triple_gap.json").read_text())
    row = local_bridge_package()
    assert saved == row
    assert (root / row["proof_note"]).is_file()
