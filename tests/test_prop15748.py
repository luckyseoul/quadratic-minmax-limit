import json
from pathlib import Path

from e1_gmin_m4_prop15748 import (
    EXPECTED_ALPHABET_SHA256,
    EXPECTED_PAYLOAD_SHA256,
    EXPECTED_Z2_SHA256,
    exact_literal_interpolation_certificate,
    literal_root_and_hard_alphabet_certificate,
    p5_excess_partition_reduction,
    proposition_15748,
    write_evidence,
)


ROOT = Path(__file__).resolve().parents[1]


def test_literal_roots_and_hard_alphabet_make_five_roots_impossible():
    row = literal_root_and_hard_alphabet_certificate()
    assert row["proved"] is True
    assert row["literal_star_even_degrees"] == [2, 4, 6]
    assert row["literal_star_power_sum_residue_sets"] == {
        "2": [0],
        "4": [0],
        "6": [0],
    }
    assert row["every_exact_literal_direction_is_a_common_M2_M4_M6_root"]
    assert row["hard_moment_alphabet_size"] == 69
    assert row["hard_moment_alphabet_sha256"] == EXPECTED_ALPHABET_SHA256
    assert row["hard_fourth_moment_value_set"] == list(range(1, 13))
    assert row["hard_fourth_moment_never_zero"] is True
    assert row["hard_fourth_values_when_N2_zero"] == [7, 8, 11]
    assert row["five_direction_subsets_checked"] == 2002
    assert row["five_root_quartic_evaluation_rank_set"] == [5]
    assert row["z_at_least_five_impossible"] is True


def test_z4_and_z3_are_empty_but_z2_has_exact_survivors_for_both_signs():
    row = exact_literal_interpolation_certificate()
    assert row["proved"] is True
    assert row["raw_payload_sha256"] == EXPECTED_PAYLOAD_SHA256
    assert row["z4_empty_for_both_hard_signs"] is True
    assert row["z3_empty_for_both_hard_signs"] is True
    assert row["z2_survivors_per_hard_sign"] == 336
    assert row["all_z2_z3_z4_cases_empty"] is False
    assert row["raw_runner_proved_flag_is_false_because_z2_survives"] is True

    for hard_sign in (-1, 1):
        sign = row["sign_rows"][str(hard_sign)]
        assert sign["proved"] is True
        assert sign["z4_parameter_cases_before_alphabet_filter"] == 420
        assert sign["z4_survivor_count"] == 0
        assert sign["z3_parameter_cases_before_alphabet_filter"] == 5880
        assert sign["z3_survivor_count"] == 0
        assert sign["z2_M2_M4_candidate_count_after_alphabet_filter"] == 1554
        assert sign["z2_N6_vectors_checked"] == 2688
        assert sign["z2_survivor_count"] == 336
        assert sign["z2_unique_survivor_count"] == 336
        assert sign["z2_survivor_catalog_sha256"] == EXPECTED_Z2_SHA256[hard_sign]
        assert sign["every_z2_survivor_independently_replayed"] is True


def test_only_the_five_ones_excess_partition_survives():
    row = p5_excess_partition_reduction()
    assert row["proved"] is True
    assert row["prior_prop_15747_forces_every_minimum_Q3_cell_literal"] is True
    assert row["opposite_excess_sum"] == 5
    assert row["prior_lower_bound_on_z"] == 2
    assert row["z_at_least_5_excluded_by_M4_root_count"] is True
    assert row["z4_excluded_by_exact_interpolation"] is True
    assert row["z3_excluded_by_exact_interpolation"] is True
    assert row["z2_moment_level_survivors_per_hard_sign"] == 336
    assert row["forced_z"] == 2
    assert row["positive_opposite_excess_count"] == 5
    assert row["only_remaining_opposite_excess_partition"] == [1] * 5
    assert row["moment_level_survivors_are_not_common_graph_realizations"] is True
    assert row["P5_branch_closed"] is False
    assert row["result_status"] == "proved open reduction"


def test_prop15748_evidence_records_reduction_not_closure(tmp_path):
    row = proposition_15748()
    assert row["proved"] is True
    assert row["result_status"] == (
        "exhaustive finite interpolation certificate and proved open reduction"
    )
    assert row["p13_t4_u4_P5_branch_closed"] is False
    assert row["p13_t4_u4_closed"] is False
    assert row["p13_k_eq_60_closed"] is False
    assert row["remaining_p13_t4_residues"] == [4, 6]
    assert row["residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved_means_exact_reduction_verified_not_all_cases_empty"] is True

    expected = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15748.json").read_text()
    )
    assert expected == row

    replay = tmp_path / "prop15748.json"
    assert write_evidence(replay) == replay
    assert json.loads(replay.read_text()) == row
    assert not list(tmp_path.glob("*.tmp"))
