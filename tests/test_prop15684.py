import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15684 import (
    conic_core_repair_lemma,
    line_pair_slack,
    low_value_cube_certificate,
    main as prop15684_main,
    p23_arc_envelope_exclusion,
    p23_complete_arc_classification,
    p23_endpoint_residue_ledger,
    p23_low_slack_conic_exclusion,
    p23_reduction_theorem,
    p23_residue_zero_profile_census,
    p23_small_mass_exclusion,
    p23_u9_open_profile,
)


ROOT = Path(__file__).resolve().parents[1]


def test_low_value_cube_certificate_has_the_needed_three_eighths_floor():
    row = low_value_cube_certificate()
    assert row["nonzero_degree_two_mean_floor"] == Fraction(1, 4)
    assert row["first_bit_degree"] == 2
    assert row["second_bit_degree"] == 4
    assert row["second_bit_support_floor"] == Fraction(1, 16)
    assert row["value_three_mean_floor"] == Fraction(3, 8)
    assert row["proved"] is True


def test_scaled_masses_twelve_and_sixteen_are_both_excluded():
    row = p23_small_mass_exclusion()
    assert row["stabilizer_maximum_inequality"] == "c>=4H"
    assert row["paired_cube_dimension"] == 11
    assert row["excluded_scaled_masses"] == [12, 16]
    cube_cases = {
        (item["scaled_mass"], item["maximum_height"])
        for item in row["case_rows"]
        if item["excluded_by_cube"]
    }
    assert cube_cases == {
        (12, 1),
        (12, 2),
        (12, 3),
        (16, 1),
        (16, 2),
        (16, 3),
    }
    shell = row["height_four_shell_case"]
    assert shell["domain_v2_dimension"] == 253
    assert shell["restriction_image_dimension"] == 230
    assert shell["restriction_kernel_dimension"] == 23
    assert shell["factor_subspace_dimension"] == 23
    assert "nonzero modulo 3" in shell["integrality_obstruction"]
    assert shell["proved"] is True


def test_exact_p23_residue_ledger_restores_open_u9():
    row = p23_endpoint_residue_ledger()
    assert [item["u0"] for item in row["pair_survivors"]] == [
        0,
        2,
        3,
        4,
        5,
        6,
        8,
        9,
    ]
    positive = {item["u0"]: item for item in row["positive_residue_rows"]}
    assert positive[6]["scaled_mass_c"] == 12
    assert positive[8]["scaled_mass_c"] == 16
    assert all(positive[u]["excluded"] for u in (2, 3, 4, 5, 6, 8))
    assert positive[9]["scaled_mass_c"] == 18
    assert positive[9]["excluded"] is False
    assert row["positive_residues_all_excluded"] is False
    assert row["open_positive_residues"] == [9]
    assert row["residue_zero_remains"] is True
    witness = p23_u9_open_profile()
    assert witness["total_deficit"] == 380
    assert witness["pair_slack"] == 0
    assert witness["restored_floor_plus_two_cell"]["reduced_size"] == 3
    assert witness["restored_floor_plus_two_cell"]["floor_plus_two_forbidden"] is False


def test_exact_residue_zero_census_is_stable_and_fingerprinted():
    row = p23_residue_zero_profile_census()
    assert row["phase_zero_row_count"] == 426
    assert row["phase_one_row_count"] == 11
    assert row["profile_count"] == 1247
    assert row["distinct_global_shape_count"] == 485
    assert row["pair_slack_histogram"] == {
        0: 363,
        4: 264,
        8: 189,
        12: 136,
        16: 94,
        20: 68,
        24: 49,
        28: 35,
        32: 21,
        36: 13,
        40: 7,
        44: 4,
        48: 1,
        52: 1,
        56: 1,
        60: 1,
    }
    assert (
        row["canonical_profile_sha256"]
        == "19ea72e792303d42863d327114eea6edde0abb3039a578b991387ead83fa5cc0"
    )


def test_degree_ten_envelope_excludes_all_arc_profiles():
    row = p23_arc_envelope_exclusion()
    assert row["profile_count"] == 363
    assert row["arc_size"] == 20
    assert row["tau"] == 5
    assert row["envelope_degree"] == 10
    assert row["minimum_high_direction_count"] == 3
    assert row["profiles_with_high_edge"] == 320
    assert row["zero_high_edge_profiles"] == 43
    assert row["complete_arc_classification"]["no_complete_arc_sizes"] == [
        18,
        19,
        20,
        21,
        22,
        23,
    ]
    assert row["excluded"] is True


def test_conic_core_repair_has_slack_twenty_four_off_conic_floor():
    assert line_pair_slack(3) == 4
    assert line_pair_slack(4) == 8
    assert line_pair_slack(5) == 16
    row = conic_core_repair_lemma()
    assert row["classification_threshold"] == 18
    assert row["off_conic_full_secant_minimum"] == 11
    assert min(
        item["pair_slack_floor"]
        for item in row["off_conic_count_rows"].values()
    ) == 24
    assert row["proved"] is True


def test_complete_arc_classification_and_low_slack_counts():
    classification = p23_complete_arc_classification()
    assert classification["doi"] == "10.1002/jcd.20211"
    assert classification["complete_arc_sizes_pg2_23"] == [
        10,
        12,
        13,
        14,
        15,
        16,
        17,
        24,
    ]
    assert classification["complete_arc_counts"][17] == 5
    row = p23_low_slack_conic_exclusion()
    observed = {
        item["pair_slack"]: (
            item["profile_count"],
            item["excluded_profile_count"],
            item["remaining_profile_count"],
        )
        for item in row["rows"]
    }
    assert observed == {
        4: (264, 264, 0),
        8: (189, 189, 0),
        12: (136, 135, 1),
        16: (94, 93, 1),
    }
    assert row["excluded_profile_count"] == 681
    assert row["remaining_low_slack_profile_count"] == 2


def test_theorem_retracts_only_residue_zero_but_retains_conditional_reduction():
    row = p23_reduction_theorem()
    assert row["positive_residues_excluded"] == [2, 3, 4, 5, 6, 8]
    assert row["open_positive_residues"] == [9]
    assert row["only_residue_zero_remains"] is False
    assert row["record_status"] == "OPEN_RETRACTED_REDUCTION"
    assert row["former_only_residue_zero_claim_retracted"] is True
    assert row["residue_zero_reduction_proved_conditionally"] is True
    assert row["residue_zero_profile_count_before"] == 1247
    assert row["residue_zero_profiles_excluded"] == 1044
    assert row["residue_zero_profile_count_after"] == 203
    assert row["remaining_pair_slack_histogram"] == {
        12: 1,
        16: 1,
        20: 68,
        24: 49,
        28: 35,
        32: 21,
        36: 13,
        40: 7,
        44: 4,
        48: 1,
        52: 1,
        56: 1,
        60: 1,
    }
    assert len(row["exceptional_low_slack_profiles"]) == 2
    assert row["p23_second_all_finite_endpoint_closed"] is False
    assert row["remaining_same_boundary_primes"] == [17, 19, 23]
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is False


def test_open_main_is_non_raising():
    assert prop15684_main()["record_status"] == "OPEN_RETRACTED_REDUCTION"


def test_canonical_evidence_matches_retracted_source_record():
    source = p23_reduction_theorem()
    canonical = json.loads(
        (ROOT / "evidence/e1_gmin_m4_prop15684.json").read_text()
    )
    assert canonical["record_status"] == source["record_status"]
    assert canonical["proved"] is source["proved"] is False
    assert canonical["open_positive_residues"] == source["open_positive_residues"]
    assert canonical["only_residue_zero_remains"] is False
    assert canonical["residue_zero_subledger"] == {
        "status": "CONDITIONAL_REUSABLE_REDUCTION",
        "profile_count_before": source["residue_zero_profile_count_before"],
        "profiles_excluded": source["residue_zero_profiles_excluded"],
        "profile_count_after": source["residue_zero_profile_count_after"],
    }
    historical = ROOT / canonical["historical_payload"]
    assert historical.exists()
    assert json.loads(historical.read_text())["record_status"] == (
        "HISTORICAL_RETRACTED_PAYLOAD"
    )
    assert canonical["theorem"]["boundary_gate_status"].startswith(
        "SUPERSEDED_AND_EXCLUDED_BY_15.721"
    )
