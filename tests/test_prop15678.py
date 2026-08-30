import json
from pathlib import Path

from e1_gmin_m4_prop15678 import (
    conic_secant_survival_ledger,
    endpoint_profiles,
    main as prop15678_main,
    p17_arc_classification_ledger,
    theorem_record,
    three_undetermined_direction_contradiction,
    type_residue_ledger,
    u2_lift_exclusion,
    u3_coefficient_exclusion,
)


ROOT = Path(__file__).resolve().parents[1]


def test_exact_p17_residue_minima_and_local_exclusions():
    ledger = type_residue_ledger()
    zero = ledger["phase_residue_rows"]["0"]
    one = ledger["phase_residue_rows"]["1"]
    assert {u: row["minimum_deficit"] for u, row in zero.items()} == {
        "0": 68,
        "2": 82,
        "3": 84,
        "4": 96,
        "5": 98,
        "6": 108,
        "7": 112,
        "8": 112,
    }
    assert {u: row["minimum_deficit"] for u, row in one.items()} == {"8": 96}
    assert u2_lift_exclusion()["nonzero_B_scaled_cost_floor"] == 6
    assert u2_lift_exclusion()["excluded"] is True
    u3 = u3_coefficient_exclusion()
    assert [row["ell"] for row in u3["rows"]] == [0, 2, 4]
    assert all(row["excluded_by_l1"] for row in u3["rows"])


def test_corrected_pair_slack_census_retains_the_two_legacy_arc_profiles():
    row = endpoint_profiles()
    assert row["phase_zero_profile_count_under_pair_cap"] == 50
    assert row["phase_one_profile_count_under_pair_cap"] == 8
    assert row["candidate_count"] == 108
    assert row["pair_slack_histogram"] == {0: 47, 4: 32, 8: 18, 12: 8, 16: 3}
    assert row["arc_profile_count"] == 47
    assert row["three_undetermined_arc_profile_count"] == 14
    assert row["legacy_common_global_secant_distribution"] == {
        "0": 3,
        "1": 1,
        "6": 8,
        "7": 6,
    }
    assert row["legacy_common_undetermined_directions"] == 3
    assert row["legacy_arc_profiles"]["A"]["phase_profiles_b"] == {
        "0": {"0": 6, "14": 3},
        "1": {"2": 8, "12": 1},
    }
    assert row["legacy_arc_profiles"]["B"]["phase_profiles_b"] == {
        "0": {"0": 6, "12": 1, "14": 2},
        "1": {"2": 8, "14": 1},
    }
    assert all(
        profile["total_deficit"] == 182 and profile["arc"] is True
        for profile in row["legacy_arc_profiles"].values()
    )
    example = row["newly_admitted_floor_plus_two_example"]
    assert example["phase_profiles_b"] == {
        "0": {"0": 4, "2": 1, "12": 4},
        "1": {"2": 8, "12": 1},
    }
    assert example["phase_deficits"] == {"0": 76, "1": 98}
    assert example["pair_slack"] == 8
    assert row["all_profiles_are_arcs"] is False


def test_imported_16_arc_class_and_conic_secant_exit_are_explicit():
    classification = p17_arc_classification_ledger()
    assert classification["external_dependency"] is True
    assert classification["pgl_classes_in_pg2_17"] == {
        "14": 4,
        "15": 1,
        "16": 1,
        "17": 1,
        "18": 1,
    }
    assert classification["every_16_arc_is_conic_contained"] is True
    secants = conic_secant_survival_ledger()
    assert secants["external_point"]["remaining_S_secants_at_least"] == 4
    assert secants["internal_point"]["remaining_S_secants_at_least"] == 5
    geometry = three_undetermined_direction_contradiction()
    assert geometry["corrected_arc_profiles_excluded"] == 14
    assert geometry["corrected_endpoint_ledger_excluded"] is False
    assert geometry["excluded"] is True


def test_scope_is_explicitly_open_and_canonical_evidence_is_not_regenerated():
    row = theorem_record()
    assert row["proved"] is False
    assert row["record_status"] == "OPEN_RETRACTED_REDUCTION"
    assert row["former_claim_retracted"] is True
    assert row["retained_sublemmas_proved"] is True
    assert row["profiles_not_covered_by_retained_geometry"] == 94
    assert row["theorem"]["p17_first_all_finite_survivor"] == "OPEN"
    assert row["theorem"]["all_odd_primes_p_at_least_17"] == (
        "NOT_PROVED_BY_THIS_CHAIN"
    )
    assert row["theorem"]["later_all_finite_boundary_sizes"] == "OPEN"
    assert row["theorem"]["general_residual_ii"] is False
    assert row["theorem"]["R1"] is False


def test_open_main_is_non_raising():
    assert prop15678_main()["record_status"] == "OPEN_RETRACTED_REDUCTION"


def test_canonical_evidence_matches_retracted_source_record():
    source = theorem_record()
    canonical = json.loads(
        (ROOT / "evidence/e1_gmin_m4_prop15678.json").read_text()
    )
    assert canonical["record_status"] == source["record_status"]
    assert canonical["proved"] is source["proved"] is False
    assert canonical["former_claim_retracted"] is True
    assert canonical["corrected_census"] == {
        "compatible_profile_count": source["endpoint_profiles"]["candidate_count"],
        "arc_profile_count": source["endpoint_profiles"]["arc_profile_count"],
        "retained_geometry_excluded_arc_profiles": source[
            "geometry_exclusion_of_three_undetermined_arcs"
        ]["corrected_arc_profiles_excluded"],
        "profiles_not_covered_by_retained_geometry": source[
            "profiles_not_covered_by_retained_geometry"
        ],
    }
    historical = ROOT / canonical["historical_payload"]
    assert historical.exists()
    assert json.loads(historical.read_text())["record_status"] == (
        "HISTORICAL_RETRACTED_PAYLOAD"
    )
    assert canonical["theorem"]["boundary_gate_status"].startswith(
        "SUPERSEDED_AND_EXCLUDED_BY_15.721"
    )
