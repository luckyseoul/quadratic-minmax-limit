import json
from pathlib import Path

from e1_gmin_m4_prop15678 import (
    conic_secant_survival_ledger,
    endpoint_profiles,
    p17_arc_classification_ledger,
    theorem_record,
    three_undetermined_direction_contradiction,
    type_residue_ledger,
    u2_lift_exclusion,
    u3_coefficient_exclusion,
)


def test_exact_p17_residue_minima_and_local_exclusions():
    ledger = type_residue_ledger()
    zero = ledger["phase_residue_rows"]["0"]
    one = ledger["phase_residue_rows"]["1"]
    assert {u: row["minimum_deficit"] for u, row in zero.items()} == {
        "0": 84,
        "2": 82,
        "3": 84,
        "4": 96,
        "5": 98,
        "6": 110,
        "7": 112,
        "8": 112,
    }
    assert {u: row["minimum_deficit"] for u, row in one.items()} == {"8": 96}
    assert u2_lift_exclusion()["nonzero_B_scaled_cost_floor"] == 6
    assert u2_lift_exclusion()["excluded"] is True
    u3 = u3_coefficient_exclusion()
    assert [row["ell"] for row in u3["rows"]] == [0, 2, 4]
    assert all(row["excluded_by_l1"] for row in u3["rows"])


def test_pair_slack_leaves_exactly_two_arc_profiles():
    row = endpoint_profiles()
    assert row["common_global_secant_distribution"] == {
        "0": 3,
        "1": 1,
        "6": 8,
        "7": 6,
    }
    assert row["common_undetermined_directions"] == 3
    assert row["profiles"]["A"]["phase_profiles_b"] == {
        "0": {"0": 6, "14": 3},
        "1": {"2": 8, "12": 1},
    }
    assert row["profiles"]["B"]["phase_profiles_b"] == {
        "0": {"0": 6, "12": 1, "14": 2},
        "1": {"2": 8, "14": 1},
    }
    assert all(
        profile["total_deficit"] == 182 and profile["arc"] is True
        for profile in row["profiles"].values()
    )


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
    assert three_undetermined_direction_contradiction()["excluded"] is True


def test_scope_and_generated_evidence_match():
    row = theorem_record()
    assert row["proved"] is True
    assert row["theorem"]["p17_first_all_finite_survivor"] == "EXCLUDED_HERE"
    assert row["theorem"]["all_odd_primes_p_at_least_17"] == (
        "FIRST_SURVIVOR_EXCLUDED"
    )
    assert row["theorem"]["later_all_finite_boundary_sizes"] == "OPEN"
    assert row["theorem"]["general_residual_ii"] is False
    assert row["theorem"]["R1"] is False

    root = Path(__file__).resolve().parents[1]
    stored = json.loads(
        (root / "evidence" / "e1_gmin_m4_prop15678.json").read_text()
    )
    assert stored == row
