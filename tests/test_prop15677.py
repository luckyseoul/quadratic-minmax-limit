import json
from pathlib import Path

from e1_gmin_m4_prop15677 import (
    arithmetic_rows,
    fibre_profile_classification,
    first_survivor_exclusion,
    phase_residue_reduction,
    remaining_regimes,
    theorem_record,
    zero_quotient_lift_exclusion,
)


def test_residue_reduction_retains_only_u2_and_the_p1_u3_row():
    for p in (23, 31, 47, 71, 79):
        row = phase_residue_reduction(p)
        assert row["phase_zero_candidate_residues"] == [2]
        assert row["phase_one_only_residue"] == (p - 1) // 2
        assert row["proved"] is True
    for p in (41, 73, 89, 97):
        row = phase_residue_reduction(p)
        assert row["phase_zero_candidate_residues"] == [2, 3]
        assert row["u3_increment"] == 2
        assert row["proved"] is True


def test_xnor_arithmetic_leaves_exactly_l2_and_l4():
    for p in (23, 31, 41, 47, 71, 73, 79, 89):
        row = remaining_regimes(p)
        assert row["normal_form"] == [
            {"u0": 2, "l": 2, "I": 2 * p, "E": 2 * p + 1},
            {"u0": 2, "l": 4, "I": p - 1, "E": 3 * p + 2},
        ]
        assert row["proved"] is True


def test_extra_p1_residue_is_excluded_in_all_three_rows():
    for p in (41, 73, 89, 97):
        row = arithmetic_rows(p, 3)
        assert row["j_candidates"] == [3]
        assert [item["l"] for item in row["rows"]] == [0, 2, 4]
        assert all(item["excluded_by_l1"] for item in row["rows"])


def test_remaining_infinity_fibre_profiles_are_near_perfect():
    for p in (23, 31, 41, 73, 103):
        l2 = fibre_profile_classification(p, 2)
        assert l2["allowed_fibre_count_histograms"] == [
            {"2": p},
            {"1": 1, "2": p - 2, "3": 1},
        ]
        l4 = fibre_profile_classification(p, 4)
        assert l4["allowed_fibre_count_histograms"] == [
            {"0": 1, "1": p - 1},
            {"0": 2, "1": p - 3, "2": 1},
        ]
        assert l2["proved"] is l4["proved"] is True


def test_nonzero_lift_mass_closes_every_pre_lift_regime():
    for p in (23, 31, 41, 47, 71, 73, 79, 89):
        assert zero_quotient_lift_exclusion(p, 2)["excluded"] is True
    for p in (41, 73, 89, 97):
        assert zero_quotient_lift_exclusion(p, 3)["excluded"] is True

    for p in (23, 31, 41, 47, 71, 73, 79, 89):
        assert first_survivor_exclusion(p)["excluded"] is True


def test_scope_is_the_first_survivor_only_and_evidence_matches():
    row = theorem_record()
    assert row["proved"] is True
    assert row["theorem"]["p_mod_8_1_7_p_at_least_23"] == "CLOSED_HERE"
    assert row["theorem"]["all_odd_primes_p_at_least_19"] == (
        "FIRST_SURVIVOR_EXCLUDED"
    )
    assert row["theorem"]["p17_endpoint"] == "OPEN_ADDITIONAL_U0_ZERO_ROW"
    assert row["theorem"]["general_residual_ii"] is False
    assert row["theorem"]["R1"] is False

    root = Path(__file__).resolve().parents[1]
    stored = json.loads(
        (root / "evidence" / "e1_gmin_m4_prop15677.json").read_text()
    )
    assert stored == row
