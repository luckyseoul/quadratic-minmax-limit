import json
from pathlib import Path

from e1_gmin_m4_prop15676 import (
    canonical_conic_profile_audit,
    conic_classification_ledger,
    external_conic_floor_exclusion,
    tangent_conic_exclusion,
    tangent_type_residue_audit,
    theorem_record,
)


def test_canonical_conics_have_the_two_classified_profiles():
    for p in (17, 19, 23, 29, 31):
        row = canonical_conic_profile_audit(p)
        m = (p + 1) // 2
        assert row["tangent_line_at_infinity"]["profile"] == {"1": p, str(p): 1}
        assert row["external_line_at_infinity"]["profile"] == {
            "1": m + 1,
            "3": m - 1,
        }
        assert row["tangent_line_at_infinity"]["pair_deficit"] == p * (p - 1)
        assert row["external_line_at_infinity"]["pair_deficit"] == p * (p - 1)


def test_external_conic_profiles_exceed_the_type_capacity():
    for p in (17, 19, 23, 29, 101):
        zero = external_conic_floor_exclusion(p, 0)
        one = external_conic_floor_exclusion(p, 1)
        assert zero["maximum_b3_per_type"] == 0
        assert one["maximum_b3_per_type"] == 1
        assert zero["global_b3_directions"] > zero["two_type_capacity"]
        assert one["global_b3_directions"] > one["two_type_capacity"]


def test_tangent_conic_four_phase_rows_are_excluded():
    for p in (17, 19, 23, 29, 31, 37, 41, 101):
        for phase in (0, 1):
            row = tangent_conic_exclusion(p, phase)
            assert row["excluded"] is True
            if p % 4 == 1 and phase == 0:
                assert row["method"] == "type floor exceeds budget"
            else:
                assert row["substituted_congruences"] == (
                    "q divides x and q divides y"
                )
                assert row["candidates"] == [
                    {
                        "x": 0,
                        "y": 0,
                        "E": 2 * phase,
                        "I": 4 * p + 1 - 2 * phase,
                        "support_upper": p + 4 * phase,
                        "support_contradiction": True,
                    }
                ]


def test_tangent_type_residues_retain_b1_baselines():
    for p in (17, 19, 23, 29, 101):
        m = (p + 1) // 2
        for phase in (0, 1):
            ordinary = tangent_type_residue_audit(p, phase, False)
            assert [row["u"] for row in ordinary["feasible_residue_rows"]] == [
                0 if phase == 0 else m - 1
            ]
            assert ordinary["minimum_b1_baselines"] >= m - 1

            exceptional = tangent_type_residue_audit(p, phase, True)
            if p % 4 == 1 and phase == 0:
                assert exceptional["feasible_residue_rows"] == []
                assert exceptional["minimum_b1_baselines"] is None
            else:
                assert [
                    row["u"] for row in exceptional["feasible_residue_rows"]
                ] == [0 if phase == 0 else m - 1]
                assert exceptional["minimum_b1_baselines"] >= m - 2


def test_scope_is_equality_only_and_evidence_matches():
    assert conic_classification_ledger()["exhaustive"] is True
    row = theorem_record()
    assert row["proved"] is True
    assert row["theorem"]["pair_deficit_equality"] == (
        "EXCLUDED_FOR_BOTH_PRODUCT_SIGNS"
    )
    assert row["theorem"]["strict_pair_deficit_branch"] == "OPEN"
    assert row["theorem"]["general_residual_ii"] is False
    assert row["theorem"]["R1"] is False

    root = Path(__file__).resolve().parents[1]
    stored = json.loads(
        (root / "evidence" / "e1_gmin_m4_prop15676.json").read_text()
    )
    assert stored == row
