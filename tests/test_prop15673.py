import itertools
import json
from pathlib import Path

from e1_gmin_m4_prop15673 import (
    branch_arithmetic,
    coefficient_ledger,
    endpoint_branch_exclusion,
    endpoint_floor_ledger,
    endpoint_type_form_ledger,
    p17_xnor_l1_ledger,
    same_type_normal_form,
    symbolic_range_ledger,
    theorem_record,
)


def test_four_endpoint_floor_patterns_and_lift_quantum():
    expected = {
        (19, 0): (20, 20),
        (19, 1): (18, 18),
        (17, 0): (18, 16),
        (17, 1): (16, 18),
    }
    for key, floors in expected.items():
        row = endpoint_floor_ledger(*key)
        assert (row["b1_floor"], row["b_p_minus_2_floor"]) == floors
        assert row["minimum_nonzero_lift_scaled_cost"] >= 4


def test_same_type_normal_forms_use_arc_and_collinear_exits():
    for p, phase in ((19, 0), (17, 1), (19, 1), (17, 0)):
        row = same_type_normal_form(p, phase)
        assert row["proved"] is True
        assert row["forbidden_two_unit_lift"] is True
        assert row["collinear_case"] == "CLOSED_BY_15.671_AND_15.672"
        assert "Segre" in row["equality_case"]
        assert "three collinear points at infinity" in row["equality_case"]


def test_geometry_forces_opposite_endpoint_baseline_types():
    expected = {
        (19, 0): (0, [10]),
        (17, 1): (1, [8, 9]),
        (19, 1): (2, [9, 10, 11]),
        (17, 0): (1, [9, 10]),
    }
    for key, (offset, r_values) in expected.items():
        row = endpoint_type_form_ledger(*key)
        assert row["proved"] is True
        assert row["forced_opposite_baseline_kinds"] is True
        assert row["finite_edge_offset"] == offset
        assert row["admissible_type_pairs"] == [
            {
                "baseline_kinds": ["b=1", "b=p-2"],
                "R_values": r_values,
                "finite_edge_offset": offset,
            }
        ]


def test_baseline_coefficient_targets():
    assert coefficient_ledger(19, 0)["b1_divisibility"] == (
        "q divides I+P_d-5"
    )
    assert coefficient_ledger(19, 1)["b1_divisibility"] == (
        "q divides I+P_d-3"
    )
    assert coefficient_ledger(17, 0)["complement_target"] == "4 + z_a*z_b"
    assert coefficient_ledger(17, 1)["complement_target"] == "4 - z_a*z_b"


def test_large_range_candidate_reductions():
    assert branch_arithmetic(19, 0)["candidates"] == [
        {
            "x": 1,
            "y": 0,
            "E": 10,
            "I": 67,
            "I_boundary_upper": 37,
            "boundary_contradiction": True,
        }
    ]
    assert branch_arithmetic(17, 1)["candidates"] == [
        {
            "x": 0,
            "y": 1,
            "E": 10,
            "I": 59,
            "I_boundary_upper": 35,
            "boundary_contradiction": True,
        }
    ]
    assert branch_arithmetic(19, 1)["candidates"] == []
    assert branch_arithmetic(29, 0)["candidates"] == []
    assert branch_arithmetic(17, 0)["candidates"] == [
        {
            "x": 0,
            "y": 7,
            "E": 64,
            "I": 5,
            "I_boundary_upper": 143,
            "boundary_contradiction": False,
        }
    ]


def test_p17_l1_minimum_independent_composition_check():
    ledger = p17_xnor_l1_ledger()
    assert ledger["exact_minimum"] == 75
    assert ledger["transverse_edge_budget"] == 57
    assert ledger["contradiction"] is True

    direct = None
    # Weak compositions of five into seventeen labelled fibres.  The
    # distinguished pair is (0,1); symmetry makes this lossless.
    for bars in itertools.combinations(range(21), 16):
        previous = -1
        counts = []
        for bar in bars + (21,):
            counts.append(bar - previous - 1)
            previous = bar
        value = sum(
            abs(
                1
                - counts[first]
                - counts[second]
                + int((first, second) == (0, 1))
            )
            for first in range(17)
            for second in range(first + 1, 17)
        )
        direct = value if direct is None else min(direct, value)
    assert direct == 75


def test_both_phases_excluded_uniformly_in_samples():
    for p in (17, 19, 23, 29, 31, 37, 41, 101):
        for phase in (0, 1):
            assert endpoint_branch_exclusion(p, phase)["excluded"] is True
    assert symbolic_range_ledger()["covers_every_odd_p_at_least_17"] is True


def test_theorem_scope_and_generated_evidence():
    row = theorem_record()
    assert row["proved"] is True
    assert row["theorem"]["all_odd_primes_p_at_least_17"] == (
        "EXCLUDED_FOR_BOTH_PRODUCT_SIGNS"
    )
    assert row["theorem"]["nonendpoint_directional_profiles"] == "OPEN"
    assert row["theorem"]["general_residual_ii"] is False
    assert row["theorem"]["R1"] is False

    root = Path(__file__).resolve().parents[1]
    stored = json.loads(
        (root / "evidence" / "e1_gmin_m4_prop15673.json").read_text()
    )
    assert stored == row
