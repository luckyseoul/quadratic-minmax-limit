from src.e1_gmin_m4_prop15648 import (
    p13_l1_minimum,
    p13_zero_baseline_l1,
    theorem_finite_negative_profiles,
    type_preserving_exception_pair_orbits,
)


def test_p13_zero_baseline_l1_exceeds_edge_budget():
    values = [p13_zero_baseline_l1(u) for u in range(1, 11)]
    assert values == [174, 152, 132, 114, 98, 84, 72, 62, 54, 48]
    result = p13_l1_minimum()
    assert result["minimum"] == 48
    assert result["minimizers"] == [10]
    assert result["transverse_edge_budget"] == 44
    assert result["contradiction"] is True


def test_exception_pair_symmetry_orbits_cover_every_opposite_pair():
    expected = {
        7: [([0, 1], 8), ([0, 3], 8)],
        11: [([0, 1], 12), ([0, 2], 12), ([0, 3], 12)],
        13: [([0, 1], 14), ([0, 4], 14), ([0, 5], 14), ([0, 13], 7)],
    }
    for p, rows_expected in expected.items():
        rows = type_preserving_exception_pair_orbits(p)
        assert [(row["representative"], row["size"]) for row in rows] == rows_expected
        pairs = [tuple(pair) for row in rows for pair in row["pairs"]]
        assert len(pairs) == len(set(pairs)) == ((p + 1) // 2) ** 2


def test_theorem_records_only_the_finite_cases_actually_closed():
    theorem = theorem_finite_negative_profiles()
    assert theorem["proved"] is True
    assert theorem["p13_negative_two_point_closed"] is True
    assert theorem["p11_negative_two_point_finitely_certified"] is True
    assert theorem["p7_unbalanced_profiles_finitely_certified"] == [
        [0, 3],
        [0, 6],
        [3, 0],
        [6, 0],
    ]
    assert theorem["remaining_negative_two_point_cases"] == [
        "p=5 (no guaranteed baseline per type)",
        "p=7 baseline counts (3,3)",
    ]
    assert theorem["closes_negative_product_infinity_point_branch_all_primes"] is False
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert theorem["L_status"] == "OPEN"
