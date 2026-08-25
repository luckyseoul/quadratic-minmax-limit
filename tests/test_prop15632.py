import random
from fractions import Fraction

from src.e1_gmin_m4_prop15632 import (
    AFFINE_P5_G,
    affine_boundary_budget,
    affine_p5_counterexample_to_affine_close,
    direct_affine_slack_audit,
    eulerian_residual_type_budget_gap,
    middle_weight,
    parity_majorant_floor,
    scaled_direction_floor,
    theorem_affine_parity_budget,
)


def test_exact_quadratic_parity_majorants_and_integral_floors():
    expected = {
        (7, 3, 0): Fraction(4, 7),
        (7, 3, 1): Fraction(1),
        (11, 3, 0): Fraction(8, 11),
        (11, 5, 1): Fraction(9, 11),
        (17, 5, 0): Fraction(1),
        (17, 5, 1): Fraction(1),
    }
    for key, value in expected.items():
        result = parity_majorant_floor(*key)
        assert result["proved"] is True
        assert result["value"] == value
        assert scaled_direction_floor(*key) == 2 * ((key[0] * value).__ceil__())


def test_majorant_complement_symmetry():
    for p in (5, 7, 11, 13):
        phase_shift = middle_weight(p) & 1
        for b in range(p + 1):
            for phase in (0, 1):
                lhs = parity_majorant_floor(p, p - b, phase)["value"]
                rhs = parity_majorant_floor(p, b, phase ^ phase_shift)["value"]
                assert lhs == rhs


def test_directional_sum_and_evenness_on_arbitrary_odd_edge_sets():
    rng = random.Random(15632)
    for p, h in ((5, 9), (5, 21), (7, 15), (7, 29)):
        n = p * p + 1
        edges = [(a, b) for a in range(n) for b in range(a + 1, n)]
        H = tuple(rng.sample(edges, h))
        result = affine_boundary_budget(p, H)
        assert result["sum_identity"] is True
        assert result["type_sum_identity"] is True
        assert result["all_a_even"] is True
        assert result["sum_a"] == (p + 1) * (h - 3 * p)
        assert result["type_sums"] == {
            "-1": (p + 1) * (h - 3 * p) // 2,
            "1": (p + 1) * (h - 3 * p) // 2,
        }


def test_boundary_parity_and_mean_are_pointwise_exact():
    candidate = AFFINE_P5_G + ((0, 1),)
    candidate_audit = direct_affine_slack_audit(5, candidate)
    assert candidate_audit["all_slacks_integral"] is True
    assert candidate_audit["all_boundary_parities"] is True
    assert candidate_audit["all_mean_identities"] is True

    # The identities do not rely on nonnegativity or on the separator.
    rng = random.Random(632)
    for p, h in ((5, 11), (7, 15)):
        n = p * p + 1
        edges = [(a, b) for a in range(n) for b in range(a + 1, n)]
        arbitrary_H = tuple(rng.sample(edges, h))
        arbitrary_audit = direct_affine_slack_audit(p, arbitrary_H)
        assert arbitrary_audit["all_slacks_integral"] is True
        assert arbitrary_audit["all_boundary_parities"] is True
        assert arbitrary_audit["all_mean_identities"] is True


def test_p5_affine_solution_prevents_false_close():
    result = affine_p5_counterexample_to_affine_close()
    assert result["proved"] is True
    assert result["a_by_direction"] == [12, 4, 0, 6, 10, 4]
    assert result["a_by_quadratic_type"] == {"-1": 18, "1": 18}
    assert result["slack_supports_by_direction"] == [
        (0, 1, 2, 3, 4),
        (0, 2),
        (0,),
        (0, 2),
        (0, 2, 4),
        (0, 2),
    ]
    assert result["boundary_is_infinity_plus_affine_line"] is True
    assert result["direct_affine_audit"]["all_slacks_nonnegative"] is True
    budget = affine_boundary_budget(5, AFFINE_P5_G + ((0, 1),))
    assert budget["sum_a"] == 36
    assert budget["expected_sum"] == 36
    assert budget["parity_floor_sum"] == 4
    assert budget["type_parity_floor_sums"] == {"-1": 4, "1": 0}
    assert budget["type_parity_budgets_hold"] is True


def test_empty_boundary_residual_branch_is_excluded():
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        assert eulerian_residual_type_budget_gap(p) == (p * p - 1) // 2
        assert eulerian_residual_type_budget_gap(p) > 0

    # A concrete odd Eulerian graph at residual size fails the necessary
    # type-split budget, independently checking the symbolic obstruction.
    p = 5
    H = tuple((i, i + 1) for i in range(20)) + ((0, 20),)
    budget = affine_boundary_budget(p, H)
    assert budget["h"] == 4 * p + 1
    assert budget["rows"][0]["boundary"] == ()
    assert budget["type_parity_budgets_hold"] is False


def test_theorem_is_structural_and_keeps_live_gates_open():
    theorem = theorem_affine_parity_budget()
    assert theorem["proved"] is True
    assert theorem["eulerian_residual_boundary_excluded"] is True
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert theorem["L_status"] == "OPEN"
