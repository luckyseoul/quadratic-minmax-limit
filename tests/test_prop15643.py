from src.e1_gmin_m4_prop15643 import (
    populated_direction_necessary,
    positive_product_arithmetic,
    theorem_positive_product_boundary,
)


def test_populated_direction_inequality_dies_at_17():
    for p in range(17, 202, 2):
        for k0 in range(9):
            for kd in range(1, 9 - k0):
                assert populated_direction_necessary(p, k0, kd) is False


def test_only_endpoint_is_forbidden_all_infinity_star():
    for p in range(17, 202, 2):
        result = positive_product_arithmetic(p)
        assert result["parallel_counts_are_q_multiples"] is True
        assert result["no_populated_direction"] is True
        assert result["only_arithmetic_endpoint"]["k0"] == 8
        assert result["endpoint_is_all_infinity_star"] is True
        assert result["endpoint_boundary_size"] == 4 * p + 2
        assert result["positive_product_infinity_point_boundary_excluded"] is True


def test_small_ranges_are_not_soft_closed():
    for p in (5, 7, 11, 13):
        result = positive_product_arithmetic(p)
        assert result["positive_product_infinity_point_boundary_excluded"] is False
        assert any(row["allowed_positive_kd"] for row in result["rows"])


def test_theorem_keeps_remaining_gates_open():
    theorem = theorem_positive_product_boundary()
    assert theorem["proved"] is True
    assert theorem["p_5_7_11_13_status"] == "OPEN"
    assert theorem["negative_product_status"] == "OPEN"
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert theorem["L_status"] == "OPEN"
