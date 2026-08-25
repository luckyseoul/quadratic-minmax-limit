from fractions import Fraction

from src.e1_gmin_m4_prop15642 import (
    certificate_is_exact,
    infinity_finite_boundary_consequence,
    nonbaseline_scaled_cost_floor,
    nonzero_quadratic_mass_floor,
    polynomial_distance_support_floor,
    stabilizer_mass_certificate,
    theorem_quadratic_lift_mass,
)


def test_stabilizer_moment_certificate_is_exact():
    for p in range(3, 104, 2):
        assert certificate_is_exact(p)
        result = stabilizer_mass_certificate(p)
        if p % 4 == 3:
            assert result["value"] == Fraction(1, p)
        else:
            r = (p - 1) // 4
            assert result["value"] == Fraction(r, (r + 1) * p)


def test_exact_polynomial_distance_support_and_combined_mass_floor():
    for p in range(5, 104, 2):
        assert polynomial_distance_support_floor(p) == Fraction(
            p * p - 1, 16 * p * (p - 2)
        )
        assert nonzero_quadratic_mass_floor(p) >= Fraction(
            p * p - 1, 16 * p * (p - 2)
        )


def test_even_scaled_cost_floor_and_uniform_exception_bound():
    expected = {3: 4, 5: 2, 7: 4, 11: 4, 13: 4, 17: 6, 19: 6,
                23: 8, 29: 8, 31: 10, 41: 12, 101: 26}
    for p, value in expected.items():
        assert nonbaseline_scaled_cost_floor(p) == value


def test_infinity_finite_boundary_rigidity_and_exception_count():
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        plus = infinity_finite_boundary_consequence(p, 1)
        assert plus["baseline"] == "x_s"
        assert plus["type_surplus"] == 0
        assert plus["maximum_nonzero_lifts_per_type"] == 0
        assert plus["pointwise_rigid"] is True

        minus = infinity_finite_boundary_consequence(p, -1)
        assert minus["baseline"] == "1-x_s"
        assert minus["type_surplus"] == p + 1
        assert minus["maximum_nonzero_lifts_per_type"] == (
            (p + 1) // nonbaseline_scaled_cost_floor(p)
        )
        if p >= 5:
            assert minus["maximum_nonzero_lifts_per_type"] <= 3
        assert minus["pointwise_rigid"] is False


def test_theorem_keeps_live_gates_open():
    theorem = theorem_quadratic_lift_mass()
    assert theorem["proved"] is True
    assert theorem["closes_infinity_finite_boundary"] is False
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert theorem["L_status"] == "OPEN"
