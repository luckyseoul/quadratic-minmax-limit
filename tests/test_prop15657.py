from fractions import Fraction

from src.e1_gmin_m4_prop15632 import scaled_direction_floor
from src.e1_gmin_m4_prop15657 import (
    finite_fibre_partition_rows,
    infinity_size_six_exclusion,
    large_boundary_floor_certificate,
    minimum_type_deficit,
    no_infinity_size_six_exclusion,
    pair_deficit_budget,
    size_six_floor_table,
    theorem_size_six_boundary_pge11,
)


def test_large_boundary_floors_have_exact_positive_quadrature_certificates():
    for p in range(11, 100, 2):
        for b in (5, 6):
            for phase in (0, 1):
                row = large_boundary_floor_certificate(p, b, phase)
                assert row["exact_positive_quadrature_certificate"] is True
                assert row["scaled_floor"] == scaled_direction_floor(p, b, phase)
                assert sum(row["quadrature_weights"], Fraction(0)) == 1


def test_size_six_floor_table_has_the_symbolic_formulas():
    for p in range(11, 100, 2):
        phase_one_large = 3 * (p - 5) if p <= 15 else 2 * p
        assert size_six_floor_table(p) == {
            0: (0, 2 * p),
            1: (p + 1, p - 1),
            2: (p + 1, p - 1),
            3: (2 * p - 6, 2 * p),
            4: (2 * p - 6, 2 * p),
            5: (2 * p, phase_one_large),
            6: (2 * p, phase_one_large),
        }


def test_pair_deficit_bound_is_verified_for_every_fibre_partition():
    for s, expected_budget in ((5, 20), (6, 30)):
        rows = finite_fibre_partition_rows(s)
        assert rows
        assert all(row["inequality_verified"] for row in rows)
        result = pair_deficit_budget(s)
        assert result["total_deficit_upper_bound"] == expected_budget
        assert result["pair_count"] * 2 == expected_budget


def test_p11_type_split_requires_too_much_pair_deficit():
    phase_zero = minimum_type_deficit(11, 0, 6)
    phase_one = minimum_type_deficit(11, 1, 6)
    assert phase_zero["minimum_deficit"] == 18
    assert phase_one["minimum_deficit"] == 20
    for c_h in (-1, 1):
        result = no_infinity_size_six_exclusion(11, c_h)
        assert result["required_total_deficit"] == 38
        assert result["pair_deficit_budget"] == 30
        assert result["contradiction_gap"] == 8
        assert result["excluded"] is True


def test_analytic_pair_budget_gaps_close_every_larger_prime_regime():
    for p in (13, 15, 17, 19, 23, 29, 31, 37, 41):
        for c_h in (-1, 1):
            no_infinity = no_infinity_size_six_exclusion(p, c_h)
            assert no_infinity["contradiction_gap"] == p * p - 12 * p + 7
            assert no_infinity["excluded"] is True
    for p in (11, 13, 15, 17, 19, 23, 29, 31, 37, 41):
        for c_h in (-1, 1):
            with_infinity = infinity_size_six_exclusion(p, c_h)
            assert with_infinity["contradiction_gap"] == p * p - 9 * p + 10
            assert with_infinity["excluded"] is True


def test_theorem_closes_exactly_the_claimed_boundary_shell():
    theorem = theorem_size_six_boundary_pge11()
    assert theorem["proved"] is True
    assert theorem["six_point_boundary_all_odd_primes_p_at_least_11"] == "CLOSED"
    assert theorem["p5_size_six"] == "OPEN"
    assert theorem["p7_size_six"] == "OPEN"
    assert theorem["larger_boundary_sizes"] == "OPEN"
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert theorem["L_status"] == "OPEN"
