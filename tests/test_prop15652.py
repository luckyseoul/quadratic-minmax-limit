from fractions import Fraction

from src.e1_gmin_m4_prop15632 import scaled_direction_floor
from src.e1_gmin_m4_prop15652 import (
    four_finite_partition_rows,
    infinity_size_four_exclusion,
    no_infinity_size_four_exclusion,
    parity_floor_certificate,
    small_boundary_floor_table,
    theorem_size_four_boundary,
    three_finite_partition_rows,
)


def test_small_boundary_floors_have_exact_positive_quadrature_certificates():
    for p in range(7, 100, 2):
        expected = {
            0: (0, 2 * p),
            1: (p + 1, p - 1),
            2: (p + 1, p - 1),
            3: (2 * p - 6, 2 * p),
            4: (2 * p - 6, 2 * p),
        }
        assert small_boundary_floor_table(p) == expected
        for b in range(5):
            for phase in (0, 1):
                row = parity_floor_certificate(p, b, phase)
                assert row["exact_positive_quadrature_certificate"] is True
                assert row["scaled_floor"] == scaled_direction_floor(p, b, phase)
                assert sum(row["quadrature_weights"], Fraction(0)) == 1


def test_boundary_partition_collision_bookkeeping_is_complete():
    four = four_finite_partition_rows()
    assert {row["partition"] for row in four} == {
        (1, 1, 1, 1),
        (2, 1, 1),
        (2, 2),
        (3, 1),
        (4,),
    }
    assert {row["b"] for row in four} == {0, 2, 4}
    assert all(
        row["pair_collisions"] == sum(n * (n - 1) // 2 for n in row["partition"])
        for row in four
    )
    three = three_finite_partition_rows()
    assert {row["partition"] for row in three} == {
        (1, 1, 1),
        (2, 1),
        (3,),
    }
    assert {row["b"] for row in three} == {1, 3}


def test_four_finite_points_are_excluded_for_every_prime_regime():
    p11 = no_infinity_size_four_exclusion(11)
    p13 = no_infinity_size_four_exclusion(13)
    p17 = no_infinity_size_four_exclusion(17)
    assert p11["reason"] == "at_most_one_good_collision_then_good_type_exceeds_budget"
    assert p11["contradiction_gap"] == 20
    assert p13["reason"] == "all_pairs_consumed_then_good_type_exceeds_budget"
    assert p13["contradiction_gap"] == 42
    assert p17["reason"] == "required_b2_directions_exceed_six_pair_directions"
    assert p17["contradiction_gap"] == 2
    assert all(
        no_infinity_size_four_exclusion(p)["excluded"]
        for p in (11, 13, 17, 19, 23, 29, 31, 37)
    )


def test_infinity_plus_three_finite_points_and_exception_scope():
    assert infinity_size_four_exclusion(7, -1)["excluded"] is True
    assert infinity_size_four_exclusion(7, 1)["excluded"] is False
    for p in (11, 13, 17, 19, 23, 29, 31):
        assert infinity_size_four_exclusion(p, -1)["excluded"] is True
        assert infinity_size_four_exclusion(p, 1)["excluded"] is True


def test_theorem_closes_exactly_the_claimed_boundary_shell():
    theorem = theorem_size_four_boundary()
    assert theorem["proved"] is True
    assert theorem["four_point_boundary_all_odd_primes_p_at_least_11"] == "CLOSED"
    assert theorem["first_open_boundary_size_for_p_at_least_11"] == 6
    assert theorem["p5_size_four"] == "OPEN"
    assert theorem["p7_size_four"].startswith("OPEN")
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
