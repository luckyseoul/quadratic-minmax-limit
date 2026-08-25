from src.e1_gmin_m4_prop15647 import (
    baseline_count_candidates,
    baseline_exists_per_type,
    theorem_negative_two_point_all_prime,
    unique_exception_from_divisibility,
)


def test_sparsity_leaves_a_baseline_in_each_type():
    assert baseline_exists_per_type(5) is False
    for p in range(7, 202, 2):
        assert baseline_exists_per_type(p) is True


def test_same_type_divisibility_forces_one_exception():
    for p in range(7, 202, 2):
        row = unique_exception_from_divisibility(p)
        assert row["positive_partitions"] == [(p + 1,)]
        assert row["unique_exception"] is True
        assert row["exception_a"] == 2 * p


def test_only_small_prime_baseline_count_candidates_survive():
    assert baseline_count_candidates(7) == [
        {
            "positive_baseline": 0,
            "negative_baseline": 3,
            "positive_exception": 1,
            "negative_exception": 4,
            "finite_edges": 14,
            "infinity_edges": 15,
        },
        {
            "positive_baseline": 0,
            "negative_baseline": 6,
            "positive_exception": 1,
            "negative_exception": 7,
            "finite_edges": 26,
            "infinity_edges": 3,
        },
        {
            "positive_baseline": 3,
            "negative_baseline": 0,
            "positive_exception": 4,
            "negative_exception": 1,
            "finite_edges": 14,
            "infinity_edges": 15,
        },
        {
            "positive_baseline": 3,
            "negative_baseline": 3,
            "positive_exception": 4,
            "negative_exception": 4,
            "finite_edges": 26,
            "infinity_edges": 3,
        },
        {
            "positive_baseline": 6,
            "negative_baseline": 0,
            "positive_exception": 7,
            "negative_exception": 1,
            "finite_edges": 26,
            "infinity_edges": 3,
        },
    ]
    assert baseline_count_candidates(11) == [
        {
            "positive_baseline": 0,
            "negative_baseline": 5,
            "positive_exception": 1,
            "negative_exception": 6,
            "finite_edges": 32,
            "infinity_edges": 13,
        },
        {
            "positive_baseline": 5,
            "negative_baseline": 0,
            "positive_exception": 6,
            "negative_exception": 1,
            "finite_edges": 32,
            "infinity_edges": 13,
        },
    ]
    assert baseline_count_candidates(13) == [
        {
            "positive_baseline": 0,
            "negative_baseline": 6,
            "positive_exception": 1,
            "negative_exception": 7,
            "finite_edges": 44,
            "infinity_edges": 9,
        },
        {
            "positive_baseline": 6,
            "negative_baseline": 0,
            "positive_exception": 7,
            "negative_exception": 1,
            "finite_edges": 44,
            "infinity_edges": 9,
        },
    ]
    for p in range(17, 202, 2):
        assert baseline_count_candidates(p) == []


def test_theorem_closes_negative_two_point_only():
    theorem = theorem_negative_two_point_all_prime()
    assert theorem["proved"] is True
    assert theorem["all_odd_primes_at_least_17"] is True
    assert theorem["uses_asymptotic_distance_theorem"] is False
    assert theorem["remaining_two_point_negative_primes"] == [5, 7, 11, 13]
    assert theorem["closes_all_infinity_point_boundaries_p_ge_17"] is False
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert theorem["L_status"] == "OPEN"
