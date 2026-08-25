from src.e1_gmin_m4_prop15632 import field_direction_data, projective_directions
from src.e1_gmin_m4_prop15646 import (
    additive_matrix_total,
    normal_form_signed_audit,
    theorem_negative_branch_exclusion,
)
from src.minmax_quadratic import paley_conference_prime_power


def test_parallel_finite_edge_sign_is_direction_type():
    for p in (5, 7):
        C = paley_conference_prime_power(p)
        for direction in projective_directions(p):
            eps, labels = field_direction_data(p, direction)
            signs = {
                int(round(C[u + 1, v + 1]))
                for u in range(p * p)
                for v in range(u + 1, p * p)
                if labels[u] == labels[v]
            }
            assert signs == {eps}


def test_additive_matrix_total_vanishes_for_zero_sum_deviations():
    for p in range(3, 30, 2):
        ideal = (0,) * p
        transfer = (1, -1) + (0,) * (p - 2)
        for eps in (-1, 1):
            assert additive_matrix_total(ideal, eps) == 0
            assert additive_matrix_total(transfer, eps) == 0


def test_both_exception_splits_contradict_a_baseline_type():
    for p in range(3, 202, 2):
        positive_heavy = normal_form_signed_audit(p, 3, 1)
        assert positive_heavy["total_finite_signed_sum"] == 2
        assert positive_heavy["baseline_transverse_signed_sum"]["-1"] == 4
        assert positive_heavy["contradictory_baseline_type"] == -1
        assert positive_heavy["contradiction"] is True

        negative_heavy = normal_form_signed_audit(p, 1, 3)
        assert negative_heavy["total_finite_signed_sum"] == -2
        assert negative_heavy["baseline_transverse_signed_sum"]["+1"] == -4
        assert negative_heavy["contradictory_baseline_type"] == 1
        assert negative_heavy["contradiction"] is True


def test_theorem_closes_only_the_asymptotic_negative_two_point_branch():
    theorem = theorem_negative_branch_exclusion()
    assert theorem["proved"] is True
    assert theorem["all_sufficiently_large_odd_primes"] is True
    assert theorem["closes_negative_product_infinity_point_branch"] is True
    assert theorem["closes_all_infinity_point_boundaries"] is False
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert theorem["L_status"] == "OPEN"
