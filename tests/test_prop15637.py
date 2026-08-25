import pytest

from src.e1_gmin_m4_prop15637 import (
    candidate_energy,
    candidate_common_sum_magnitudes,
    candidate_scaled_norm,
    dense_branch_recurrence_certificate,
    dense_candidate_branches_excluded,
    doubled_energy_six_cubic_defect,
    doubled_energy_six_cubic_factor,
    energy_four_cubic_defect,
    energy_four_cubic_factor,
    mds_allowed_active_counts_before_one_profile_kill,
    one_profile_candidate_excluded,
    one_profile_mass_patterns,
    profile_balancing_gap,
    remaining_candidate_active_counts,
    two_double_ode_descent,
    zero_common_sum_candidate_excluded,
    zero_common_sum_gap_theorem,
)


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29, 31, 37, 41, 43])
def test_only_two_dense_active_counts_survive_at_candidate_energy(p):
    k = (p - 1) // 2
    R = k + 1
    assert candidate_energy(p) == p + 3
    assert candidate_scaled_norm(p) == 2 * (p + 3)
    assert mds_allowed_active_counts_before_one_profile_kill(p) == (
        1,
        R - 1,
        R,
    )
    assert one_profile_candidate_excluded(p)
    assert remaining_candidate_active_counts(p) == (R - 1, R)


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23])
def test_one_profile_multiplicity_dichotomy(p):
    one_double, two_double = one_profile_mass_patterns(p)
    assert one_double["energy_defect"] == 2
    assert one_double["magnitude_two_entries"] == 1
    assert one_double["distinct_support"] == p
    assert one_double["root_polynomial_difference_degree_at_most"] == 1
    assert two_double["energy_defect"] == 4
    assert two_double["magnitude_two_entries"] == 2
    assert two_double["distinct_support"] == p - 3
    assert two_double["root_polynomial_difference_degree_at_most"] == 0


@pytest.mark.parametrize("p", [11, 13, 17, 19, 29])
def test_two_double_formal_series_descent_indices(p):
    cert = two_double_ode_descent(p)
    k = (p - 1) // 2
    assert cert["known_zero_gap"] == [k + 1, 2 * k - 1]
    assert cert["highest_series_index_used"] == k + 4
    assert cert["gap_reaches_first_equation"]
    assert cert["nonzero_descent_multipliers"]
    assert cert["conclusion"] == "u_1=...=u_k=0, hence N=D"


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29])
def test_dense_branch_zero_counts_force_the_moment_identities(p):
    cert = dense_branch_recurrence_certificate(p)
    assert cert["full_h_R"]["forced_identity"]
    assert cert["near_full_h_R_minus_1_two_energy_four"]["forced_identity"]
    assert cert["near_full_h_R_minus_1_one_energy_six"][
        "forced_cubic_identity"
    ]
    assert cert["near_full_h_R_minus_1_one_energy_six"][
        "forced_quartic_identity"
    ]
    assert cert["nonzero_characteristic_factors"]
    assert cert["newton_degree_four_valid"]
    assert dense_candidate_branches_excluded(p)
    assert zero_common_sum_candidate_excluded(p)


def test_exceptional_profile_defect_factorizations():
    for values in ((0, 1, 2, 3), (1, 4, 7, 9), (-3, 2, 5, 11)):
        assert energy_four_cubic_defect(*values) == energy_four_cubic_factor(
            *values
        )
    for values in ((0, 1, 2), (3, -2, 7), (11, 4, -5)):
        assert doubled_energy_six_cubic_defect(
            *values
        ) == doubled_energy_six_cubic_factor(*values)


def test_theorem_scope_kills_zero_sum_but_keeps_other_tail_open():
    theorem = zero_common_sum_gap_theorem()
    assert theorem["proved"]
    assert "no zero-common-sum profile exists" in theorem["scope"]
    assert all(
        row["allowed_h_after_dense_branch_kills"] == ()
        for row in theorem["rows"].values()
    )


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29, 31, 43])
def test_balancing_leaves_exactly_three_nonzero_common_sums(p):
    assert candidate_common_sum_magnitudes(p) == (0, 2, p - 1, p + 1)
    assert profile_balancing_gap(p, 2 * p) > candidate_scaled_norm(p)
