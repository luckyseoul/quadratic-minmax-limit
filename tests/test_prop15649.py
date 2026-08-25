import itertools
from math import comb

from src.e1_gmin_m4_prop15649 import (
    CERTIFICATE_ARCHIVE_SHA256,
    fixed_star_certificate,
    lift_classification,
    reconstruction_holds,
    theorem_balanced_p7_exclusion,
)


def test_pair_incidence_reconstruction_on_support_five_lift():
    support = {
        (0, 1, 2, 4),
        (0, 2, 3, 4),
        (0, 2, 4, 5),
        (0, 2, 4, 6),
        (1, 3, 5, 6),
    }
    values = {
        X: (2 if X in support else 0)
        for X in itertools.combinations(range(7), 4)
    }
    assert reconstruction_holds(values)


def test_complete_lift_histograms():
    lift = lift_classification()
    assert lift["possible_support_sizes"] == [5, 8, 9, 10]
    assert lift["maximum_value"] == 2
    assert lift["total_labelled_vectors"] == 1764
    assert lift["enumeration_complete"] is True


def test_fixed_star_coverage_is_exact():
    stars = fixed_star_certificate()
    assert stars["all_stars_per_pair"] == comb(49, 3) == 18424
    assert stars["total_orbits"] == 6076
    assert stars["main_infeasible"] == 6049
    assert stars["retried_infeasible"] == 27
    assert stars["final_unknown"] == stars["final_feasible"] == 0
    assert stars["complete"] is True


def test_theorem_closes_p7_but_not_larger_problem():
    theorem = theorem_balanced_p7_exclusion()
    assert theorem["proved"] is True
    assert theorem["all_p7_negative_two_point_profiles_closed"] is True
    assert theorem["remaining_negative_two_point_cases"] == ["p=5"]
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert len(CERTIFICATE_ARCHIVE_SHA256) == 64
