from collections import Counter

from scripts.p5_negative_profile_cpsat import count_candidates
from scripts.p5_negative_symmetry_certificate import CASES, validate_case_cover
from src.e1_gmin_m4_prop15650 import (
    CERTIFICATE_ARCHIVE_SHA256,
    p5_arithmetic_profiles,
    p5_placement_orbit_count,
    p5_type_excess_profiles,
    theorem_p5_negative_two_point_exclusion,
)


def test_two_and_only_two_type_excess_profiles():
    assert p5_type_excess_profiles() == [(0, 0, 6), (2, 2, 2)]


def test_arithmetic_profiles_match_finite_model():
    profiles = p5_arithmetic_profiles()
    assert profiles == count_candidates()
    assert len(profiles) == 24
    kinds = Counter(
        (row["positive_profile"], row["negative_profile"]) for row in profiles
    )
    assert kinds == {
        ("unique", "unique"): 9,
        ("unique", "distributed"): 5,
        ("distributed", "unique"): 5,
        ("distributed", "distributed"): 5,
    }


def test_placement_orbit_cover():
    orbit_count = p5_placement_orbit_count()
    assert orbit_count["placement_orbits"] == len(CASES) == 33
    assert validate_case_cover() is True


def test_theorem_closes_negative_two_point_branch_only():
    theorem = theorem_p5_negative_two_point_exclusion()
    assert theorem["proved"] is True
    assert theorem["remaining_negative_two_point_cases"] == []
    assert theorem["closes_negative_product_infinity_point_branch_all_primes"] is True
    assert theorem["closes_all_infinity_point_boundaries"] is False
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert len(CERTIFICATE_ARCHIVE_SHA256) == 64
