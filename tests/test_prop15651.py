import itertools

from scripts.p7_positive_star_classify import classify
from scripts.positive_two_point_additive_cpsat import (
    allowed_k0_values,
    exact_l1_star_profiles,
)
from src.e1_gmin_m4_prop15651 import (
    CERTIFICATE_ARCHIVE_SHA256,
    P7_ORBIT_CERTIFICATE_SHA256,
    k0_zero_type_capacity,
    p7_finite_coverage,
    small_positive_k0_values,
    theorem_positive_two_point_all_primes,
)


def test_small_prime_arithmetic_cover_is_complete():
    expected = small_positive_k0_values()
    assert expected == {
        5: [0, 1, 2, 3, 4, 5, 8],
        7: [0, 2, 4, 8],
        11: [0, 2, 8],
        13: [0, 1, 8],
    }
    assert {p: allowed_k0_values(p) for p in expected} == expected


def test_exact_l1_profile_eliminations():
    assert exact_l1_star_profiles(7, 4, 1) == []
    assert exact_l1_star_profiles(11, 2, 1) == []
    assert exact_l1_star_profiles(13, 1, 0) == []
    assert len(exact_l1_star_profiles(13, 1, 1)) == 3


def test_additive_coefficients_reconstruct_every_middle_slice_score():
    for p, k0, kd in ((5, 2, 1), (7, 0, 2), (11, 0, 1), (13, 0, 1)):
        q = (p - 1) // 2
        profiles = exact_l1_star_profiles(p, k0, kd)
        assert profiles
        special_count, other_counts = profiles[0]
        counts = (special_count, *other_counts)
        for eps in (-1, 1):
            K = {
                (s, t): eps
                * (
                    k0 + kd + int(s == 0) + int(t == 0)
                    - counts[s] - counts[t]
                )
                for s, t in itertools.combinations(range(p), 2)
            }
            for chosen in itertools.combinations(range(p), (p + 1) // 2):
                z = [1 if s in chosen else -1 for s in range(p)]
                reconstructed = (
                    q * kd
                    + sum(counts[s] * z[s] for s in range(p))
                    + eps
                    * sum(
                        K[s, t] * z[s] * z[t]
                        for s, t in itertools.combinations(range(p), 2)
                    )
                )
                assert reconstructed == 4 + z[0]


def test_p7_rigid_star_enumeration_and_orbits():
    for populated_type in (-1, 1):
        result = classify(populated_type)
        assert result["generated_candidates"] == 238644
        assert result["survivor_count"] == 2250
        assert result["stabilizer_size"] == 48
        assert result["orbit_count"] == 56
        assert sum(orbit["size"] for orbit in result["orbits"]) == 2250
        assert all(orbit["contains_zero"] for orbit in result["orbits"])


def test_k0_zero_capacity_and_p7_coverage():
    capacity = k0_zero_type_capacity()
    assert capacity["p11"]["excluded"] is True
    assert capacity["p13"]["excluded"] is True
    assert capacity["p7_dichotomy"]["no_zero_direction"] == [1] * 8
    coverage = p7_finite_coverage()
    assert coverage["rigid_type_split"]["infeasible_fixed_star_orbits"] == 112
    assert coverage["all_kd_one"]["infeasible"] == 3


def test_theorem_closes_only_the_infinity_point_boundary():
    theorem = theorem_positive_two_point_all_primes()
    assert theorem["proved"] is True
    assert theorem["closes_infinity_plus_point_boundary_all_primes"] is True
    assert theorem["closes_other_boundary_shapes"] is False
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert len(CERTIFICATE_ARCHIVE_SHA256) == 64
    assert len(P7_ORBIT_CERTIFICATE_SHA256) == 64
