from e1_gmin_m4_prop15692 import (
    affine_binary_radon_isomorphism,
    p19_inverse_radon_profile_reduction,
)


def test_binary_affine_radon_dimensions_and_inverse():
    for p in (3, 5, 19):
        row = affine_binary_radon_isomorphism(p)
        assert row["incidence_gram_over_F2"] == "A^T A = I + J"
        assert row["source_dimension"] == p * p - 1
        assert row["target_dimension"] == p * p - 1
        assert row["inverse"] == "x = A^T r"
        assert row["proved"] is True


def test_p19_fourteen_profiles_all_reduce_to_inverse_weight():
    row = p19_inverse_radon_profile_reduction()
    assert row["profile_count"] == 14
    assert row["additional_linear_compatibility_conditions"] == 0
    assert row["all_profiles_pass_inverse_weight_mod_four"] is True
    assert {item["pair_slack"] for item in row["rows"]} == {16, 20, 24, 28, 32}


def test_p19_second_moment_relaxation_has_even_witnesses():
    row = p19_inverse_radon_profile_reduction()
    for item in row["rows"]:
        witness = item["all_even_second_moment_witness"]
        assert set(witness) == {4, 6, 8}
        assert all(value >= 0 for value in witness.values())
        assert sum(witness.values()) == 1
        assert item["second_moments_force_positive_odd_density"] is False


def test_prop15692_changes_no_top_level_gate():
    row = p19_inverse_radon_profile_reduction()
    assert row["p19_second_all_finite_endpoint_closed"] is False
    assert row["closes_residual_ii"] is False
    assert row["closes_R1"] is False
    assert row["closes_type_I"] is False
    assert row["L_status"] == "OPEN"
