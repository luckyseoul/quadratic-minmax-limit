import pytest

from e1_gmin_m4_inversion_symmetric_lattice import (
    even_moment_cokernel,
    exact_symmetric_mod2_rank,
    mobius_central_remainder_box,
    symmetric_dimensions,
    symmetric_integral_lattice_theorem,
    symmetric_mod2_decomposition,
    theorem_record,
)


def test_symmetric_source_target_and_kernel_ranks():
    for p in (3, 5, 7, 11, 31, 43):
        h = (p - 1) // 2
        d = p + 1
        row = symmetric_dimensions(p)
        assert row["proved"]
        assert row["source_fixed_antipodal_edge_rank"] == d * h
        assert row["target_plus_rank"] == d * h * (h + 1)
        assert row["kernel_plus_rank"] == d * p * h * h


def test_even_moments_are_exactly_the_symmetric_cokernel():
    for p in (3, 5, 7, 11, 31, 43):
        h = (p - 1) // 2
        row = even_moment_cokernel(p)
        assert row["proved"]
        assert row["even_moment_rank_closed"] == (
            (h - 1) * (2 * h * h + 5 * h + 6) // 6
        )
        assert (
            row["even_moment_rank_direct"] + row["odd_moment_rank_direct"]
            == row["full_cokernel_rank"]
        )


def test_symbolic_mod2_fixed_nonfixed_split_is_surjective():
    for p in (3, 5, 7, 11, 31):
        h = (p - 1) // 2
        d = p + 1
        row = symmetric_mod2_decomposition(p)
        assert row["proved"]
        assert row["fixed_edge_map_injective"]
        assert row["fixed_edge_map_rank"] == d * h
        assert row["nonfixed_pair_map_surjective"]
        assert row["paired_nonfixed_target_rank"] == d * h * h
        assert row["symmetric_mod2_image_rank"] == d * h * (h + 1)
        assert row["symmetric_mod2_surjective"]


def test_exact_small_mod2_matrices_match_both_summands():
    for p in (3, 5, 7, 11):
        row = exact_symmetric_mod2_rank(p)
        assert row["proved"]
        assert row["full_symmetric_map_rank"] == (
            row["fixed_map_rank"] + row["nonfixed_pair_map_rank"]
        )
        assert (
            row["full_symmetric_map_rank"]
            == row["expected_symmetric_target_rank"]
        )


def test_integral_cokernel_argument_preserves_boolean_scope():
    row = symmetric_integral_lattice_theorem(31)
    assert row["proved"]
    assert row["cokernel"] == "(Z/31Z)^1239"
    assert row["cokernel_injects_into_full_edge_Radon_cokernel"]
    assert row["image_in_full_cokernel"] == "the +1 inversion eigenspace"
    assert row["unrestricted_signed_integral_central_lift_proved"]
    assert not row["restricted_Boolean_central_lift_proved"]


def test_mobius_remainder_box_freezes_used_orbits():
    row = mobius_central_remainder_box(31, 16)
    assert row["proved"]
    assert row["mobius_used_nonfixed_orbits"] == 480
    assert row["antisymmetric_identity"] == "q_U-Jq_U=z"
    assert row["central_target_formula"] == "T_U=Y-Rq_U=(Y+IY-RC_U)/2"
    assert row["remaining_source_box"] == {
        "used_nonfixed_orbits_after_subtraction": "{0}",
        "unused_nonfixed_orbits": "{0,tau_O*(e+Je)}",
        "fixed_antipodal_edges": "{0,tau_f*f}",
    }
    assert row["unrestricted_integral_lift_if_even_moments_vanish"]
    assert row["mod2_central_lift_always_exists_for_compatible_target"]
    assert not row["same_lift_integral_and_in_box_proved"]
    assert not row["residual_ii_closed"]


def test_theorem_record_keeps_only_the_central_box_open():
    record = theorem_record(31)
    assert record["proved_all"]
    assert record["proved"]["symmetric_cokernel_is_even_moment_eigenspace"]
    assert record["proved"]["symmetric_mod2_map_surjective"]
    assert record["proved"]["unrestricted_signed_integral_central_lift_characterized"]
    assert not record["proved"]["restricted_central_Boolean_box_nonempty"]
    assert not record["proved"]["residual_ii_closed"]
    assert record["L_status"] == "OPEN"


def test_parameter_guards():
    with pytest.raises(ValueError):
        symmetric_dimensions(9)
    with pytest.raises(ValueError):
        exact_symmetric_mod2_rank(13)
    with pytest.raises(ValueError):
        mobius_central_remainder_box(31, 17)
