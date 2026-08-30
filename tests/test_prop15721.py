"""Prop. 15.721: signed PSL transport removes the all-finite ladder."""
from __future__ import annotations

from e1_gmin_m4_prop15675 import first_even_survivor
from e1_gmin_m4_prop15721 import (
    minus_one_is_square_in_prime_square,
    mobius_boundary_normalization,
    old_all_finite_ladder_coverage,
    signed_relative_flip_transport,
    theorem_boundary_transport_floor,
    transported_boundary_exclusion,
    universal_boundary_transport_certificate,
)


def test_explicit_mobius_normalization_lies_in_psl():
    for p in (5, 7, 11, 17, 19, 101):
        assert minus_one_is_square_in_prime_square(p) is True
        row = mobius_boundary_normalization(p)
        assert row["determinant"] == -1
        assert row["lies_in_PSL_2_q"] is True
        assert row["selected_boundary_point_maps_to_infinity"] is True


def test_signed_conjugation_only_permutes_the_relative_flip_mask():
    row = signed_relative_flip_transport()
    assert row["proved"] is True
    assert row["diagonal_signs_cancel_from_relative_mask"] is True
    assert row["flip_set_size_preserved"] is True
    assert row["odd_degree_boundary_is_permuted"] is True
    assert row["both_separation_inequalities_preserved"] is True


def test_every_even_boundary_through_p_minus_one_is_excluded():
    for p in (17, 19, 23, 29, 31, 41):
        for d in range(0, p, 2):
            row = transported_boundary_exclusion(p, d)
            assert row["excluded"] is True, (p, d, row["method"])
        first_open = transported_boundary_exclusion(p, p + 1)
        assert first_open["excluded"] is False
        assert "strict pair deficit" in first_open["remaining_at_first_open_size"]


def test_old_first_survivor_is_transport_excluded():
    for p in (19, 23, 29, 31, 37, 41, 43, 101):
        old = first_even_survivor(p)
        assert transported_boundary_exclusion(p, old)["excluded"] is True
    assert transported_boundary_exclusion(17, 14)["excluded"] is True


def test_old_second_shells_are_transport_excluded_too():
    rows = old_all_finite_ladder_coverage()
    assert rows["proved"] is True
    for row in rows["rows"].values():
        assert row["first_transport_excluded"] is True
        assert row["second_transport_excluded"] is True


def test_theorem_stops_honestly_at_total_boundary_p_plus_one():
    out = theorem_boundary_transport_floor()
    assert out["proved"] is True
    assert out["universal_certificate"]["proved"] is True
    assert out["first_boundary_size_not_excluded"] == "p+1"
    assert out["remaining"]["residual_ii"] is False
    assert out["remaining"]["type_I"] is False
    assert out["remaining"]["limit_exists"] is False


def test_universal_scope_uses_the_exact_even_size_partition():
    row = universal_boundary_transport_certificate()
    assert row["scope"] == "every prime p>=17"
    assert row["partition_disjoint_and_exhaustive"] is True
    assert row["middle_range_matches_prop_15_669"] is True
    assert row["endpoint_p_minus_1_maps_to_infinity_plus_p_minus_2"] is True
    assert row["proved"] is True
