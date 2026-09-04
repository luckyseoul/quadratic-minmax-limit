import pytest

from e1_gmin_m4_p31_top_nonorigin_collision import (
    opposite_direct_collision,
    opposite_swapped_collision,
    pair_collision_locus_theorem,
    theta_pair_intersection,
    theorem_record,
    top_parallel_collision_ledger,
    triple_overlap_countermechanism,
)


def test_opposite_direct_locus_pins_one_nonorigin_orbit() -> None:
    result = opposite_direct_collision(31, alpha=0, beta=2, q=2)
    assert result["proved"] is True
    assert result["locus_equations"] == ["q*r=1", "alpha+beta=2"]
    assert result["auxiliary_directions_distinct"] is True
    assert result["matching"] == "opposite-direct"
    assert result["shared_inversion_orbits"] == 1
    assert result["common_edge_is_nonorigin"] is True


def test_opposite_swapped_quadratic_has_one_replayed_root() -> None:
    result = opposite_swapped_collision(31, alpha=4, beta=9, z=6)
    assert result["proved"] is True
    assert result["parameters"]["q"] == 2
    assert result["parameters"]["r"] == 3
    assert result["parameters"]["A"] == 6
    assert result["parameters"]["B"] == 12
    assert result["matching"] == "opposite-swapped"
    assert result["at_most_two_roots"] is True
    assert result["shared_inversion_orbits"] == 1
    assert result["common_edge_is_nonorigin"] is True


def test_distinct_auxiliary_pair_theorem_excludes_rigid_double_point() -> None:
    result = pair_collision_locus_theorem()
    assert result["proved"] is True
    assert result["distinct_auxiliary_condition"] == "alpha*beta!=1"
    assert result[
        "distinct_auxiliaries_force_at_most_one_opposite_shared_orbit_per_pair"
    ] is True
    assert result["pair_only_description_of_top_cancellation_is_complete"] is False
    assert "M1=M2" in result["rigid_two_overlap_point_excluded"]


@pytest.mark.parametrize(("theta", "phi"), ((1, 2), (1, 3), (2, 3)))
def test_theta_family_pairs_share_exactly_the_prescribed_edge(
    theta: int, phi: int
) -> None:
    result = theta_pair_intersection(31, theta, phi)
    assert result["proved"] is True
    assert result["common_inversion_orbits"] == 1
    assert result["same_orientation_common_edges"] == 1
    assert result["common_edge"] == [[0, 1], [1, 0]]


def test_explicit_three_half_overlap_is_ternary_with_one_cancellation() -> None:
    result = triple_overlap_countermechanism()
    assert result["proved"] is True
    assert result["thetas"] == [1, 2, 3]
    assert result["centers"] == [1, 1, -1]
    assert result["target_signs"] == [1, 1, 1]
    assert result["projective_auxiliaries"] == [[1, 2], [1, 3], [1, 4]]
    assert result["auxiliary_signs"] == [1, 1, -1]
    assert result["auxiliary_directions_distinct"] is True
    assert result["pairwise_shared_orbit_counts"] == [1, 1, 1]
    assert result["common_orientation_coefficients"] == [-1, -1, 1]
    assert result["common_orientation_multiplicities"] == {
        "negative": 2,
        "positive": 1,
    }
    assert result["raw_orbit_occurrences"] == 90
    assert result["final_support_orbits"] == 88
    assert result["cancellation_units"] == 1
    assert result["full_three_trade_sum_is_ternary"] is True
    assert result["common_edge_is_nonorigin"] is True
    assert result["cancellation_direction"] == [1, 1]
    assert result["cancellation_direction_type"] == "hard"
    assert result["full_p31_top_lift_constructed"] is False


@pytest.mark.parametrize(
    ("fixed_type", "cancellation_type", "trace"),
    (
        ("hard", "hard", -26),
        ("hard", "opposite", -30),
        ("opposite", "hard", -24),
        ("opposite", "opposite", -28),
    ),
)
def test_pair_and_triple_have_the_same_parallel_pullback(
    fixed_type: str, cancellation_type: str, trace: int
) -> None:
    result = top_parallel_collision_ledger(fixed_type, cancellation_type)
    assert result["proved"] is True
    assert result["weighted_trace"] == trace
    assert result["identity_applies_to_pair_1_to_1_overlap"] is True
    assert result["identity_applies_to_triple_2_to_1_overlap"] is True


def test_theorem_record_marks_pair_only_models_incomplete_and_gate_open() -> None:
    result = theorem_record()
    assert result["proved"] is True
    assert result["pair_only_top_models_are_complete"] is False
    assert result["full_top_common_graph_constructed"] is False
    assert result["branch_c_closed"] is False
    assert result["residual_ii_closed"] is False


def test_degenerate_pair_and_triple_parameters_are_rejected() -> None:
    with pytest.raises(ValueError):
        opposite_direct_collision(31, alpha=1, beta=1, q=2)
    with pytest.raises(ValueError):
        opposite_swapped_collision(31, alpha=4, beta=9, z=7)
    with pytest.raises(ValueError):
        theta_pair_intersection(31, 1, 1)
