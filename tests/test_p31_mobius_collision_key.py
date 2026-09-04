from itertools import combinations

import pytest

from e1_gmin_m4_inversion_antisymmetric_radon import (
    _negative_edge,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import mobius_parameter_edges
from e1_gmin_m4_p31_direct_mobius_parallel_design import (
    _oriented_orbit_coefficients,
)
from e1_gmin_m4_p31_mobius_collision_key import (
    CollisionKey,
    DIRECTIONS,
    EXPECTED_SHARED_KEYS,
    FROZEN_LABELLED_KEY_CATALOG_SHA256,
    FROZEN_CORRECTION_SUPPORT,
    FROZEN_HALF_CHOICES,
    P,
    REQUIRED_CANCELLATION_DIRECTION_INDEX,
    centred_edges_share_orbit,
    collision_homothety_scalar,
    frozen_witness_collision_key_replay,
    mobius_half_collision_keys,
    nonorigin_collision_key,
    prescribed_endpoint_candidate,
    target_options_for_collision_key,
    theorem_record,
)


def _half(row: tuple[int, int, int]):
    target_index, auxiliary_index, scale = row
    target = DIRECTIONS[target_index]
    auxiliary_row = DIRECTIONS[auxiliary_index]
    auxiliary = (
        scale * auxiliary_row[0] % P,
        scale * auxiliary_row[1] % P,
    )
    return target, auxiliary


def test_formula_matches_every_nonorigin_edge_in_frozen_family() -> None:
    for row in FROZEN_HALF_CHOICES:
        target, auxiliary = _half(row)
        keys = mobius_half_collision_keys(P, target, auxiliary)
        assert len(keys) == P - 2
        assert set(keys) == set(range(1, P - 1))


def test_equal_key_is_necessary_and_sufficient_for_free_centres() -> None:
    replay = frozen_witness_collision_key_replay()
    for record in replay["shared_collision_keys"]:
        (first_half, first_parameter), (second_half, second_parameter) = record[
            "owners_half_index_parameter"
        ]
        first_target, first_auxiliary = _half(FROZEN_HALF_CHOICES[first_half])
        second_target, second_auxiliary = _half(FROZEN_HALF_CHOICES[second_half])
        scale = collision_homothety_scalar(
            P,
            first_target,
            first_auxiliary,
            first_parameter,
            second_target,
            second_auxiliary,
            second_parameter,
        )
        assert scale not in (None, 0)

        shared_count = 0
        cancelling_count = 0
        for second_center in range(1, P):
            for first_center in range(1, P):
                relation = centred_edges_share_orbit(
                    P,
                    first_target,
                    first_auxiliary,
                    first_parameter,
                    first_center,
                    second_target,
                    second_auxiliary,
                    second_parameter,
                    second_center,
                )
                shared_count += int(relation["same_inversion_orbit"])
                cancelling_count += int(relation["cancellation_orientation"])
        assert shared_count == 2 * (P - 1)
        assert cancelling_count == P - 1

        cancelling_first_center = -scale % P
        relation = centred_edges_share_orbit(
            P,
            first_target,
            first_auxiliary,
            first_parameter,
            cancelling_first_center,
            second_target,
            second_auxiliary,
            second_parameter,
            1,
        )
        assert relation["opposite_physical_edges"] is True
        first_edge = mobius_parameter_edges(
            P, first_target, first_auxiliary, cancelling_first_center
        )[first_parameter]
        second_edge = mobius_parameter_edges(
            P, second_target, second_auxiliary, 1
        )[second_parameter]
        assert first_edge == _negative_edge(P, second_edge)
        first_coefficients = _oriented_orbit_coefficients(
            first_target, first_auxiliary, cancelling_first_center
        )
        second_coefficients = _oriented_orbit_coefficients(
            second_target, second_auxiliary, 1
        )
        orbit = min(first_edge, second_edge)
        assert first_coefficients[orbit] == -second_coefficients[orbit]


def test_unequal_key_has_no_centre_solution() -> None:
    first_target, first_auxiliary = _half(FROZEN_HALF_CHOICES[0])
    second_target, second_auxiliary = _half(FROZEN_HALF_CHOICES[1])
    first_parameter = 1
    first_key = nonorigin_collision_key(
        P, first_target, first_auxiliary, first_parameter
    )
    second_parameter = next(
        parameter
        for parameter in range(1, P - 1)
        if nonorigin_collision_key(
            P, second_target, second_auxiliary, parameter
        )
        != first_key
    )
    assert (
        collision_homothety_scalar(
            P,
            first_target,
            first_auxiliary,
            first_parameter,
            second_target,
            second_auxiliary,
            second_parameter,
        )
        is None
    )
    for first_center, second_center in combinations(range(1, 6), 2):
        relation = centred_edges_share_orbit(
            P,
            first_target,
            first_auxiliary,
            first_parameter,
            first_center,
            second_target,
            second_auxiliary,
            second_parameter,
            second_center,
        )
        assert relation["same_inversion_orbit"] is False


def test_prescribed_endpoint_formula_pins_unique_parameter_and_bucket() -> None:
    candidates = []
    for half_index, row in enumerate(FROZEN_HALF_CHOICES):
        target, auxiliary = _half(row)
        candidate = prescribed_endpoint_candidate(
            P, target, auxiliary, FROZEN_CORRECTION_SUPPORT
        )
        if candidate is not None:
            candidates.append((half_index, *candidate))
    assert len(candidates) == 1
    half_index, parameter, key = candidates[0]
    assert FROZEN_HALF_CHOICES[half_index] == (29, 31, 27)
    assert parameter == 12
    assert key.as_tuple() == (8, 2, 23)
    assert key.spatial_direction != REQUIRED_CANCELLATION_DIRECTION_INDEX


def test_inverse_key_formula_leaves_at_most_two_options_per_target() -> None:
    # Check every target against the five physically shared frozen classes.
    for raw_key in EXPECTED_SHARED_KEYS:
        key = CollisionKey(raw_key[0], raw_key[1:])
        for target_index in range(P + 1):
            options = target_options_for_collision_key(P, target_index, key)
            assert len(options) <= 2
            for auxiliary_index, scale, parameter, _fixed_endpoint in options:
                auxiliary_row = DIRECTIONS[auxiliary_index]
                auxiliary = (
                    scale * auxiliary_row[0] % P,
                    scale * auxiliary_row[1] % P,
                )
                assert (
                    nonorigin_collision_key(
                        P, DIRECTIONS[target_index], auxiliary, parameter
                    )
                    == key
                )


def test_frozen_class_join_replaces_the_108000_centre_scan() -> None:
    result = frozen_witness_collision_key_replay()
    assert result["proved"] is True
    assert result["shared_collision_key_count"] == 5
    assert (
        result["labelled_key_catalog_sha256"]
        == FROZEN_LABELLED_KEY_CATALOG_SHA256
    )
    assert tuple(
        record["key"] for record in result["shared_collision_keys"]
    ) == EXPECTED_SHARED_KEYS
    assert result["distinct_half_pair_key_count"] == 5
    assert result["derived_shared_orbit_center_incidences"] == 300
    assert result["derived_cancelling_center_incidences"] == 150
    assert result["prescribed_endpoint_candidate_count"] == 1
    assert result["required_direction_prescribed_collision_exists"] is False
    assert result["frozen_sixteen_half_family_physically_excluded"] is True
    assert result["centre_enumeration_needed"] is False


def test_theorem_scope_does_not_claim_global_closure() -> None:
    result = theorem_record()
    assert result["proved"] is True
    assert result["equal_keys_iff_unit_edges_are_homothetic"] is True
    assert result["same_orbit_center_condition"] == "c1=+/-lambda*c2"
    assert result["cancellation_center_condition"] == "c1=-lambda*c2"
    assert result["residual_ii_closed"] is False


def test_degenerate_parameters_and_endpoint_requests_are_rejected() -> None:
    target, auxiliary = _half(FROZEN_HALF_CHOICES[0])
    with pytest.raises(ValueError):
        nonorigin_collision_key(P, target, auxiliary, 0)
    with pytest.raises(ValueError):
        nonorigin_collision_key(P, target, auxiliary, -1)
    with pytest.raises(ValueError):
        prescribed_endpoint_candidate(P, target, auxiliary, (2, 2))
    with pytest.raises(ValueError):
        prescribed_endpoint_candidate(P, target, auxiliary, (2, P + 1))
    with pytest.raises(ValueError):
        CollisionKey(2, (2, 23))


def test_direction_indexing_agrees_with_repository_order() -> None:
    assert DIRECTIONS == tuple(projective_functionals(P))
