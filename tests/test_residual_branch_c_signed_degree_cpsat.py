import pytest

from e1_gmin_m4_mobius_half_symmetric import (
    mobius_parameter_edges,
    paley_direction_sign,
    paley_edge_sign,
)
from scripts.residual_branch_c_signed_degree_cpsat import (
    build_geometry,
    build_model,
    geometry_identity_self_test,
    mobius_top_auxiliary_indicator,
    mobius_top_origin_profile,
    realise_six_triangle_occurrences,
    target_row_sum_certificate,
)


EXPECTED_MODEL_HASHES = {
    "any": "cea8f8dc2c4cab04aa12bc8dd9752cdd533e572cb38c91f1e4dce88257534c55",
    "hard": "7b26cef3cb83a225fcdb343e8ecbfcb7683eefd60c7f805a778fb4afd752365c",
    "opposite": "2d5381ea2c4d3e459faa13e032291046390046441173d56b0baa9ca3ea3465b9",
}

EXPECTED_MOBIUS_TOP_MODEL_HASHES = {
    "any": "048a5c1d3baf369c03ff29712cc2cf773566c3af5be2809c1cc19363375fe880",
    "hard": "d08289a60a60846e550a6ab06e42a310dc82453c503a17848e6d67eed430b79f",
    "opposite": "4ab7ae87fec5031c2c9125cc65dc157bcb194168a44a54f8d3ff1c052c160bb8",
}


def test_actual_edge_geometry_replays_paley_signs_and_line_inverse() -> None:
    geometry = build_geometry()
    replay = geometry_identity_self_test(geometry)
    assert len(geometry.points) == 961
    assert len(geometry.directions) == 32
    assert geometry.epsilon.count(1) == geometry.epsilon.count(-1) == 16
    assert len(geometry.antipodal_representatives) == 480
    assert replay["proved"] is True
    assert replay["sample_edge_count"] == 37
    assert replay["raw_signed_degree_sum"] == replay["twice_edge_sign_sum"]


@pytest.mark.parametrize("fixed_edge_sign", ("any", "hard", "opposite"))
def test_full_and_fixed_edge_sign_shard_models_validate_and_hash(
    fixed_edge_sign: str,
) -> None:
    model, _, metadata = build_model(build_geometry(), fixed_edge_sign)
    assert model.validate() == ""
    assert metadata["cp_sat_variable_count"] == 7395
    assert metadata["cp_sat_constraint_count"] == 5513 + int(fixed_edge_sign != "any")
    assert metadata["cp_sat_textproto_sha256"] == EXPECTED_MODEL_HASHES[fixed_edge_sign]
    assert metadata["fixed_edge_sign_shard"] == fixed_edge_sign
    assert metadata["hard_center_domain"] == "F_31"
    assert metadata["selected_edge_count"] == 479
    assert metadata["raw_signed_degree_total"] == -54


def test_two_fixed_edge_sign_shards_are_an_exact_partition() -> None:
    geometry = build_geometry()
    by_sign = {1: 0, -1: 0}
    for representative in geometry.antipodal_representatives:
        direction = geometry.spatial_direction[representative]
        by_sign[geometry.epsilon[direction]] += 1
    assert by_sign == {1: 240, -1: 240}


@pytest.mark.parametrize("fixed_edge_sign", ("any", "hard", "opposite"))
def test_mobius_top_origin_slice_validates_and_hashes(fixed_edge_sign: str) -> None:
    model, _, metadata = build_model(build_geometry(), fixed_edge_sign, True)
    assert model.validate() == ""
    assert metadata["mobius_top_origin_slice"] is True
    assert metadata["hard_center_domain"] == "F_31^*"
    assert metadata["cp_sat_variable_count"] == 7427
    assert metadata["cp_sat_constraint_count"] == 5596 + int(
        fixed_edge_sign != "any"
    )
    assert (
        metadata["cp_sat_textproto_sha256"]
        == EXPECTED_MOBIUS_TOP_MODEL_HASHES[fixed_edge_sign]
    )


def test_mobius_t0_edge_sign_is_exactly_auxiliary_direction_sign() -> None:
    geometry = build_geometry()
    hard_target = next(
        direction
        for direction, epsilon in zip(geometry.directions, geometry.epsilon)
        if epsilon == 1
    )
    for auxiliary in geometry.directions:
        if auxiliary == hard_target:
            continue
        t0_edge = mobius_parameter_edges(31, hard_target, auxiliary, 1)[0]
        assert paley_edge_sign(31, t0_edge) == paley_direction_sign(31, auxiliary)


def test_mobius_top_origin_profiles_use_the_verified_sign_convention() -> None:
    assert mobius_top_origin_profile("hard") == {
        "positive_degree": 13,
        "negative_degree": 3,
        "signed_degree": 10,
        "ordinary_degree": 16,
    }
    assert mobius_top_origin_profile("opposite") == {
        "positive_degree": 14,
        "negative_degree": 2,
        "signed_degree": 12,
        "ordinary_degree": 16,
    }
    with pytest.raises(ValueError):
        mobius_top_origin_profile("any")


@pytest.mark.parametrize("fixed_sign", (1, -1))
def test_direction_resolved_mobius_auxiliary_word(fixed_sign: int) -> None:
    geometry = build_geometry()
    hard = [index for index, eps in enumerate(geometry.epsilon) if eps == 1]
    opposite = [index for index, eps in enumerate(geometry.epsilon) if eps == -1]
    e_values = {
        direction: 12 if position < 2 else 11
        for position, direction in enumerate(hard)
    }
    q_values = {
        direction: 16 if position < 13 else 15
        for position, direction in enumerate(opposite)
    }
    v_support = [
        direction
        for direction in range(32)
        if (
            (12 - e_values[direction])
            if direction in e_values
            else (16 - q_values[direction])
        )
    ]
    fixed_direction = next(
        direction
        for direction in v_support
        if geometry.epsilon[direction] == fixed_sign
    )
    z = mobius_top_auxiliary_indicator(
        geometry, e_values, q_values, fixed_direction
    )
    assert sum(z) == 16
    assert z[fixed_direction] == 0
    assert all(
        value
        == (
            (12 - e_values[direction])
            if direction in e_values
            else (16 - q_values[direction])
        )
        - int(direction == fixed_direction)
        for direction, value in enumerate(z)
    )
    positive = sum(z[direction] for direction in hard)
    negative = sum(z[direction] for direction in opposite)
    expected = mobius_top_origin_profile("hard" if fixed_sign == 1 else "opposite")
    assert (positive, negative) == (
        expected["positive_degree"],
        expected["negative_degree"],
    )


def test_hard_and_opposite_target_degree_sums_have_common_line_total() -> None:
    certificate = target_row_sum_certificate()
    assert certificate["proved"] is True
    assert certificate["common_raw_line_sum"] == -54
    assert certificate["hard_rows"] == {
        "11": {
            "transverse_degree_sum": -82,
            "parallel_quota": 14,
            "raw_line_sum": -54,
        },
        "12": {
            "transverse_degree_sum": -84,
            "parallel_quota": 15,
            "raw_line_sum": -54,
        },
    }
    assert certificate["opposite_rows"] == {
        "15": {
            "transverse_degree_sum": 24,
            "parallel_quota": 15,
            "raw_line_sum": -54,
        },
        "16": {
            "transverse_degree_sum": 22,
            "parallel_quota": 16,
            "raw_line_sum": -54,
        },
    }


@pytest.mark.parametrize(
    "counts",
    (
        [6, 6, 6] + [0] * 28,
        [6, 6, 5, 1] + [0] * 27,
        [3, 3, 3, 3, 3, 3] + [0] * 25,
        [1] * 18 + [0] * 13,
    ),
)
def test_opposite_occurrence_vector_has_independent_triangle_replay(
    counts: list[int],
) -> None:
    triangles = realise_six_triangle_occurrences(counts)
    assert len(triangles) == 6
    assert all(len(triangle) == len(set(triangle)) == 3 for triangle in triangles)
    replay = [0] * 31
    for triangle in triangles:
        for label in triangle:
            replay[label] += 1
    assert replay == counts


@pytest.mark.parametrize(
    "counts",
    ([7, 6, 5] + [0] * 28, [6, 6, 5] + [0] * 28),
)
def test_triangle_replay_rejects_nonrealisable_occurrence_vectors(
    counts: list[int],
) -> None:
    with pytest.raises(ValueError):
        realise_six_triangle_occurrences(counts)


def test_model_rejects_nonpartition_fixed_edge_shard() -> None:
    with pytest.raises(ValueError):
        build_model(build_geometry(), "square")
