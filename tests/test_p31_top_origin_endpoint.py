import pytest

from e1_gmin_m4_p31_top_origin_endpoint import (
    origin_line_injectivity_certificate,
    p31_top_origin_endpoint,
    theorem_record,
)


def test_equal_origin_orbits_force_equal_auxiliary_directions_and_sign() -> None:
    result = origin_line_injectivity_certificate(31)
    assert result["proved"] is True
    assert result["projective_auxiliary_directions"] == 32
    assert result["distinct_projective_kernel_directions"] == 32
    assert result["auxiliary_to_kernel_direction_bijection"] is True
    assert result["equal_origin_orbits_force_equal_auxiliary_direction"] is True
    assert result["origin_edge_sign_equals_auxiliary_sign"] is True


@pytest.mark.parametrize(
    ("fixed_type", "sign_counts", "raw_signed", "opposite_minus_hard"),
    (
        ("hard", {"hard": 13, "opposite": 3}, 10, -10),
        ("opposite", {"hard": 14, "opposite": 2}, 12, -12),
    ),
)
def test_top_endpoint_forces_distinct_auxiliaries_and_nonorigin_cancellation(
    fixed_type: str,
    sign_counts: dict[str, int],
    raw_signed: int,
    opposite_minus_hard: int,
) -> None:
    result = p31_top_origin_endpoint(fixed_type)
    assert result["proved"] is True
    assert result["v_weight"] == 17
    assert result["fixed_direction_lies_in_v_support"] is True
    assert result["auxiliary_parity_weight"] == 16
    assert result["auxiliary_direction_occurrences"] == 16
    assert result["auxiliary_directions_are_distinct"] is True
    assert result["cancellation_units"] == 1
    assert result["origin_cancellation_excluded"] is True
    assert result["sole_cancellation_is_nonorigin"] is True
    assert result["fixed_antipodal_edge_avoids_origin"] is True
    assert result["origin_unsigned_degree"] == 16
    assert result["auxiliary_sign_counts"] == sign_counts
    assert result["origin_raw_signed_degree_hard_minus_opposite"] == raw_signed
    assert result["origin_signed_degree_opposite_minus_hard"] == opposite_minus_hard


def test_theorem_record_keeps_scope_open() -> None:
    result = theorem_record()
    assert result["proved"] is True
    assert result["origin_unsigned_degree"] == 16
    assert result["origin_raw_signed_degree_cases"] == [10, 12]
    assert result["origin_cancellation_possible"] is False
    assert result["residual_ii_closed"] is False


def test_invalid_fixed_edge_type_is_rejected() -> None:
    with pytest.raises(ValueError):
        p31_top_origin_endpoint("square")
