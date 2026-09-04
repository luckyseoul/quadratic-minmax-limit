import pytest

from e1_gmin_m4_p31_row_atom_incremental_bound import AtomRowSpec
from e1_gmin_m4_p31_semimetric_row_separator import (
    EDGES,
    generated_atom_row_validation,
    learned_cut_score,
    normalized_cut_metric,
    separate_semimetric_cone,
    theorem_record,
    three_label_cut_strict_witness,
)


def test_generated_hard_and_opposite_atom_rows_satisfy_real_oracle() -> None:
    row = generated_atom_row_validation()
    assert row["proved"] is True
    assert abs(row["hard_margin"]) <= 1e-8
    assert abs(row["opposite_margin"]) <= 1e-8
    assert row["hard_generated_from_eleven_identical_compact_atoms"] is True
    assert row["opposite_generated_from_six_triangles_and_six_compact_atoms"] is True


def test_full_semimetric_lp_strictly_advances_the_496_cut_bank() -> None:
    row = three_label_cut_strict_witness()
    assert row["proved"] is True
    assert row["singleton_two_label_bank_cost"] == 0
    assert row["singleton_two_label_edit_lower_bound"] == 0
    assert row["violated_cut_subset"] == [0, 1, 2]
    assert row["cut_numerator"] == 6
    assert row["cut_metric_normalization_denominator"] == 84
    assert row["exact_normalized_cut_margin"] == "1/14"
    assert row["lp_separation_margin"] == pytest.approx(1 / 14, abs=1e-10)
    assert row["lp_sparse_weight_count"] == 84
    assert row["strictly_stronger_than_496_cut_bank"] is True


def test_returned_lp_weights_are_a_normalized_semimetric_gpu_cut() -> None:
    coefficients = three_label_cut_strict_witness()["coefficients"]
    result = separate_semimetric_cone(coefficients, AtomRowSpec.hard(11))
    diagnostics = result["semimetric_diagnostics"]
    assert result["proved_numeric"] is True
    assert result["edge_sum_matches_atom_count"] is True
    assert result["strictly_separates_real_atom_relaxation"] is True
    assert len(result["normalized_semimetric_weights"]) == len(EDGES) == 465
    assert diagnostics["normalization"] == pytest.approx(1.0, abs=2e-7)
    assert diagnostics["minimum_weight"] >= -1e-10
    assert diagnostics["maximum_triangle_excess"] <= 2e-7
    assert result["separation_margin"] == pytest.approx(1 / 14, abs=1e-10)
    assert result["integral_decomposition_certified"] is False


def test_learned_cut_cell_delta_is_one_multiply_add() -> None:
    coefficients = three_label_cut_strict_witness()["coefficients"].copy()
    weights = normalized_cut_metric((0, 1, 2))
    old = learned_cut_score(coefficients, weights, rhs=0.0)

    crossing = (0, 9)
    delta = -3
    coefficients[crossing] = coefficients.get(crossing, 0) + delta
    new = learned_cut_score(coefficients, weights, rhs=0.0)
    assert new == pytest.approx(
        max(0.0, old + delta * weights[crossing]), abs=1e-12
    )

    internal = (0, 1)
    old = learned_cut_score(coefficients, weights, rhs=0.0)
    coefficients[internal] += 7
    new = learned_cut_score(coefficients, weights, rhs=0.0)
    assert weights[internal] == 0.0
    assert new == pytest.approx(old, abs=1e-12)


def test_invalid_oracle_requests_are_rejected() -> None:
    with pytest.raises(ValueError):
        normalized_cut_metric(())
    with pytest.raises(ValueError):
        normalized_cut_metric(tuple(range(31)))
    with pytest.raises(ValueError):
        separate_semimetric_cone({}, AtomRowSpec(1, 1))


def test_theorem_record_does_not_promote_real_to_integral_sufficiency() -> None:
    row = theorem_record()
    assert row["proved"] is True
    assert row["integral_normality_asserted"] is False
    assert row["learned_gpu_cut_cell_update"] == "lhs += delta*d_uv"
    assert row["residual_ii_closed"] is False
