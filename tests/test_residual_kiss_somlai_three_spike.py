import numpy as np
import pytest

from residual_kiss_somlai_three_spike import (
    BASE_SPECIAL_DIRECTIONS,
    find_square_triangle_linear_map,
    fourier_direction_ledger,
    kiss_somlai_direction_audit,
    representative_audit,
    signed_pgl_triangle_transitivity,
    triangular_augmented_function,
    triangular_three_spike_certificate,
)


@pytest.mark.parametrize("p", [5, 7, 11, 13])
def test_kiss_somlai_triangle_has_exactly_three_special_directions(p):
    row = kiss_somlai_direction_audit(p)
    assert row["proved_for_this_prime"]
    assert {tuple(direction) for direction in row["special_spatial_directions"]} == set(
        BASE_SPECIAL_DIRECTIONS
    )
    assert row["nonspecial_line_sum"] == (p - 1) // 2


@pytest.mark.parametrize("p", [5, 7, 11, 13])
def test_augmented_triangle_has_exactly_two_double_points(p):
    f0 = triangular_augmented_function(p)
    values, counts = np.unique(f0, return_counts=True)
    assert dict(zip(values.tolist(), counts.tolist())) == {
        0: (p * p - 3 * p + 4) // 2,
        1: p * (p + 3) // 2 - 4,
        2: 2,
    }
    assert f0[(p - 2) * p] == 2
    assert f0[(p - 2) * p + (p - 1)] == 2


@pytest.mark.parametrize("p", [5, 7, 11, 13])
def test_deterministic_map_and_fourier_sign_chain(p):
    image_a, image_b = find_square_triangle_linear_map(p)
    assert image_a != image_b
    ledger = fourier_direction_ledger(
        p,
        (
            image_a,
            image_b,
            (image_a % p + image_b % p) % p
            + ((image_a // p + image_b // p) % p) * p,
        ),
    )
    assert ledger["proved_for_this_prime"]
    assert all(row["spatial_direction_character"] == 1 for row in ledger["rows"])
    assert all(row["finite_Paley_multiplier_sign"] == 1 for row in ledger["rows"])
    if p % 4 == 1:
        assert ledger["gauss_sum_over_Fp2"] == "-1*p"
        assert all(row["annihilator_direction_character"] == -1 for row in ledger["rows"])
    else:
        assert ledger["gauss_sum_over_Fp2"] == "1*p"
        assert all(row["annihilator_direction_character"] == 1 for row in ledger["rows"])


@pytest.mark.parametrize("p", [5, 7, 11])
def test_positive_signed_triples_form_one_signed_psl_orbit(p):
    orbit = signed_pgl_triangle_transitivity(p)
    assert orbit["switching_factor_exact_on_all_edges"]
    assert orbit["ordered_positive_triangles"] == orbit["PSL_2_q_order"]
    assert orbit["positive_unordered_triangles"] == orbit["negative_unordered_triangles"]
    assert orbit["positive_support_triangles_single_PSL_orbit"]
    assert orbit["positive_signed_triples_single_signed_PSL_orbit"]
    assert orbit["triangular_completion_exists_for_every_signed_triple_datum"]
    assert not orbit["all_completions_of_fixed_datum_classified"]


@pytest.mark.parametrize("p", [5, 7, 11, 13])
def test_integral_eigenvector_has_positive_all_bad_three_spike_endpoint(p):
    row = triangular_three_spike_certificate(p)
    assert row["proved_for_this_prime"]
    assert row["Qf_equals_pf_minus_constant"]
    assert row["C_y_equals_p_y"]
    assert row["spike_edge_values"] == [1, 1, 1]
    assert row["shadow_is_boolean"]
    assert row["signed_triple_values_on_support"] == [-1, -1, -1]
    assert row["all_three_shell_signs_bad"]
    assert row["boolean_shadow_defect"] == 6 * p - 12
    assert len(row["spike_indices"]) == 3
    assert not row["residual_ii_closed"]


def test_representative_audit_keeps_the_global_gate_open():
    out = representative_audit()
    assert out["all_checks"]
    assert set(out["rows"]) == {"5", "7", "11", "13"}
    assert not out["residual_ii_closed"]


@pytest.mark.parametrize("p", [2, 3, 4, 9])
def test_construction_rejects_out_of_scope_parameters(p):
    with pytest.raises(ValueError, match="odd prime p>=5"):
        triangular_three_spike_certificate(p)
