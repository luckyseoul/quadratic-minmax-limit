import pytest

from e1_gmin_m4_symmetric_fixed_parallel_flow import (
    alternating_cycle_steering_theorem,
    clean_block_mixed_moment_theorem,
    even_channel_midpoint_value,
    exact_p7_flow_replay,
    fixed_parallel_flow_theorem,
    four_cycle_monomial_delta,
    theorem_record,
)


def test_fixed_parallel_flow_dimensions_and_scope():
    out = fixed_parallel_flow_theorem(31)
    assert out["projective_directions"] == 32
    assert out["fixed_antipodal_cell_nodes"] == 480
    assert out["projective_midpoint_difference_blocks"] == 1024
    assert out["variables_per_projective_block"] == 225
    assert out["nonfixed_orbit_variables"] == 230400
    assert out["network_matrix"] is True
    assert out["totally_unimodular"] is True
    assert out["arbitrary_used_column_deletion_preserves_TU"] is True
    assert out["unit_box_integrality"] is True
    assert out["nonfixed_transverse_cells_solved"] is False
    assert out["residual_ii_closed"] is False
    assert out["proved"] is True


def test_alternating_cycle_scope_is_transverse_steering_only():
    out = alternating_cycle_steering_theorem(31)
    assert out["radial_classes_per_projective_direction"] == 15
    assert "circulation" in out["flow_fibre_difference"]
    assert "branch-C Paley sign" in out["radial_cycle_sign_hypothesis"]
    assert any(
        "fixed-edge word" in item
        for item in out["parallel_arc_exchange_preserves"]
    )
    assert "four-cycles alone need not connect" in out["punctured_binary_warning"]
    assert "no completion-side cycle" in out["top_endpoint_warning"]
    assert out["proved"] is True


def test_even_channel_midpoint_formula_against_endpoint_definition():
    p = 31
    x = 7
    y = 11
    left = (x - y) % p
    right = (x + y) % p
    alpha = x * x % p
    beta = y * y % p
    for n in range(1, 5):
        degree = 2 * n
        for k in range(n):
            direct = (
                pow(left - right, 2, p)
                * pow(left * right, k, p)
                * pow(left + right, degree - 2 - 2 * k, p)
            ) % p
            assert even_channel_midpoint_value(p, n, k, alpha, beta) == direct


def test_four_cycle_mixed_monomial_factorization():
    for r in range(5):
        for s in range(5):
            out = four_cycle_monomial_delta(31, r, s, 1, 4, 9, 16)
            assert out["identity_holds"] is True
            if r == 0 or s == 0:
                assert out["direct_delta"] == 0


def test_clean_blocks_span_degree_six_and_eight_mixed_forms():
    out = clean_block_mixed_moment_theorem(31)
    assert out["projective_blocks"] == 1024
    assert out["used_orbit_cap"] == 480
    assert out["clean_projective_blocks_at_least"] == 544
    assert out["degree_six_nonzero_biform_zero_bound"] == 184
    assert out["degree_six_spanning_margin"] == 360
    assert out["degree_eight_nonzero_biform_zero_bound"] == 244
    assert out["degree_eight_spanning_margin"] == 300
    assert out["mixed_degree_six_global_forms_spanned"] is True
    assert out["mixed_degree_eight_global_forms_spanned"] is True
    assert out["pure_radial_margins_spanned"] is False
    assert out["conformal_Boolean_sequence_proved"] is False
    assert out["full_transverse_cells_solved"] is False
    assert out["proved"] is True


def test_exact_p7_signed_incidence_and_four_cycle_replay():
    out = exact_p7_flow_replay()
    assert out["nonfixed_columns_replayed"] == 576
    assert out["every_divided_fixed_column_is_a_signed_network_arc"] is True
    assert out["radial_half_edge_columns"] == 72
    assert out["two_endpoint_arc_columns"] == 504
    assert out["parallel_arc_bins"] == 192
    assert out["parallel_arcs_per_bin"] == 3
    assert out["all_bins_have_expected_parallel_multiplicity"] is True
    assert out["four_cycle_replay"]["fixed_projection_zero"] is True
    assert out["four_cycle_replay"]["transverse_projection_nonzero"] is True
    assert out["mixed_monomial_factorization_replay"]["identity_holds"] is True
    assert "not a prime census" in out["role"]
    assert out["proved"] is True


def test_record_keeps_residual_open():
    out = theorem_record(31)
    assert out["proved"] is True
    assert out["common_simple_graph_constructed"] is False
    assert out["nonfixed_transverse_cells_solved"] is False
    assert out["residual_ii_closed"] is False


@pytest.mark.parametrize("p", [2, 9, 15])
def test_rejects_non_odd_primes(p):
    with pytest.raises(ValueError):
        fixed_parallel_flow_theorem(p)


def test_clean_block_degree_eight_threshold_is_explicit():
    with pytest.raises(ValueError):
        clean_block_mixed_moment_theorem(11)
