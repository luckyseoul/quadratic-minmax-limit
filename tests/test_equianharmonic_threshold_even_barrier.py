from fractions import Fraction

from e1_gmin_m4_equianharmonic_threshold_even_barrier import (
    CHANNELS,
    JACOBIAN_DETERMINANT,
    JACOBIAN_FACTORS,
    U_AFFINE_COEFFICIENTS,
    U_AFFINE_CONSTANT,
    equianharmonic_threshold_even_barrier_certificate,
    threshold_excess_assembly_certificate,
    trade_deviation_polynomials,
)


def test_threshold_excess_assembly_is_exhaustive_in_required_range():
    record = threshold_excess_assembly_certificate()
    assert record["proved"] is True
    assert record["global_component_excess"] == "b-2*Delta=1"
    assert record["positive_mass_negative_mass_possibilities"] == [
        [1, 0],
        [2, 1],
        [3, 2],
    ]
    assert record["positive_block_assemblies"] == {
        "1": [["cap"], ["F"]],
        "2": [["HH"], ["F", "cap"]],
        "3": [["HH", "cap"]],
    }
    assert record["zero_excess_tuples_K_AE_cycle_rank_caps"] == [
        (0, 1, 0, 0),
        (2, 0, 0, 1),
        (4, 0, 0, 0),
    ]
    assert len(record["minus_one_tuples_K_AE_cycle_rank_caps"]) == 5
    assert len(record["minus_two_tuples_K_AE_cycle_rank_caps"]) == 6


def test_u_v_polynomials_and_characteristic_zero_barrier():
    record = equianharmonic_threshold_even_barrier_certificate()
    assert record["proved"] is True
    assert record["channels"] == list(CHANNELS)
    assert len(trade_deviation_polynomials("U")) == 7
    assert len(trade_deviation_polynomials("V")) == 7
    assert record["U_odd_edge_orbit_chain_exact"] is True
    assert record["V_odd_edge_orbit_chain_exact"] is True
    assert record["U_affine_syndrome_invariant"] == {
        "coefficients_in_channel_order": [
            str(value) for value in U_AFFINE_COEFFICIENTS
        ],
        "constant_per_trade": str(U_AFFINE_CONSTANT),
        "proved": True,
    }
    assert all(isinstance(value, Fraction) for value in U_AFFINE_COEFFICIENTS)
    jacobian = record["mixed_jacobian"]
    assert jacobian["determinant"] == JACOBIAN_DETERMINANT
    assert jacobian["factorization"] == JACOBIAN_FACTORS
    product = 1
    for prime, exponent in JACOBIAN_FACTORS.items():
        product *= prime**exponent
    assert product == JACOBIAN_DETERMINANT
    assert jacobian["nonzero_characteristic_zero"] is True
    assert jacobian["dominant_outside_displayed_characteristics"] is True


def test_p31_p43_atom_replays_match_symbolic_trade_families():
    replays = equianharmonic_threshold_even_barrier_certificate()[
        "witness_replays"
    ]
    assert replays["p31_U_x18"]["degree_six_deviation"] == [3, 9, 8]
    assert replays["p31_U_x18"]["degree_eight_deviation"] == [20, 14, 21, 4]
    assert replays["p43_U_x38"]["degree_six_deviation"] == [35, 23, 14]
    assert replays["p43_U_x38"]["degree_eight_deviation"] == [36, 12, 8, 25]
    assert replays["p43_V_x7"]["degree_six_deviation"] == [36, 21, 18]
    assert replays["p43_V_x7"]["degree_eight_deviation"] == [8, 14, 20, 40]
    assert all(
        replay["symbolic_matches_atom_witness"]
        for replay in replays.values()
    )


def test_barrier_does_not_claim_the_open_lifts():
    record = equianharmonic_threshold_even_barrier_certificate()
    assert record["finite_field_rational_zero_syndrome_trade_matching_constructed"] is False
    assert record["uniform_zero_degree_six_eight_exclusion_proved"] is False
    assert record["common_global_form_lift_constructed"] is False
    assert record["Boolean_lift_constructed"] is False
    assert record["residual_ii_closed"] is False
