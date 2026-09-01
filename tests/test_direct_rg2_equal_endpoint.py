"""Exact tests for the direct equal-endpoint skew RG2 algebra."""
from __future__ import annotations

from fractions import Fraction

import pytest

from direct_rg2_equal_endpoint import (
    equal_endpoint_block,
    equal_endpoint_k_by_cuts,
    equal_endpoint_k_by_frames,
    equal_endpoint_k_by_pairs,
    exhaustive_local_block_certificate,
    exhaustive_small_order_certificate,
    hereditary_endpoint_certificate,
    n5_cycle_chord_signing,
    n5_disk_counterexample,
    n5_maximizer_basis,
    quadratic_energy,
    simple_walk_absolute_mean,
    skew_norm_floor_certificate,
    skew_norm_floor_exact,
    skew_signings,
    symmetric_signings,
    verify_equal_endpoint_algebra,
    zero_error_disk_holds,
)


def test_local_hadamard_block_formula_is_exhaustive() -> None:
    row = exhaustive_local_block_certificate()
    assert row == {
        "cases_checked": 64,
        "both_endpoints_equal_A": True,
        "mixed_states_are_skew_cuts": True,
        "proved": True,
    }
    for a in (-1, 1):
        for r in (-1, 1):
            block = equal_endpoint_block(a, r)
            assert block[0][0] * block[0][1] * block[1][0] * block[1][1] == -1
            assert sum(value * value for row in block for value in row) == 4
            assert sum(block[0][j] * block[1][j] for j in range(2)) == 0


def test_every_A_R_pair_through_order_three_satisfies_all_exact_forms() -> None:
    row = exhaustive_small_order_certificate(3)
    assert row["A_R_pairs_by_order"] == {1: 1, 2: 4, 3: 64}
    assert row["A_R_pairs_checked"] == 69
    assert row["local_block_cases"] == 64
    assert row["exact_integer_arithmetic"] is True
    assert row["proved_for_enumerated_orders"] is True


def test_a_nontrivial_order_four_pair_has_identical_frame_cut_pair_minimax() -> None:
    a = tuple(symmetric_signings(4))[37]
    r = tuple(skew_signings(4))[22]
    row = verify_equal_endpoint_algebra(a, r)
    assert row["cut_cases_checked"] == 16 * 16
    assert row["lift_cases_checked"] == 16 * 16
    assert row["endpoint_zero_equals_A"] is True
    assert row["endpoint_one_equals_A"] is True
    assert row["K_by_frames"] == row["K_by_cuts"] == row["K_by_pairs"]
    assert equal_endpoint_k_by_frames(a, r) == equal_endpoint_k_by_cuts(a, r)
    assert equal_endpoint_k_by_cuts(a, r) == equal_endpoint_k_by_pairs(a, r)


def test_equal_endpoints_make_every_hereditary_cut_bound_automatic() -> None:
    # This is a theorem for every A; replay all 64 order-four signings.
    for a in symmetric_signings(4):
        row = hereditary_endpoint_certificate(a)
        assert row["max_P_A_T_plus_P_A_Tc"] <= row["Phi_A"]
        assert row["max_N_A_T_plus_N_A_Tc"] <= row["Phi_A"]
        assert row["max_cut_abs_internal_plus_abs_cross"] <= row["Phi_A"]
        assert row["hereditary_endpoint_bounds_automatic"] is True
        assert row["proved_for_input"] is True


@pytest.mark.parametrize(
    ("length", "expected"),
    [
        (0, Fraction(0)),
        (1, Fraction(1)),
        (2, Fraction(1)),
        (3, Fraction(3, 2)),
        (4, Fraction(3, 2)),
        (5, Fraction(15, 8)),
        (6, Fraction(15, 8)),
    ],
)
def test_simple_walk_absolute_mean_closed_formula(
    length: int, expected: Fraction
) -> None:
    assert simple_walk_absolute_mean(length) == expected


def test_skew_infinity_to_one_floor_is_the_exact_random_sign_average() -> None:
    # The average is independent of the orientation.  Check every order-four R.
    assert skew_norm_floor_exact(4) == 6
    for r in skew_signings(4):
        row = skew_norm_floor_certificate(r)
        assert row["exact_average_over_y"] == 6
        assert row["infinity_to_one_norm"] >= 6
        assert row["floor_verified"] is True


def test_n5_cycle_chord_matrix_refutes_only_the_zero_error_disk() -> None:
    row = n5_disk_counterexample()
    assert row["Phi_A"] == 4
    assert row["energy_alphabet"] == [-4, 0, 4]
    assert row["E_Q_squared"] == 10
    assert row["moment_and_parity_force_Phi_at_least_4"] is True
    assert row["maximizer_basis_energies"] == [4] * 5
    assert row["determinant_V"] == 16
    assert row["skew_signings_checked"] == 2**10
    assert row["skew_signings_passing_zero_error_anchor_constraints"] == 0
    assert row["minimum_max_anchor_abs_x_R_y"] == 4
    assert row["minimum_anchor_disk_value_squared"] == 20
    assert row["disk_radius_squared"] == 16
    assert row["zero_error_disk_impossible"] is True
    assert row["refutes_only_zero_error_disk"] is True
    assert row["does_not_refute_asymptotic_o_n_cubed_error"] is True

    a = n5_cycle_chord_signing()
    basis = n5_maximizer_basis()
    assert [quadratic_energy(a, vector) for vector in basis] == [4] * 5
    assert all(not zero_error_disk_holds(a, r) for r in skew_signings(5))


@pytest.mark.parametrize("bad_order", [0, 5])
def test_small_order_replay_rejects_accidental_large_or_empty_census(
    bad_order: int,
) -> None:
    with pytest.raises(ValueError, match="1<=max_order<=4"):
        exhaustive_small_order_certificate(bad_order)
