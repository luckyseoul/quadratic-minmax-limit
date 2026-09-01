"""Fail-when-wrong tests for Proposition 15.752."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
from e1_gmin_m4_prop15734 import BRANCH_B2, BRANCH_P1_LAST, BRANCH_P3_LAST
from e1_gmin_m4_prop15752 import (
    band_arithmetic,
    band_branch_exclusion,
    band_hard_residue_certificate,
    band_maximum_t,
    p19_sharp_mechanism_witness,
    p_plus_nine_boolean_certificate,
    p_plus_nine_height_certificate,
    p_plus_nine_local_exclusion,
    proposition_15752,
    residual_band_exclusion,
    write_evidence,
)


ROOT = Path(__file__).resolve().parents[1]


def test_p_plus_nine_height_contradiction_covers_both_congruence_classes() -> None:
    for p in (23, 29, 31, 37, 43):
        row = p_plus_nine_height_certificate(p)
        assert Fraction(row["height_lower_bound"]) > 3
        assert Fraction(row["paired_cube_average_upper_bound"]) < Fraction(3, 4)
        assert row["half_mean_cube_maximum_upper_bound"] == 3
        assert row["contradiction"] is True


def test_p_plus_nine_boolean_branch_reuses_only_the_fixed_four_bit_catalog() -> None:
    for p in (23, 29, 31, 37, 43, 47):
        row = p_plus_nine_boolean_certificate(p)
        low, high = map(Fraction, row["target_strict_bracket"])
        target = Fraction(row["target_density"])
        assert Fraction(row["largest_zero_influence_class_complement_bound"]) < 7
        assert row["junta_coordinates_at_most"] == 6
        assert row["cube_coordinates_actually_needed_at_most"] == 4
        assert low < target < high
        assert row["target_absent"] is True
        assert row["seven_gap_second_derivative_at_p_23"] > 0
        assert p_plus_nine_local_exclusion(p)["excluded"] is True


def test_p19_is_a_real_local_threshold_not_an_unchecked_small_prime() -> None:
    row = p19_sharp_mechanism_witness()
    assert row["layer_values"] == [3, 1, 0, 0, 1]
    assert row["mean"] == "7/19"
    assert row["scaled_mass_4p_E_B"] == 28
    assert row["equals_p_plus_9"] is True
    assert row["is_only_a_local_quadratic_not_a_residual_graph"] is True


def test_band_endpoints_are_exact_and_keep_an_isolated_chart() -> None:
    expected = {23: 8, 29: 10, 31: 12, 37: 14, 43: 18}
    for p, t_max in expected.items():
        assert band_maximum_t(p) == t_max
        row = band_arithmetic(p, t_max)
        assert row["closed_t_interval"] == [4, t_max]
        assert row["guaranteed_isolated_vertices"] > 0
        assert row["maximum_low_phase_one_mean"] < row[
            "next_nonbaseline_phase_one_floor"
        ]
        assert row["maximum_endpoint_baseline_lift_excess"] < p - 3
        assert row["boundary_size_hypothesis_used"] is False


def test_hard_residues_reduce_to_only_the_three_named_branches() -> None:
    p1 = band_hard_residue_certificate(29, 10)
    assert p1["possible_branches"] == [BRANCH_B2, BRANCH_P1_LAST]
    assert p1["u_0_through_t_rows"][0]["surviving_branch"] == BRANCH_P1_LAST
    assert all(
        row["surviving_branch"] is None for row in p1["u_0_through_t_rows"][1:]
    )

    p3 = band_hard_residue_certificate(23, 8)
    assert p3["possible_branches"] == [BRANCH_B2, BRANCH_P3_LAST]
    assert all(row["surviving_branch"] is None for row in p3["u_0_through_t_rows"])
    assert p3["equal_mean_p3_endpoint_cells_cannot_mix"] is True


def test_each_branch_forces_a_forbidden_next_mass() -> None:
    for p, t in ((23, 4), (23, 8), (31, 12)):
        for branch in (BRANCH_B2, BRANCH_P3_LAST):
            row = band_branch_exclusion(p, t, branch)
            assert row["forced_next_scaled_mean"] == p + 9
            assert row["next_cell_forced_to_b_zero"] is True
            assert row["local_mass_exclusion"]["proved"] is True
            assert row["branch_excluded"] is True

    for branch, target in ((BRANCH_B2, 38), (BRANCH_P1_LAST, 36)):
        row = band_branch_exclusion(29, 4, branch)
        assert row["forced_next_scaled_mean"] == target
        assert row["branch_excluded"] is True


def test_fifth_shell_and_full_claimed_band_are_closed_boundary_independently() -> None:
    for p in (23, 29, 31, 37, 43):
        for t in (4, band_maximum_t(p)):
            row = residual_band_exclusion(p, t)
            assert row["original_k"] == 4 * p + 2 * t
            assert row["all_boundary_sizes_excluded"] is True
            assert row["finite_prime_or_configuration_census_used"] is False
            assert row["residual_ii_layer_excluded"] is True


def test_prop15752_flips_a_real_infinite_layer_but_not_the_global_gate(
    tmp_path: Path,
) -> None:
    row = proposition_15752()
    assert row["k_eq_4p_plus_8_closed_for_every_prime_p_ge_23"] is True
    assert row["boundary_size_hypothesis_used"] is False
    assert row["finite_prime_or_configuration_census_used"] is False
    assert row["residual_ii_k_ge_4p_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert residual_ii_k_ge_4p_ND_closed() is False

    expected = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15752.json").read_text()
    )
    assert expected == row
    replay = tmp_path / "prop15752.json"
    assert write_evidence(replay) == replay
    assert json.loads(replay.read_text()) == row
    assert not list(tmp_path.glob("*.tmp"))


@pytest.mark.parametrize("bad_p", [13, 17, 19, 21, 25, 27])
def test_local_theorem_rejects_out_of_scope_primes(bad_p: int) -> None:
    with pytest.raises(ValueError):
        p_plus_nine_local_exclusion(bad_p)


@pytest.mark.parametrize("p,t", [(23, 3), (23, 9), (29, 11), (31, 13)])
def test_band_api_rejects_layers_outside_the_proved_interval(p: int, t: int) -> None:
    with pytest.raises(ValueError):
        residual_band_exclusion(p, t)
