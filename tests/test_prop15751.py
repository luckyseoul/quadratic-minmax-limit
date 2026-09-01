"""Fail-when-wrong tests for Proposition 15.751."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
from e1_gmin_m4_prop15751 import (
    EXPECTED_HISTOGRAM,
    EXPECTED_SHA256,
    affine_slice_support_certificate,
    cube_half_mean_height_certificate,
    density_profile_certificate,
    exact_four_cube_catalog,
    height_at_least_two_certificate,
    height_one_junta_certificate,
    proposition_15751,
    write_evidence,
)


ROOT = Path(__file__).resolve().parents[1]


def test_dimension_free_half_mean_height_bound_is_sharp() -> None:
    row = cube_half_mean_height_certificate()
    assert row["maximum_upper_bound"] == 3
    assert row["sharp_example"]["layer_values"] == [3, 1, 0, 0, 1]
    assert row["sharp_example"]["mass"] == 8
    assert row["proved"] is True


def test_corrected_influence_normalization_forces_six_coordinate_junta() -> None:
    for p in (29, 37, 41, 53, 61):
        support = affine_slice_support_certificate(p)
        junta = height_one_junta_certificate(p)
        assert support["minimum_support_density"] == str(
            Fraction(p - 3, 2 * (p - 2))
        )
        assert junta["dictator_normalization_check"] == (
            "sum I_ij=(p/2)*Var(x_i)"
        )
        assert Fraction(junta["junta_coordinate_bound_exact"]) < 7
        assert junta["junta_coordinates_at_most"] == 6
        assert junta["cube_coordinates_actually_needed_at_most"] == 4


def test_height_at_least_two_branch_contradicts_cube_height_bound() -> None:
    for p in (29, 37, 41):
        row = height_at_least_two_certificate(p)
        assert Fraction(row["height_lower_bound"]) > 3
        assert Fraction(row["paired_cube_average_upper_bound"]) < Fraction(3, 4)
        assert row["contradiction"] is True


def test_exact_cpu_replay_matches_all_accelerator_constants() -> None:
    row = exact_four_cube_catalog()
    assert row["tables_checked"] == 65_536
    assert row["valid_tables"] == 222
    assert row["valid_table_signature_sha256"] == EXPECTED_SHA256
    assert row["packed_layer_signature_histogram"] == {
        str(key): value for key, value in sorted(EXPECTED_HISTOGRAM.items())
    }
    assert len(row["profiles"]) == 14


def test_all_four_bit_density_profiles_miss_the_branch_target() -> None:
    for p in (29, 37, 41, 53, 61):
        row = density_profile_certificate(p)
        low, high = map(Fraction, row["nearest_strict_bracket"])
        target = Fraction(row["target_density"])
        assert low < target < high
        assert row["target_absent"] is True


def test_prop15751_closes_exactly_the_fourth_shell_front(tmp_path: Path) -> None:
    row = proposition_15751()
    assert row["generic_branch_B_t3_p_ge_29_closed"] is True
    assert row["finite_prime_census_used"] is False
    assert row["fixed_four_bit_certificate_used"] is True
    assert row["residual_ii_k_ge_4p_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert residual_ii_k_ge_4p_ND_closed() is False

    expected = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15751.json").read_text()
    )
    assert expected == row
    replay = tmp_path / "prop15751.json"
    assert write_evidence(replay) == replay
    assert json.loads(replay.read_text()) == row
    assert not list(tmp_path.glob("*.tmp"))


@pytest.mark.parametrize("bad_p", [13, 17, 25, 31, 45])
def test_uniform_branch_api_rejects_out_of_scope_values(bad_p: int) -> None:
    with pytest.raises(ValueError, match="prime p>=29"):
        height_one_junta_certificate(bad_p)
