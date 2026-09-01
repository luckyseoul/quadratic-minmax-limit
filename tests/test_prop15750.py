"""Focused tests for Proposition 15.750."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
from e1_gmin_m4_prop15275 import type_I_multilevel_bad_case_ND_closed
from e1_gmin_m4_prop15750 import (
    central_krawtchouk,
    distinguished_edge_normalization_certificate,
    general_prime_row,
    johnson_pair_constancy_certificate,
    nonsquare_parallel_two_certificate,
    parity_bias_exact_replay,
    proposition_15750,
    square_direction_kappa_candidates,
    square_direction_rigidity_certificate,
    type_I_multilevel_bad_case_closed_all_primes,
    uniform_general_prime_certificate,
    write_evidence,
)
from e1_main_chain_status import four_e1_units_closed


ROOT = Path(__file__).resolve().parents[1]


def test_square_direction_rigidity_and_doubled_p7_endpoint() -> None:
    assert square_direction_kappa_candidates(
        7, remove_doubled_edge_endpoint=False
    ) == [-1, 0]
    assert square_direction_kappa_candidates(7) == [0]
    for p in (11, 13, 17, 101):
        assert square_direction_kappa_candidates(p) == [0]
        johnson = johnson_pair_constancy_certificate(p)
        square = square_direction_rigidity_certificate(p)
        assert johnson["all_pair_coefficients_equal"] is True
        assert johnson["proved"] is True
        assert square["parallel_W_multiplicity"] == 3
        assert square["all_off_fibre_signed_block_sums"] == 0
        assert square["proved"] is True


def test_small_prime_canonical_edge_normalizes_by_signed_psl() -> None:
    for p in (5, 7):
        row = distinguished_edge_normalization_certificate(p)
        assert row["PSL_2_is_two_transitive"] is True
        assert row["ordered_edge_target"] == ["infinity", 0]
        assert row["canonical_matrix_edge_indices"] == [0, 1]
        assert row["signed_lift_preserves_bad_box_system"] is True
        assert row["proved"] is True


def test_central_krawtchouk_parity_bound_replays_exactly() -> None:
    for p in (5, 7, 11, 13, 17, 31):
        replay = parity_bias_exact_replay(p)
        assert replay["K_1"] == replay["K_2"]
        assert replay["recurrence_exact"] is True
        assert (
            replay["max_p_times_absolute_K"]
            <= replay["bias_bound_denominator"]
        )
        assert replay["proved"] is True
    assert central_krawtchouk(11, 1) == -42
    assert central_krawtchouk(11, 2) == -42


def test_nonsquare_parallel_two_has_only_means_four_and_six() -> None:
    positive = nonsquare_parallel_two_certificate(11, 1)
    negative = nonsquare_parallel_two_certificate(11, -1)
    assert positive["some_direction_has_P"] == 2
    assert positive["a_equals_2p_E_T"] == 4
    assert positive["mass_shortfall_from_all_P_at_least_3"] == 3
    assert negative["some_direction_has_P"] == 2
    assert negative["a_equals_2p_E_T"] == 6
    assert negative["mass_shortfall_from_all_P_at_least_3"] == 4


def test_general_rows_clear_parity_and_lift_floors() -> None:
    for p in (11, 13, 17, 31):
        for c_e, expected_mean in ((1, 4), (-1, 6)):
            row = general_prime_row(p, c_e)
            assert row["a_equals_2p_E_T"] == expected_mean
            assert row["parity_gap"] > 0
            assert row["T_is_constantly_even"] is True
            assert row["B_equals_T_over_2_is_nonzero"] is True
            assert row["lift_floor_gap"] > 0
            assert row["closed"] is True


def test_uniform_certificate_is_theorem_not_finite_census() -> None:
    row = uniform_general_prime_certificate()
    assert row["domain"] == "every prime p>=11"
    assert row["proof_kind"] == "uniform theorem; samples are replays only"
    assert row["bad_case_setup"]["W"] == "G+2e as an edge multiset"
    assert min(row["minimum_symbolic_gaps_at_p_11"].values()) > 0
    assert row["proved"] is True


def test_prop15750_flips_only_type_i_and_matches_evidence(
    tmp_path: Path,
) -> None:
    row = proposition_15750()
    assert row["type_I_multilevel_bad_case_ND_closed"] is True
    assert row["residual_ii_closed"] is False
    assert row["E1_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["L_status"] == "OPEN"
    assert row["finite_graph_census_used"] is False
    assert row["scipy_or_optimizer_theorem_dependency"] is False
    assert row["proved"] is True
    assert type_I_multilevel_bad_case_closed_all_primes() is True
    assert type_I_multilevel_bad_case_ND_closed() is True
    assert residual_ii_k_ge_4p_ND_closed() is False

    expected = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15750.json").read_text()
    )
    assert expected == row

    replay = tmp_path / "prop15750.json"
    assert write_evidence(replay) == replay
    assert json.loads(replay.read_text()) == row
    assert not list(tmp_path.glob("*.tmp"))


def test_main_chain_remains_open_only_because_other_units_are_open() -> None:
    units = four_e1_units_closed()
    assert units["type_I_multilevel"] is True
    assert units["residual_ii_k_ge_4p"] is False
    assert units["closed"] is False


@pytest.mark.parametrize("bad_p", [9, 15, 2])
def test_general_theorem_rejects_out_of_domain_inputs(bad_p: int) -> None:
    with pytest.raises(ValueError, match="prime p>=11"):
        general_prime_row(bad_p, 1)
