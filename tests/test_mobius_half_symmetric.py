from __future__ import annotations

from e1_gmin_m4_mobius_half_symmetric import (
    branch_c_capacity_ledger,
    forced_symmetric_certificate,
    ordinary_half_certificate,
    paley_tau_certificate,
    two_trade_origin_cancellation_certificate,
)


def test_ordinary_half_has_exact_parallel_formula_and_central_other_rows():
    for p in (7, 11):
        out = ordinary_half_certificate(p)
        assert out["edge_count"] == p - 1
        assert out["sum_parallel_counts"] == p - 1
        assert out["parallel_counts"][0] == 1
        assert out["all_other_K_rows_centrally_symmetric"] is True
        assert out["proved"] is True


def test_paley_tau_is_the_exact_even_quartic_trace():
    for p in (7, 11):
        out = paley_tau_certificate(p)
        assert out["tau_1"] == out["epsilon_L"] == 1
        assert out["domain_sum_S"] == (
            out["complete_quartic_character_sum"] - out["epsilon_L"]
        )
        assert out[
            "tau_t_equals_tau_minus_t_when_both_parameters_are_in_domain"
        ] is True
        assert out["proved"] is True


def test_forced_pair_total_and_selected_half_have_exact_rows():
    for p in (7, 11):
        out = forced_symmetric_certificate(p)
        assert out["used_inversion_orbits"] == p - 1
        assert out["selected_graph_edges"] == p - 1
        assert out["forced_symmetric_K_rows_central"] is True
        assert out["selected_nonchosen_K_rows_central"] is True
        assert out["selected_antisymmetric_image"] == "A_j in L and zero elsewhere"
        assert out["coupled_symmetric_target_realized"] is False
        assert out["proved"] is True


def test_two_arbitrary_nonzero_stars_can_cancel_exactly_one_origin_orbit():
    cases = (
        (7, 0, 1, 1, 2),
        (7, 0, 7, 2, 3),
        (11, 2, 11, 3, 5),
    )
    for args in cases:
        out = two_trade_origin_cancellation_certificate(*args)
        p = args[0]
        assert out["shared_inversion_orbits"] == 1
        assert out["same_sign_shared_orbits"] == 0
        assert out["opposite_sign_shared_orbits"] == 1
        assert out["nonzero_orbits_after_cancellation"] == 2 * (p - 1) - 2
        assert out["source_stays_ternary"] is True
        assert out["both_direction_targets_stay_exact"] is True
        assert out["proved"] is True


def test_upper_endpoint_one_edge_gap_is_defeated_not_promoted_to_a_theorem():
    out = branch_c_capacity_ledger(31)
    assert out["all_centers_nonzero_disjoint_trade_edges"] == 480
    assert out["H_edge_count_interval"] == [261, 479]
    assert out["one_cancelled_pair_trade_edges"] == 478
    assert out["universal_support_lower_bound_proved"] is False
    assert out["capacity_contradiction_proved"] is False
    assert out["proved"] is True
