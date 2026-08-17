"""Tests for Prop 15.307 — signatures classify H_+; E[M²] mixture."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15307 import (
    classify_p,
    main,
    prove_halfspaces,
    prove_open,
    prove_p5_single_NL,
    prove_p7_orbits,
)


def test_halfspace_count_is_n1d():
    A = prove_halfspaces()
    assert A["proved"] is True
    for p in (3, 5, 7):
        assert A["by_p"][str(p)]["n_hs"] == n_1d(p)


def test_p5_NL_is_p2_times_p_minus_1():
    B = prove_p5_single_NL()
    assert B["proved"] is True
    assert B["p5_NL"] == 5 * 5 * 4
    assert B["p5_allR"] == 100


def test_p7_is_not_single_NL_type():
    B = prove_p5_single_NL()
    assert B["p7_n_D_types"] >= 6
    assert abs(B["p7_one_NL_formula"] - B["p7_live_EM2"]) > 1.0


def test_p7_orbit_sizes_sum():
    C = prove_p7_orbits()
    assert C["proved"] is True
    assert C["Nplus"] == 5726
    assert C["Nplus_named_sum"] == 5726
    assert C["unit_p2_pm1"] == 294


def test_EM2_from_signatures():
    for p in (5, 7):
        r = classify_p(p)
        assert r["E_M2_match"] is True
        assert abs(r["E_M"] - r["E_M_named"]) < 1e-9


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Hplus_named"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert D["lambda_ge_6"] is False
    assert D["Q_tau_named_in_p"] is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["halfspaces_are_n1d"] is True
    assert out["proved"]["p5_NL_is_p2_pm1"] is True
    assert out["proved"]["p7_orbit_sum"] is True
    assert out["proved"]["Hplus_named"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
