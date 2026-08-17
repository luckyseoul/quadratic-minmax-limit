"""Tests for Prop 15.324 — fiber-average W; U-Parseval; G_U on Ω."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15324 import (
    EW_named,
    main,
    prove_GU_spread,
    prove_fiber_avg,
    prove_open,
    prove_parseval,
    sum_W_squares_named,
)


def test_fiber_average_is_EW_but_W_c_varies():
    A = prove_fiber_avg()
    assert A["proved"] is True
    assert A["W_c_constant_false"] is True
    assert sum_W_squares_named(5) == 60
    assert sum_W_squares_named(5) == 2 * EW_named(5)
    assert sum_W_squares_named(7) == 3 * EW_named(7)
    assert A["by_p"]["5"]["err"] < 1e-6
    assert A["by_p"]["7"]["err"] < 1e-6
    assert A["by_p"]["5"]["per_c"][1]["nuniq"] > 1


def test_odd_U_hats_vanish_parseval():
    B = prove_parseval()
    assert B["proved"] is True
    assert B["odd_hats_false"] is True
    assert B["by_p"]["5"]["odd_max"] < 1e-8
    assert abs(B["by_p"]["5"]["even_hats"][2] - 750 / 13) < 1e-4
    assert abs(B["by_p"]["5"]["even_hats"][2] - B["by_p"]["5"]["even_hats"][4]) < 1e-6
    assert B["by_p"]["5"]["parseval_ok"] is True
    assert B["by_p"]["7"]["parseval_ok"] is True


def test_GU_spreads_on_Omega_not_on_squares():
    C = prove_GU_spread()
    assert C["proved"] is True
    assert C["GU_const_on_Omega_false"] is True
    assert C["Omega_spread"] > 5
    assert C["square_spread"] < 1e-6


def test_majorant_and_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["A_square_named"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["fiber_avg"] is True
    assert out["proved"]["U_parseval"] is True
    assert out["proved"]["GU_spreads_on_Omega"] is True
    assert out["proved"]["A_square_named"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
