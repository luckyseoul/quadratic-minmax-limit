"""Tests for Prop 15.294 — Q_rest; familywise S≥6 is false."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import S_1d_over_q2_ntriv_on_Fp
from e1_gmin_m4_prop15294 import main, prove_familywise_floor_false, prove_open, prove_rest_Q


def test_rest_empty_p5_nonempty_p7():
    A = prove_rest_Q()
    assert A["proved"] is True
    assert A["sample"]["5"]["n_rest"] == 0
    assert A["sample"]["7"]["n_rest"] > 0
    assert A["invert_hit"] < A["invert_tot"]


def test_Q_rest_not_Q_H_and_type_constant():
    A = prove_rest_Q()
    assert A["sample"]["7"]["max_rel_rest_H"] > 1e-4
    for row in A["sample"]["7"]["Q_rest"].values():
        assert row["spread"] < 1e-8


def test_familywise_floor_is_false():
    B = prove_familywise_floor_false()
    assert B["proved"] is True
    assert B["sample"]["5"]["Smin"]["H"] < 1e-6
    assert B["sample"]["7"]["Smin"]["rest"] < 6
    assert B["sample"]["5"]["Smin"]["1d"] > 6
    assert B["sample"]["7"]["Smin"]["full"] > 6
    assert B["invert_hit"] < B["invert_tot"]


def test_S_1d_formula_still_above_six():
    assert float(S_1d_over_q2_ntriv_on_Fp(5)) > 6
    assert float(S_1d_over_q2_ntriv_on_Fp(7)) > 6


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["rest_Q_type_constant_p7"] is True
    assert out["proved"]["familywise_S_ge_6_false"] is True
    assert out["proved"]["Q_rest_named_in_p"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
