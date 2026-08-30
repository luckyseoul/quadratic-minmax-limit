"""Tests for Prop 15.481 — W_c is the E-sum; named K_++ recovers Q++."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15478 import a_E_sq
from e1_gmin_m4_prop15481 import (
    E_card_minus_1,
    Kpp_named,
    Wc_exact,
    Wc_named,
    main,
    pred_Q,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_Wc_is_E_sum():
    A = prove_A()
    assert A["proved"] is True
    F = _kernel_field(5)
    assert abs(Wc_exact(F)[0] - (2 * 5 - a_E_sq(5))) < 1e-10
    assert E_card_minus_1(5) == 31
    assert abs(Wc_exact(F)[0] - 31) > 1.0
    assert abs(Wc_exact(F) - Wc_named(F, 5)).max() < 1e-10


def test_Kpp_identity():
    B = prove_B()
    assert B["proved"] is True
    assert B["rows"]["5"]["err"] < 1e-10
    assert B["rows"]["5"]["err_drop38"] > 0.2


def test_Q_recovered():
    C = prove_C()
    assert C["proved"] is True
    F = _kernel_field(5)
    q = pred_Q(F, Kpp_named(F, 5), 5)
    assert abs(q - float(LIVE_QPP[5])) < 1e-8
    qdrop = pred_Q(F, Kpp_named(F, 5, drop_c=True), 5)
    assert abs(qdrop - float(LIVE_QPP[5])) > 0.05


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.481"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["Wc_is_E_sum"] is True
    assert out["proved"]["named_K_recovers_Qpp"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
