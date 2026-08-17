"""Tests for Prop 15.297 — Aut partners fail; Hoffman X_A vanish is p=5-only."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15297 import (
    main,
    prove_XA_vanishes_on_Hoffman,
    prove_aut_does_not_sign,
    prove_open,
)


def test_aut_pair_min_below_floor():
    A = prove_aut_does_not_sign()
    assert A["proved"] is True
    assert A["maps"]["signed_inv"]["below_floor"] is True
    assert A["invert_hit"] < A["invert_tot"]
    assert A["maps"]["signed_inv"]["pair_min"] < A["floor"]


def test_Hoffman_XA_zero_only_p5():
    B = prove_XA_vanishes_on_Hoffman()
    assert B["proved"] is True
    assert B["sample"]["5:2"]["vanishes"] is True
    assert B["sample"]["7:40"]["vanishes"] is False
    assert B["invert_hit"] < B["invert_tot"]


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
    assert out["proved"]["aut_pair_below_floor"] is True
    assert out["proved"]["Hoffman_XA_zero_only_p5"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
