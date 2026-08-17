"""Tests for Prop 15.476 — Δ² from direction type."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15476 import (
    is_triple_rainbow,
    live_delta2_and_sigs,
    main,
    prove_A,
    prove_B,
    prove_open,
    rainbow_tuple,
    regular_tuple,
)


def test_type_determines_delta2():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["split"] == 0
    assert A["rows"]["7"]["split"] == 0
    assert A["rows"]["7"]["nfull_split"] >= 1


def test_triple_rainbow():
    assert rainbow_tuple(5) == (0, 1, 2, 3, 4)
    assert regular_tuple(5) == (2, 2, 2, 2, 2)
    B = prove_B()
    assert B["proved"] is True
    assert B["rows"]["5"]["n_triple_rainbow"] == 100
    assert B["rows"]["5"]["tr_d2"] == [4500.0]
    assert B["rows"]["7"]["n_triple_rainbow"] == 1176
    assert B["rows"]["7"]["n_parunion"] == 140


def test_nfull_is_not_the_name():
    d2, sigs, n_full = live_delta2_and_sigs(7)
    ones = [round(float(d2[i]), 4) for i, n in enumerate(n_full) if n == 1]
    assert len(set(ones)) >= 2


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.476"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["delta2_fn_of_dirtype"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
