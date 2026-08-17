"""Tests for Prop 15.360 — energy conditioning; cubics exhaust p=5."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15360 import (
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_E,
    prove_open,
)


def test_energy_identity():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["id_ok"] is True
    assert A["rows"]["7"]["E_ok"] is True


def test_p5_count_is_130_not_65():
    B = prove_B()
    assert B["proved"] is True
    assert B["n_energy_p5"] == 130
    assert B["n_energy_p5"] != 65


def test_quadrics_are_exactly_1d():
    C = prove_C()
    assert C["proved"] is True
    assert C["quad_unique"]["5"] == 30
    assert C["quad_unique"]["7"] == 140


def test_cubics_exhaust_p5():
    D = prove_D()
    assert D["proved"] is True
    assert D["cubic_unique_5"] == 130
    assert D["cubic_unique_5"] != 30


def test_nl_line_sum_type():
    E = prove_E()
    assert E["proved"] is True
    assert E["const1_1d"] == [5]
    assert E["const1_nl"] == [3]


def test_floor_still_open():
    F = prove_open()
    assert F["proved"] is False
    assert F["D_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["energy_identity"] is True
    assert out["proved"]["p5_count_130"] is True
    assert out["proved"]["quadrics_are_1d"] is True
    assert out["proved"]["cubics_exhaust_p5"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
