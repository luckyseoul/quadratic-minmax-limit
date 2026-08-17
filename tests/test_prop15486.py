"""Tests for Prop 15.486 — 15.276 family is TR, not a name of c."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15486 import (
    enumerate_family,
    family_size,
    family_size_drop_pm1,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_family_size():
    A = prove_A()
    assert A["proved"] is True
    assert family_size(5) == 100
    assert family_size(7) == 1176
    assert family_size_drop_pm1(5) == 25 != 100


def test_p5_is_NL_TR():
    B = prove_B()
    assert B["proved"] is True
    e = enumerate_family(5)
    assert e["n_max"] == e["n_tr"] == e["n_unique"] == 100
    assert e["n_fullge2"] == 0


def test_p7_is_TR_not_c():
    C = prove_C()
    assert C["proved"] is True
    assert C["c_family"] == 4
    assert C["c_live"] == 19
    e = enumerate_family(7)
    assert e["n_unique"] == 1176
    assert e["n_tr"] == 1176
    assert e["n_fullge2"] == 0


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.486"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["p7_family_is_TR_not_c"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
