"""Tests for Prop 15.517 — p=5 H_+ W is 2-orbit; p=7 is not."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15323 import EW_named, Var_2orb
from e1_gmin_m4_prop15517 import main, prove_A, prove_B, prove_open, support_named


def test_p5_W_is_two_orbit():
    A = prove_A()
    assert A["proved"] is True
    assert A["fold"]["var"] == str(Var_2orb(5))
    assert A["fold"]["mean"] == str(EW_named(5))
    assert tuple(A["fold"]["support"]) == support_named(5)


def test_p7_fails_four_point_and_Var_2orb():
    B = prove_B()
    assert B["proved"] is True
    assert B["fold"]["nuniq"] != 4
    assert B["fold"]["var"] == "21504/409"
    assert B["fold"]["var"] != str(Var_2orb(7))


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["phi_F_ge_6"] is False
