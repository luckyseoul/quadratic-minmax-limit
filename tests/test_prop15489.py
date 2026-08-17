"""Tests for Prop 15.489 — leftover ridges are 2p(χ*μ)."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15488 import csp_count
from e1_gmin_m4_prop15489 import F_named, main, prove_A, prove_B, prove_open


def test_p3_p5_csp_with_chi_family():
    A = prove_A()
    assert A["proved"] is True
    assert A["N3"] == 6
    assert A["N5"] == 130
    assert csp_count(5, F_named(5), flip_T1=True) != 130
    assert csp_count(5, F_named(5, drop_level=True)) != 130


def test_p7_chi_family_necessary():
    B = prove_B()
    assert B["proved"] is True
    assert len(F_named(7)) == 638
    assert len(F_named(7, drop_chi=True)) == 78
    assert len(F_named(7, drop_chi=True)) != len(F_named(7))


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["N5"] == 130
    assert out["B"]["nF"] == 638
