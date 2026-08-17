"""Tests for Prop 15.413 — NL P(++)=4k(k−1) at p=5, dies at p=7."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import nonsquare_NN_mass
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15413 import four_k_km1, k_named, main, prove_open


def test_named_masses():
    assert k_named(5) == 10
    assert k_named(7) == 21
    assert four_k_km1(5) == 360
    assert four_k_km1(7) == 1680
    assert nonsquare_NN_mass(5) == 300
    assert four_k_km1(5) != nonsquare_NN_mass(5)


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["P_ns_is_W_mass"] is True
    assert out["proved"]["p5_NL_P_is_4kkm1"] is True
    assert out["proved"]["p7_NL_P_not_4kkm1"] is True
    assert out["algebra"]["B"]["P_pp"] == "360"
    assert out["algebra"]["B"]["L_pp"] == "24"
    assert out["algebra"]["C"]["P_pp"] == "216580/57"
    assert out["algebra"]["C"]["L_pp"] == "38143/342"
    assert Fraction(out["algebra"]["C"]["P_pp"]) != four_k_km1(7)
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.413"
