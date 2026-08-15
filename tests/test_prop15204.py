"""Prop 15.204: PSD ker characterization; residual_i stays open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15204 import (
    free_e_bound_proved_general,
    free_e_dual_Sstar_sc,
    residual_i_closed_via_204,
    theorem_eight_over_p_beats_need,
    theorem_psd_ker_characterization,
)


def test_psd_ker_proved():
    r = theorem_psd_ker_characterization()
    assert r["proved"] is True


def test_eight_over_p_beats_need():
    r = theorem_eight_over_p_beats_need()
    assert r["proved"] is True
    assert r["by_p_sample"]["5"]["beats"] is True


def test_Sstar_sc_p5_exact():
    d = free_e_dual_Sstar_sc(5)
    assert d["ok"] is True
    assert d["Sstar"] == "123/7"
    assert d["blocks_if_ker_sc"] is True
    assert d["Sstar_le_n_minus_2"] is True
    # free-e ub = α · 123/7 = 369/455
    assert Fraction(d["free_e_ub"]) == Fraction(369, 455)


def test_Sstar_sc_p7_exact():
    d = free_e_dual_Sstar_sc(7)
    assert d["ok"] is True
    assert d["Sstar"] == "3912/113"
    assert d["blocks_if_ker_sc"] is True
    assert d["Sstar_le_n_minus_2"] is True


def test_predicates_remain_open():
    assert free_e_bound_proved_general() is False
    assert residual_i_closed_via_204() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
