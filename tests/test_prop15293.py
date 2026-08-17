"""Tests for Prop 15.293 — Q_NL is not Hoffman-on-15.290-types."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import Q_1d_pm1_over_q2, Q_1d_sub_over_q2
from e1_gmin_m4_prop15293 import main, prove_open, prove_p11_spread_from_scratch


def test_Q_1d_still_named():
    assert Q_1d_pm1_over_q2(5).numerator != 0
    assert Q_1d_sub_over_q2(7).denominator == 5


def test_p11_type_does_not_determine_Q_H():
    B = prove_p11_spread_from_scratch()
    assert B["proved"] is True
    assert B["n_big_spreads"] >= 1
    # fail-when-wrong: type-constancy
    assert any(sp > 1e-4 for sp in B["spreads"].values())


def test_Q_NL_and_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["H_equals_NL_only_p5"] is True
    assert out["proved"]["p11_type_does_not_determine_Q_H"] is True
    assert out["proved"]["Q_NL_named_in_p"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    # p=7 invert: H and NL differ
    rel7 = out["algebra"]["A"]["sample"]["7"]["max_rel"]
    rel5 = out["algebra"]["A"]["sample"]["5"]["max_rel"]
    assert rel5 < 1e-8
    assert rel7 > 1e-4
