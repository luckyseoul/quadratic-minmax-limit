"""Tests for Prop 15.356 — Q_1D++ named; ensemble is the 1D/NL mix."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS, LIVE_QPP
from e1_gmin_m4_prop15355 import Q_AP_named
from e1_gmin_m4_prop15356 import (
    Q_1d_pp_named,
    Q_NL_from_live,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    w_hs,
)


def test_q1d_named_and_not_ensemble():
    assert Q_1d_pp_named(5) == Fraction(16, 9) == Q_AP_named(5)
    assert Q_1d_pp_named(7) == Fraction(104, 15)
    assert Q_1d_pp_named(7) != Q_AP_named(7)
    assert Q_1d_pp_named(7) != LIVE_QPP[7]
    assert Q_1d_pp_named(5) != LIVE_QPP[5]
    assert Q_1d_pp_named(11) == Fraction(344, 27)
    assert Q_1d_pp_named(11) != Q_AP_named(11)


def test_mix_recovers_live():
    assert prove_A()["proved"] is True
    assert prove_B()["proved"] is True
    assert prove_C()["proved"] is True
    assert w_hs(5) == Fraction(n_1d(5), HPLUS[5]) == Fraction(3, 13)
    assert w_hs(7) == Fraction(10, 409)
    assert Q_NL_from_live(5) == Fraction(64, 15)
    assert w_hs(7) != Fraction(n_1d(7), n_1d(7) + 7 * 7 * 6)


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Q_1d_pp_named_in_p"] is True
    assert out["proved"]["ensemble_is_1d_NL_mix"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
