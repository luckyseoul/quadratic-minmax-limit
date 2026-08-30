"""Tests for Prop 15.372 — Q++_form in S0 box; 6-net sanity."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15326 import Q_lo_named, even_S_p1
from e1_gmin_m4_prop15361 import n_1d_plus_hoff
from e1_gmin_m4_prop15367 import A_num
from e1_gmin_m4_prop15371 import D_form, D_form_drop_binom
from e1_gmin_m4_prop15372 import (
    Qpp_form,
    box_hi,
    box_lo,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_six_net_recovers_live_5_7():
    A = prove_A()
    assert A["proved"] is True
    assert A["p5_3net"] == 130
    assert A["p7_4net"] == 5726
    assert A["p7_4net"] != n_1d_plus_hoff(7)
    assert A["replaced_3net"] is True


def test_Qform_in_interval_and_box():
    B = prove_B()
    assert B["proved"] is True
    Q5 = Qpp_form(5)
    Q7 = Qpp_form(7)
    assert Q5 == Fraction(48, 13)
    assert Q7 == Fraction(1544, 409)
    assert Q_lo_named(5) <= Q5 <= Qpp_floor_ub(5)
    assert Q_lo_named(7) <= Q7 <= Qpp_floor_ub(7)
    assert box_lo(5) <= Q5 <= box_hi(5)
    assert box_lo(7) <= Q7 <= box_hi(7)
    # fail-when-wrong: drop binomial
    Qbad = Fraction(8 * A_num(5), D_form_drop_binom(5))
    assert Qbad == 16
    assert not (Q_lo_named(5) <= Qbad <= Qpp_floor_ub(5))
    assert D_form(5) == 13
    assert D_form_drop_binom(5) == 3


def test_p1_S0_even_ge_6():
    C = prove_C()
    assert C["proved"] is True
    min5 = min(even_S_p1(5, Qpp_form(5)).values())
    assert min5 == Fraction(80, 13)
    assert min5 >= 6
    minbad = min(even_S_p1(5, Fraction(16)).values())
    assert minbad < 6


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["D_named_in_p"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["six_net_recovers_5_7"] is True
    assert out["proved"]["Qform_in_s0_and_box"] is True
    assert out["proved"]["p1_S0_even_ge_6_at_Qform"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
