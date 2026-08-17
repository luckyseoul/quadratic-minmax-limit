"""Tests for Prop 15.430 — L_par(R) named in p."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15414 import LIVE_LPAR_NL
from e1_gmin_m4_prop15416 import LIVE_LR, Lpar_R_flag, Lpar_R_indep_wrong
from e1_gmin_m4_prop15430 import (
    Lpar_R_naive_wrong,
    Lpar_R_named,
    S_from_J,
    S_named,
    S_shift_wrong,
    main,
    prove_open,
)


def test_S_and_fail():
    assert S_named(1, 5) == S_from_J(1, 5)
    assert S_named(2, 7) == S_from_J(2, 7)
    assert not all(S_shift_wrong(s, 5) == S_from_J(s, 5) for s in range(1, 5))


def test_named_and_fail():
    assert Lpar_R_named(5) == LIVE_LR[5] == Lpar_R_flag(5) == 24
    assert Lpar_R_named(7) == LIVE_LR[7] == Lpar_R_flag(7) == 133
    assert Lpar_R_named(11) == Lpar_R_flag(11) == Fraction(5379, 5)
    assert Lpar_R_indep_wrong(5) == Fraction(49, 2) != Lpar_R_named(5)
    assert Lpar_R_naive_wrong(7) == 150 != Lpar_R_named(7)
    assert Lpar_R_named(7) != LIVE_LPAR_NL[7]


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["S_named"] is True
    assert out["proved"]["Lpar_R_named_in_p"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["B"]["by_p"]["5"]["named"] == "24"
    assert out["algebra"]["B"]["by_p"]["7"]["named"] == "133"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.430"
