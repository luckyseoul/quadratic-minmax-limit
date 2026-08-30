"""Tests for Prop 15.415 — L_spl^{1D} named; 1D Paley ++ L named."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15414 import Lpar_1d_named
from e1_gmin_m4_prop15415 import (
    L1d_pp_mix_at_p1_wrong,
    L1d_pp_named,
    Lspl_1d_named,
    main,
    prove_open,
)


def test_Lspl_1d_formula():
    assert Lspl_1d_named(7) == Fraction(539, 5)
    assert Lspl_1d_named(11) == Fraction(4477, 6)
    assert Lspl_1d_named(5) == 25
    assert Lspl_1d_named(5) != Lpar_1d_named(5)


def test_L1d_pp_named():
    assert L1d_pp_named(5) == Lpar_1d_named(5) == Fraction(100, 3)
    assert L1d_pp_named(7) == Fraction(3871, 30)
    assert L1d_pp_named(7) == (2 * Lpar_1d_named(7) + Lspl_1d_named(7)) / 3
    assert L1d_pp_mix_at_p1_wrong(5) != Fraction(100, 3)


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["Lspl_1d_named"] is True
    assert out["proved"]["L1d_pp_named"] is True
    assert out["proved"]["p7_live_Lspl"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["rows"]["7"] == "539/5"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.415"
