"""Tests for Prop 15.436 — type L_par den|12 4-points."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15414 import LIVE_LPAR_NL
from e1_gmin_m4_prop15416 import Lpar_R_flag
from e1_gmin_m4_prop15436 import (
    FAIL_SUM,
    Lpar_CRRR,
    Lpar_CRRR_dropC_wrong,
    atom_den19_wrong,
    main,
    mix_Lpar,
    prove_open,
    type_Lpar,
)


def test_type_atoms_and_hand19():
    got = type_Lpar()
    assert got["CRRR"] == Fraction(1337, 12) == Lpar_CRRR()
    assert got["O219"] == Fraction(683, 6)
    assert got["O695"] == Fraction(219, 2)
    for v in got.values():
        assert 12 % v.denominator == 0
        assert v.denominator % 19 != 0
    assert Lpar_CRRR_dropC_wrong() == Lpar_R_flag(7) == 133 != got["CRRR"]
    assert atom_den19_wrong()["CRRR"] == Fraction(1337, 19) != got["CRRR"]


def test_hoffman_mix():
    assert mix_Lpar(19) == LIVE_LPAR_NL[7] == Fraction(2113, 19)
    assert mix_Lpar(FAIL_SUM) == Fraction(2113, 18) != LIVE_LPAR_NL[7]
    assert Lpar_R_flag(7) != LIVE_LPAR_NL[7]


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["L_ns_atoms_are_Lpar"] is False
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert LIVE_QNL[7].denominator % 19 == 0
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["type_Lpar_den12"] is True
    assert out["proved"]["hoffman_mix_is_Lpar_NL"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["B"]["mix18"] == "2113/18"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.436"
