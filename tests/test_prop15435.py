"""Tests for Prop 15.435 — 30/19 and 3/38 are Hoffman mixes."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15433 import AMP_MIX7, AMP_P5, AMP_UNM7
from e1_gmin_m4_prop15434 import LIVE_J
from e1_gmin_m4_prop15435 import (
    FAIL_SUM,
    HOFFMAN_SUM,
    HOFFMAN_W,
    UNMIXED_ABS,
    atom_den19_wrong,
    atom_live,
    main,
    mix_mixed,
    mix_mixed_unsigned,
    mix_unmixed,
    prove_open,
)


def test_unmixed_mix_and_hand19():
    assert sum(HOFFMAN_W) == HOFFMAN_SUM == 19
    assert mix_unmixed(19) == AMP_UNM7 == Fraction(30, 19)
    assert mix_unmixed(FAIL_SUM) == Fraction(5, 3) != AMP_UNM7
    live = atom_live()
    wrong = atom_den19_wrong()
    assert live["p/6"] == Fraction(7, 6)
    assert wrong["p/6"] == Fraction(7, 19) != Fraction(7, 6)
    assert wrong["(p-2)/12"] == Fraction(5, 19) != Fraction(5, 12)
    assert wrong["p/3"] == Fraction(7, 19) != Fraction(7, 3)
    assert wrong["(q+4)/12"] == Fraction(53, 19) != Fraction(53, 12)
    for v in UNMIXED_ABS:
        assert 12 % v.denominator == 0
        assert v.denominator % 19 != 0


def test_mixed_mix_and_unsigned():
    assert mix_mixed(19) == AMP_MIX7 == Fraction(3, 38)
    assert mix_mixed(FAIL_SUM) == Fraction(1, 12) != AMP_MIX7
    assert mix_mixed_unsigned(19) == Fraction(21, 38) != AMP_MIX7


def test_J_delta_is_mix_image():
    assert -4 * AMP_UNM7 == LIVE_J[7]["pp_fp"] == Fraction(-120, 19)
    assert -8 * AMP_MIX7 == LIVE_J[7]["pp_U"] == Fraction(-12, 19)
    assert -4 * mix_unmixed(FAIL_SUM) != LIVE_J[7]["pp_fp"]


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["crrr_joint"] is False
    assert AMP_P5 != Fraction(7, 6)
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["unmixed_hoffman_mix"] is True
    assert out["proved"]["mixed_hoffman_mix"] is True
    assert out["proved"]["J_delta_is_mix_image"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["mix18"] == "5/3"
    assert out["algebra"]["A"]["hand19"]["p/6"] == "7/19"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.435"
