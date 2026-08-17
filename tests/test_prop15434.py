"""Tests for Prop 15.434 — J_Δ slot-constant; no formula in p."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15434 import (
    J_delta_live,
    J_zero_wrong,
    LIVE_J,
    extra_ns,
    field_fake_J,
    main,
    meanK_unmixed_fp,
    mean_J,
    prove_open,
)


def test_live_and_fail_zero():
    got5 = J_delta_live(5)
    got7 = J_delta_live(7)
    assert got5 == LIVE_J[5] == {"pp_fp": Fraction(-2), "pp_U": Fraction(2)}
    assert got7 == LIVE_J[7]
    assert got7["pp_fp"] == Fraction(-120, 19)
    assert got7["pp_U"] == Fraction(-12, 19)
    assert got5["pp_fp"] != got5["pp_U"]
    assert J_zero_wrong(5) == 0 != got5["pp_fp"]
    assert extra_ns(5) == Fraction(32, 375) != 0


def test_mean_extra():
    assert mean_J(5) == Fraction(2, 3)
    assert mean_J(7) == Fraction(-84, 19)
    assert extra_ns(5) == Fraction(32, 375)
    assert extra_ns(7) == Fraction(-192, 931)
    assert extra_ns(7).denominator % 19 == 0


def test_K_and_field_fake():
    assert meanK_unmixed_fp(5) == {2}
    assert meanK_unmixed_fp(7) == {2}
    assert meanK_unmixed_fp(11) == {1}
    assert 2 not in meanK_unmixed_fp(11)
    assert field_fake_J(5) == LIVE_J[5]
    assert field_fake_J(7)["pp_fp"] == 2 != LIVE_J[7]["pp_fp"]


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["p7_fp_den19"] is True
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["J_slot_constant"] is True
    assert out["proved"]["mean_extra_named_live"] is True
    assert out["proved"]["K_pattern_dies_p11"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["C"]["S7_has_19"] is False
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.434"
