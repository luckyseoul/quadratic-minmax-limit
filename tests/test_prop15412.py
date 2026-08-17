"""Tests for Prop 15.412 — Q_NL is not the 6:4:9 triple or a short rational."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15404 import QT_p3_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import (
    Q_NL_649,
    Q_NL_is_short_rational,
    Q_NL_p7_fit,
    deg1_hits_live,
    main,
    prove_A,
    prove_B,
    prove_open,
    u_from_649,
    weight_649,
)


def test_weight_649_hits_p7_dies_p11():
    A = prove_A()
    assert A["proved"] is True
    assert weight_649(7) == (6, 4, 9)
    assert Q_NL_649(7) == LIVE_QNL[7] == Fraction(632, 171)
    assert u_from_649(7) == 57
    assert u_from_649(11).denominator != 1
    assert u_from_649(11) == Fraction(446339, 649)


def test_Q_NL_not_tempting_or_deg1():
    B = prove_B()
    assert B["proved"] is True
    assert QT_p3_named(7) == Fraction(8, 3)
    assert QT_p3_named(7) != LIVE_QNL[7]
    assert Fraction(4 * 8, 15) != LIVE_QNL[5]
    assert Q_NL_p7_fit(7) == LIVE_QNL[7]
    assert Q_NL_p7_fit(5) != LIVE_QNL[5]
    assert deg1_hits_live() is False
    assert Q_NL_is_short_rational() is False


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False


def test_main():
    out = main()
    assert out["proved"]["weight_649_dies_p11"] is True
    assert out["proved"]["Q_NL_not_short_rational"] is True
    assert out["proved"]["Q_NL_is_short_rational"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.412"
