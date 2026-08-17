"""Tests for Prop 15.429 — C_ns=−(q−1)/4; ns+base ≠ Q_NL."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15429 import (
    C_ns_live,
    C_ns_named,
    C_ns_zero_wrong,
    L_ns_named,
    base_ns_Q,
    main,
    prove_open,
)


def test_C_ns_formula_and_fail():
    assert C_ns_named(5) == Fraction(-6)
    assert C_ns_named(7) == Fraction(-12)
    assert C_ns_named(11) == Fraction(-30)
    assert C_ns_zero_wrong(5) != C_ns_named(5)
    F = _kernel_field(5)
    assert C_ns_live(F, 2, 1) == C_ns_named(5)
    assert C_ns_live(F, 2, -1) == C_ns_named(5)


def test_base_ns_not_QNL():
    from e1_gmin_m4_prop15429 import A0_named

    assert A0_named(5) == 150
    assert A0_named(7) == 588
    assert base_ns_Q(5) == 0 != LIVE_QNL[5]
    assert base_ns_Q(7) == Fraction(-384, 49) != LIVE_QNL[7]
    assert L_ns_named(5) == Fraction(25, 2)
    assert L_ns_named(7) == Fraction(147, 2)


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
    assert out["proved"]["C_ns_named"] is True
    assert out["proved"]["base_ns_not_QNL"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["B"]["p5"] == "0"
    assert out["algebra"]["B"]["p7"] == "-384/49"
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.429"
