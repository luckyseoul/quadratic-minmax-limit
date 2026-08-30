"""Tests for Prop 15.431 — square-slot C_τ named in p."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15429 import C_ns_named
from e1_gmin_m4_prop15431 import (
    C_Nstar_halfn_wrong,
    C_live,
    C_named,
    C_pm1_allp_wrong,
    main,
    prove_open,
)


def test_named_and_fail():
    assert C_named(5, "pm1", 1) == Fraction(2, 3)
    assert C_named(5, "ns", 1) == C_ns_named(5) == Fraction(-6)
    assert C_named(7, "mm", 1) == Fraction(-20, 3)
    assert C_named(7, "N*", 1) == Fraction(-32, 3)
    assert C_pm1_allp_wrong(5, 1) == 4 != C_named(5, "pm1", 1)
    assert C_named(7, "mm", 1) != 0
    assert C_Nstar_halfn_wrong(7) == -6 != C_named(7, "N*", 1)


def test_live_match_p5p7():
    live5 = C_live(5)
    assert live5[("pm1", 1)] == C_named(5, "pm1", 1)
    assert live5[("pp_U", 1)] == C_named(5, "pp_U", 1)
    assert live5[("N*", -1)] == C_named(5, "N*", -1)
    live7 = C_live(7)
    assert live7[("pp_fp", 1)] == C_named(7, "pp_fp", 1)
    assert live7[("mm", -1)] == C_named(7, "mm", -1)


def test_flags_untouched():
    B = prove_open()
    assert B["proved"] is False
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["C_slots_named"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.431"
