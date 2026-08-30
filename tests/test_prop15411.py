"""Tests for Prop 15.411 — B_NL names p u Q_NL; c=1 misses p≡1 window."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15409 import D_form_on_lattice_general, live_u, u_window
from e1_gmin_m4_prop15411 import (
    B_NL,
    B_NL_drop_1d,
    Q_from_u,
    c1_misses_window,
    main,
    prove_A,
    prove_B,
    prove_open,
    u_c1,
)


def test_B_NL_hits_live_mix():
    A = prove_A()
    assert A["proved"] is True
    assert B_NL(5) == Fraction(128, 3)
    assert B_NL(7) == Fraction(4424, 3)
    assert B_NL_drop_1d(5) == 48
    assert B_NL_drop_1d(5) != B_NL(5)
    assert B_NL(5) == 5 * live_u(5) * Fraction(64, 15)
    assert B_NL(7) == 7 * live_u(7) * Fraction(632, 171)


def test_c1_misses_p1_window():
    B = prove_B()
    assert B["proved"] is True
    assert u_c1(5) == 2
    assert u_window(5) == (2, 2)
    assert u_c1(13) == 6
    assert u_window(13) == (3306, 3854)
    assert c1_misses_window(13)
    assert Q_from_u(13, 6) > Qpp_floor_ub(13)
    assert c1_misses_window(17)
    assert c1_misses_window(29)
    assert c1_misses_window(37)


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False


def test_main():
    out = main()
    assert out["proved"]["B_NL_equals_pu_QNL"] is True
    assert out["proved"]["c1_misses_p1_window"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.411"
