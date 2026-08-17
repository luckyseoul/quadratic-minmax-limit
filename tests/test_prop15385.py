"""Tests for Prop 15.385 — 1D square-line occupancy named in p."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15385 import (
    En4_1d,
    En4_1d_drop_p3,
    LIVE_N_EQ_P,
    P_1d,
    P_1d_fail_p_as_invp,
    P_ens_p5,
    main,
    n_eq_p_if_nl_misses_full_line,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    w_1d,
)


def test_P_1d_named():
    A = prove_A()
    assert A["proved"] is True
    assert P_1d(5, 0) == Fraction(1, 5)
    assert P_1d(5, 5) == Fraction(2, 15)
    assert P_1d(5, 2) == Fraction(4, 6)
    assert P_1d_fail_p_as_invp(5) != P_1d(5, 5)
    assert sum(P_1d(7, k) for k in range(8)) == 1


def test_En4_1d():
    B = prove_B()
    assert B["proved"] is True
    assert En4_1d(5) == 94
    assert En4_1d(7) == 318
    assert En4_1d_drop_p3(5) != 94


def test_p5_mixture_p7_obstruction():
    C = prove_C()
    assert C["proved"] is True
    assert w_1d(5) == Fraction(3, 13)
    assert P_ens_p5(0) == Fraction(1, 5)
    assert P_ens_p5(5) == Fraction(2, 65)
    assert n_eq_p_if_nl_misses_full_line(7) == 420
    assert LIVE_N_EQ_P[7] == 3948
    assert n_eq_p_if_nl_misses_full_line(7) != LIVE_N_EQ_P[7]


def test_phi_F_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["P_1d_named"] is True
    assert out["proved"]["En4_1d_named"] is True
    assert out["proved"]["p5_mixture_p7_obstruction"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
