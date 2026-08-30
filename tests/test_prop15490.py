"""Tests for Prop 15.490 — A-lattice numerators over ridge D."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import HPLUS, LIVE_QPP
from e1_gmin_m4_prop15409 import D_1d
from e1_gmin_m4_prop15457 import c_eq_named
from e1_gmin_m4_prop15490 import (
    A_mm_ns,
    A_num,
    LIVE_Q_MM_NS,
    N_from_c,
    Q_mm_ns_named,
    Q_pp_named,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_N_at_c_eq_not_c_min():
    A = prove_A()
    assert A["proved"] is True
    assert N_from_c(5, c_eq_named(5)) == HPLUS[5]
    assert N_from_c(7, c_eq_named(7)) == HPLUS[7]
    assert N_from_c(7, 18) != HPLUS[7]


def test_Qpp_is_8A_over_ridge_D():
    B = prove_B()
    assert B["proved"] is True
    assert Q_pp_named(5, HPLUS[5]) == LIVE_QPP[5]
    assert Q_pp_named(7, HPLUS[7]) == LIVE_QPP[7]
    assert Fraction(8 * A_num(5), D_1d(5)) != LIVE_QPP[5]


def test_Qmm_Ns_shift():
    C = prove_C()
    assert C["proved"] is True
    assert Q_mm_ns_named(5, HPLUS[5]) == LIVE_Q_MM_NS[5]
    assert Q_mm_ns_named(7, HPLUS[7]) == LIVE_Q_MM_NS[7]
    assert A_mm_ns(5) == A_num(5) - 2
    assert A_mm_ns(7) == A_num(7) - 6
    # drop 2(p-4)
    assert Q_pp_named(5, HPLUS[5]) != LIVE_Q_MM_NS[5]


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
