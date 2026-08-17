"""Tests for Prop 15.464 — ns hypergeometric; Paley-N*≠Qn; mean-field misses."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15298 import mu_minus
from e1_gmin_m4_prop15326 import Qn_from_Qpp
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15464 import (
    E_N_ns_drop_den,
    E_N_ns_named,
    live_QN_paley,
    live_Q_by_type,
    main,
    meanfield_Qpp,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_ns_hypergeometric():
    A = prove_A()
    assert A["proved"] is True
    assert E_N_ns_named(5) == mu_minus(5) == Fraction(5, 2)
    assert E_N_ns_named(7) == mu_minus(7) == 7
    assert E_N_ns_named(13) == mu_minus(13)
    assert E_N_ns_drop_den(5) == 10
    assert E_N_ns_drop_den(5) != mu_minus(5)


def test_paley_Nstar_not_Qn():
    B = prove_B()
    assert B["proved"] is True
    assert live_QN_paley(5) == Fraction(32, 13)
    assert live_QN_paley(7) == Fraction(1496, 409)
    assert Qn_from_Qpp(7, LIVE_QPP[7]) == Fraction(1440, 409)
    assert live_QN_paley(7) != Qn_from_Qpp(7, LIVE_QPP[7])
    Q7 = live_Q_by_type(7)
    assert Q7[(1, 1, "sub")] == LIVE_QPP[7] == Fraction(1544, 409)
    assert Q7[(-1, -1, "N1")] == Fraction(1376, 409)
    assert Q7[(-1, -1, "N*")] == Fraction(1496, 409)


def test_meanfield_misses():
    C = prove_C()
    assert C["proved"] is True
    mf5 = meanfield_Qpp(5)
    assert abs(mf5 - float(LIVE_QPP[5])) > 1e-3
    assert Fraction(mf5).limit_denominator(20) != Fraction(48, 13)


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.464"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["E_N_ns_is_mu_minus"] is True
    assert out["proved"]["paley_Nstar_ne_Qn"] is True
    assert out["proved"]["meanfield_misses_Qpp"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
