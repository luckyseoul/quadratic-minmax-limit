"""Tests for Prop 15.342 — Q++ / S_k / S^{(k)} named-formula misses."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_bound_4_minus_2_over_pplus2
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15342 import main, prove_A, prove_B, prove_C, prove_open


def test_Sk_not_jacobi_linear():
    A = prove_A()
    assert A["proved"] is True
    assert A["max_err"] > 1.0
    assert A["max_err"] > 1e-6


def test_live_Q_not_pplus2_candidate():
    B = prove_B()
    assert B["proved"] is True
    assert LIVE_QPP[5] != Qpp_bound_4_minus_2_over_pplus2(5)
    assert LIVE_QPP[7] != Qpp_bound_4_minus_2_over_pplus2(7)
    assert LIVE_QPP[5] != Fraction(26, 7)
    assert LIVE_QPP[5] == Fraction(48, 13)
    assert B["by_p"]["5"]["in_interval"] is True
    assert B["by_p"]["7"]["in_interval"] is True


def test_Sk_not_2term_chi():
    C = prove_C()
    assert C["proved"] is True
    assert C["n_2term"] == 0


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Sk_not_ReJ_linear"] is True
    assert out["proved"]["Qpp_not_pplus2_or_endpoint"] is True
    assert out["proved"]["Sk_not_2term_chi"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
