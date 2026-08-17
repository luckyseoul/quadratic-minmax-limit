"""Tests for Prop 15.383 — p≡1 Var(e) floor window; phi_F open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15326 import even_S_p1, Q_lo_named
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15383 import (
    V_e_hi,
    V_e_hi_drop_kappa,
    V_e_lo,
    floor_square_thr,
    floor_square_thr_den1,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_Ve_window():
    A = prove_A()
    assert A["proved"] is True
    assert V_e_hi(5) == Fraction(18, 7)
    assert V_e_lo(5) == Fraction(9, 8)
    assert V_e_lo(5) < Fraction(33, 13) < V_e_hi(5)
    assert V_e_hi_drop_kappa(5) == Fraction(360, 7)
    assert V_e_hi_drop_kappa(5) != V_e_hi(5)
    assert even_S_p1(5, Qpp_floor_ub(5))["sm_joth"] == 6
    assert even_S_p1(5, Q_lo_named(5))["sm_j0"] == 6


def test_square_threshold():
    B = prove_B()
    assert B["proved"] is True
    assert floor_square_thr(5) != floor_square_thr_den1(5)
    assert floor_square_thr(5) == Fraction(3, 4) * 6 * 625


def test_phi_F_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["phi_F_ge_6"] is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Ve_window"] is True
    assert out["proved"]["square_threshold"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
