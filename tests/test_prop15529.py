"""Tests for Prop 15.529 — Δ_conn is the Ω-triple A–I contraction."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15529 import (
    LIVE,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_AI_hits_live_delta_conn():
    A = prove_A()
    assert A["proved"] is True
    assert Fraction(A["rows"]["5"]["pred_frac"]) == LIVE[5]
    assert Fraction(A["rows"]["7"]["pred_frac"]) == LIVE[7]
    assert LIVE[5] == Fraction(-328, 65)
    assert LIVE[7] == Fraction(-1144, 2863)


def test_drop_q3_is_not_live():
    A = prove_A()
    for p in (5, 7):
        rec = A["rows"][str(p)]
        assert rec["drop_is_live"] is False
        assert Fraction(rec["drop_pred_frac"]) != LIVE[p]
        assert abs(rec["drop_scale"] - rec["q3"]) / rec["q3"] < 1e-6


def test_I_dbl_is_not_two_q_at_p5():
    B = prove_B()
    assert B["proved"] is True
    A = prove_A()
    assert 50 not in A["rows"]["5"]["I_dbl_unique"]
    assert 30 in A["rows"]["5"]["I_dbl_unique"]


def test_type_I_flags_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["type_I_multilevel_bad_case_ND_closed"] is True
    assert C["type_I_aut_e_3AB_positive_general"] is False
    assert C["Delta_conn_p_rational"] is False
    assert type_I_multilevel_bad_case_ND_closed() is True
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["AI_identity"] is True
    assert out["proved"]["I_dbl_not_2p2"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is True
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
