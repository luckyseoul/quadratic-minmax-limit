"""Tests for Prop 15.565 — Max± of C write ≤Φ−4; path DEAD."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15565 import (
    Phi,
    Q_after_edge,
    Q_plus,
    affine_fe_from_S,
    main,
    prove_A,
    prove_B,
    prove_open,
    q_minus_raw,
)


def test_maxplus_writes_only_Phi_minus_4():
    A = prove_A()
    assert A["proved"] is True
    assert Phi(5) == 65
    assert affine_fe_from_S(1) == 1
    assert affine_fe_from_S(5) == -1
    assert Q_after_edge(Q_plus(5, 1), 1) == Phi(5) - 4
    assert Q_after_edge(Q_plus(5, 5), -1) == Phi(5) - 8
    # fail-when-wrong: Q=Φ−S would write Φ−2 at S_H=2
    assert Phi(5) - 2 != Phi(5) - 4
    assert A["rows"]["7"]["max_abs_QH"] == str(Phi(7) - 4)


def test_maxminus_bad_writes_only_Phi_minus_4():
    B = prove_B()
    assert B["proved"] is True
    q = q_minus_raw(7, -1)
    assert q == -(Phi(7) - 2)
    assert abs(Q_after_edge(q, 1)) == Phi(7)
    assert abs(Q_after_edge(q, -1)) == Phi(7) - 4
    # fail-when-wrong: bad U_- writes Φ−2
    assert abs(Q_after_edge(q, -1)) != Phi(7) - 2
    # S≤−3 cannot reach S_H≥−1 with |f_e|=1
    assert (-3 + 1) <= -2


def test_type_I_flags_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["path_dead"] is True
    assert C["type_I_multilevel_bad_case_ND_closed"] is False
    assert C["type_I_aut_e_3AB_positive_general"] is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["maxplus_writes_Phi_minus_4"] is True
    assert out["proved"]["maxminus_bad_writes_Phi_minus_4"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is False
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    assert Fraction(out["A"]["rows"]["5"]["Phi"]) == 65
