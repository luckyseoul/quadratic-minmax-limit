"""Tests for Prop 15.452 — Type+ S-split; S_NL=2q dies at p=7."""
from __future__ import annotations

from fractions import Fraction
from math import comb

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15452 import (
    E1_S2_named,
    Q_from_ES2,
    Q_from_false_Snl_2q,
    S_star_named,
    S_star_no_quarter,
    main,
    nl_S_values,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    type_plus_S_split,
)


def test_S_star_and_split():
    rec = type_plus_S_split(5)
    assert rec["S"] == S_star_named(5) == 150
    assert rec["nz"] == comb(5, 3) == 10
    assert rec["z"] == n_1d(5) - 10
    assert S_star_no_quarter(5) == 600
    assert S_star_no_quarter(5) != S_star_named(5)
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["7"]["S"] == "588"


def test_E1_inversion_not_live():
    Q1 = Q_from_ES2(5, E1_S2_named(5))
    assert Q1 == 16
    assert Q1 != LIVE_QPP[5]
    assert Q_from_ES2(7, E1_S2_named(7)) == 20
    B = prove_B()
    assert B["proved"] is True


def test_Snl_not_2q_at_p7():
    cnt5 = nl_S_values(5)
    assert set(cnt5) == {Fraction(50)}
    cnt7 = nl_S_values(7)
    assert Fraction(98) in cnt7
    assert set(cnt7) != {Fraction(98)}
    Qf = Q_from_false_Snl_2q(7)
    assert Qf == Fraction(-332, 409)
    assert Qf != LIVE_QPP[7]
    C = prove_C()
    assert C["proved"] is True
    assert C["false_Q7"] == "-332/409"


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.452"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["S_star_named"] is True
    assert out["proved"]["E1_inversion_not_live"] is True
    assert out["proved"]["Snl_2q_dies"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["C"]["false_Q7"] == "-332/409"
