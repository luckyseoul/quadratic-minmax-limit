"""Tests for Prop 15.370 — N4 PGL split; S(x³) not Q; D unnamed."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15368 import dim_Vp
from e1_gmin_m4_prop15370 import (
    N3_LIVE,
    N4_H,
    N4_O,
    cr_type,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_Sx3_not_rational_at_p7():
    A = prove_A()
    assert A["proved"] is True
    assert A["p7_x3_in_Q"] is False
    assert A["p5_S_x3"] == 65
    assert A["p7_x3_palindromic"] is True


def test_D_not_nonlin_named():
    B = prove_B()
    assert B["proved"] is True
    assert B["n_hits_2atom"] == 0
    assert B["n_hits_3atom"] == 0
    assert dim_Vp(7) == 25
    assert dim_Vp(7) != 409


def test_N4_splits_by_cross_ratio():
    C = prove_C()
    assert C["proved"] is True
    assert C["N3"] == N3_LIVE
    assert C["N3_constant"] is True
    assert C["N4_H_pairs"] == 9
    assert C["N4_O_pairs"] == 6
    assert C["fail_N4_13_is_not_5726"] is True
    assert cr_type(1, 2) == "H"
    assert cr_type(1, 3) == "O"
    assert N4_H == 5726
    assert N4_O == 4844
    assert C["pop5_line_unions"] == 210


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["D_named_in_p"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Sx3_not_Q_at_p7"] is True
    assert out["proved"]["D_not_nonlin_named"] is True
    assert out["proved"]["N4_splits_by_PGL_type"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
