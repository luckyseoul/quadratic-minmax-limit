"""Tests for Prop 15.546 — p=7 mix Aut_e G>T."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    T_bitight,
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15546 import (
    C_full_p7,
    C_full_wrong_drop15,
    C_part_p7,
    main,
    mix_bound_p7,
    mix_bound_wrong_drop15,
    prove_A,
    prove_B,
    prove_open,
)


def test_p7_majorants_not_one_over_p2():
    A = prove_A()
    assert A["proved"] is True
    assert C_full_p7() == Fraction(29, 735)
    assert C_part_p7() == Fraction(5, 147)
    assert C_full_p7() != Fraction(1, 49)
    assert C_part_p7() != Fraction(1, 49)
    assert Fraction(A["live"]["max_full"]) == C_full_p7()
    assert Fraction(A["live"]["max_part"]) == C_part_p7()
    assert A["live"]["n_kappa1"] == 14112
    assert A["live"]["n_full_vals"] == 6


def test_p7_mix_below_T_and_drop15_fails():
    B = prove_B()
    assert B["proved"] is True
    assert mix_bound_p7() == Fraction(109, 2863)
    assert mix_bound_p7() < -T_bitight(7)
    assert mix_bound_wrong_drop15() == Fraction(283, 2863)
    assert mix_bound_wrong_drop15() > -T_bitight(7)
    assert C_full_wrong_drop15() == Fraction(29, 245)
    assert Fraction(B["live_max_mu"]) == Fraction(109, 2863)


def test_type_I_flags_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["type_I_multilevel_bad_case_ND_closed"] is False
    assert C["type_I_aut_e_3AB_positive_general"] is False
    assert C["A_full_named_p_law"] is False
    assert C["mu_full_p_law_general"] is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["mu_majorants_p7"] is True
    assert out["proved"]["mix_G_gt_T_p7"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is False
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
