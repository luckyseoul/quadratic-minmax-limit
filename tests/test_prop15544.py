"""Tests for Prop 15.544 — p=5 |μ_full|=1/p²; mix G>T."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    T_bitight,
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15543 import mu_1d_named
from e1_gmin_m4_prop15544 import (
    main,
    mix_bound_1d_only_p5,
    mix_bound_p5,
    prove_A,
    prove_B,
    prove_open,
)


def test_mu_full_abs_is_one_over_p2_at_p5():
    A = prove_A()
    assert A["proved"] is True
    live = A["live"]
    assert live["n_kappa1"] == 1800
    assert live["n_abs_1p2"] == 1800
    assert live["n_sign_opp_kappa"] > 0
    assert live["n_pos"] > 0 and live["n_neg"] > 0
    assert live["n_1d"] == 30
    assert live["n_full"] == 100


def test_p5_mix_strictly_inside_T():
    B = prove_B()
    assert B["proved"] is True
    assert mix_bound_p5() == Fraction(3, 65)
    assert mix_bound_p5() < -T_bitight(5)
    assert mix_bound_1d_only_p5() == -T_bitight(5)
    assert mix_bound_1d_only_p5() == abs(mu_1d_named(5, 1))
    # fail-when-wrong: drop full-Ω
    assert mix_bound_p5() != mix_bound_1d_only_p5()
    assert Fraction(B["live_max_mu"]) == Fraction(3, 65)


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
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["mu_full_abs_p5"] is True
    assert out["proved"]["mix_G_gt_T_p5"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is False
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
