"""Tests for Prop 15.306 — 1_D is Paley θ_+; Q++ from E[M²]; E[n⁴] unnamed."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15306 import (
    e_from_theta,
    four_design_moment,
    main,
    paley_e,
    prove_eigenfunction,
    prove_open,
    prove_pairing_2design,
    prove_Qpp_from_M,
    theta_minus,
    theta_plus,
    var_n_named,
)


def test_theta_plus_matches_e():
    A = prove_eigenfunction()
    assert A["proved"] is True
    for p in (5, 7, 11):
        assert e_from_theta(p, theta_plus(p)) == paley_e(p)


def test_theta_minus_mismatches_e():
    for p in (5, 7, 11):
        assert e_from_theta(p, theta_minus(p)) != paley_e(p)


def test_2design_pairing_is_var_n():
    B = prove_pairing_2design()
    assert B["proved"] is True
    for p in (5, 7, 11):
        assert Fraction(B["by_p"][str(p)]["2design_E_pair2"]) == var_n_named(p)


def test_4design_moment_fails_p5():
    # Real 4-design prediction is 1728/168; live central 4th moment is 110/13.
    assert abs(four_design_moment(5) - float(Fraction(1728, 168))) < 1e-9
    assert abs(four_design_moment(5) - float(Fraction(110, 13))) > 0.5


def test_Qpp_from_M_needs_antipodes():
    C = prove_Qpp_from_M()
    assert C["proved"] is True
    assert C["antipode_hit"] >= 1
    assert C["E_M2_named"] is False


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert D["lambda_ge_6"] is False
    assert D["Q_tau_named_in_p"] is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["theta_plus_eigenfunction"] is True
    assert out["proved"]["pairing_2design"] is True
    assert out["proved"]["Qpp_from_EM2"] is True
    assert out["proved"]["En4_named"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
