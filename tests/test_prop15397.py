"""Tests for Prop 15.397 — J_++(ψ_2)=−4; λ(ψ_2)=sm_joth."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15326 import Qn_from_Qpp, even_S_p1
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15355 import n_pp_named
from e1_gmin_m4_prop15397 import (
    J_pm1,
    J_type,
    lambda_psi2_drop_m4,
    lambda_psi2_from_Q,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_Jpp_is_minus_four():
    A = prove_A()
    assert A["proved"] is True
    F = _kernel_field(5)
    assert abs(J_type(F, "++", 2) + 4) < 1e-8
    assert abs(J_type(F, "N*", 2) - 2) < 1e-8
    assert abs(J_pm1(F, 2) - 2) < 1e-8
    assert abs(J_type(F, "++", 2) + n_pp_named(5)) > 0.5


def test_lambda_is_sm_joth():
    B = prove_B()
    assert B["proved"] is True
    Qpp = LIVE_QPP[5]
    Qn = Qn_from_Qpp(5, Qpp)
    assert lambda_psi2_from_Q(Qpp, Qn) == even_S_p1(5, Qpp)["sm_joth"]
    assert lambda_psi2_drop_m4(Qpp, Qn) != even_S_p1(5, Qpp)["sm_joth"]


def test_phi_F_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Jpp2_eq_m4"] is True
    assert out["proved"]["lambda_psi2_is_sm_joth"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
